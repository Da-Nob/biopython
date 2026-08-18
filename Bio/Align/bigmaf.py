# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "bigmaf" multiple alignment format.

The bigMaf format stores multiple alignments in a format compatible with
the MAF (Multiple Alignment Format) format. BigMaf files are binary and are
indexed as a bigBed file.

See https://genome.ucsc.edu/goldenPath/help/bigMaf.html
"""

import re
import struct
import zlib

import numpy as np

from Bio.Align import _aligncore  # type: ignore
from Bio.Align import Alignment
from Bio.Align import Alignments
from Bio.Align import bigbed
from Bio.Align import maf
from Bio.Align.bigbed import AutoSQLTable
from Bio.Align.bigbed import Field
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


declaration = AutoSQLTable(
    "bedMaf",
    "Bed3 with MAF block",
    [
        Field(
            as_type="string",
            name="chrom",
            comment="Reference sequence chromosome or scaffold",
        ),
        Field(
            as_type="uint",
            name="chromStart",
            comment="Start position in chromosome",
        ),
        Field(
            as_type="uint",
            name="chromEnd",
            comment="End position in chromosome",
        ),
        Field(
            as_type="lstring",
            name="mafBlock",
            comment="MAF block",
        ),
    ],
)


class AlignmentWriter(bigbed.AlignmentWriter):
    """Alignment file writer for the bigMaf file format."""

    fmt = "bigMaf"

    def __init__(
        self,
        target,
        targets=None,
        compress=True,
        blockSize=256,
        itemsPerSlot=512,
    ):
        """Create an AlignmentWriter object.

        Arguments:
         - target       - output stream or file name.
         - targets      - A list of SeqRecord objects with the chromosomes in
                          the order as they appear in the alignments. The
                          sequence contents in each SeqRecord may be undefined,
                          but the sequence length must be defined, as in this
                          example:

                          SeqRecord(Seq(None, length=248956422), id="chr1")

                          If targets is None (the default value), the alignments
                          must have an attribute .targets providing the list
                          of SeqRecord objects.
         - compress     - If True (default), compress data using zlib.
                          If False, do not compress data.
                          Use compress=False for faster searching.
         - blockSize    - Number of items to bundle in r-tree.
                          See UCSC's bedToBigBed program for more information.
                          Default value is 256.
         - itemsPerSlot - Number of data points bundled at lowest level.
                          See UCSC's bedToBigBed program for more information.
                          Use itemsPerSlot=1 for faster searching.
                          Default value is 512.
        """
        super().__init__(
            target,
            bedN=3,
            declaration=declaration,
            targets=targets,
            compress=compress,
            blockSize=blockSize,
            itemsPerSlot=itemsPerSlot,
        )

    def write_file(self, stream, alignments):
        """Write the file."""
        fixed_alignments = Alignments()

        for alignment in alignments:
            if not isinstance(alignment, Alignment):
                raise TypeError("Expected an Alignment object")

            maf_block = format(alignment, "maf")[:-1].replace("\n", ";")
            coordinates = alignment.coordinates

            if not coordinates.size:
                continue

            alignment = alignment[:2]
            reference, chromosome = alignment.target.id.split(".", 1)
            alignment.target.id = chromosome

            assert coordinates[0, 0] < coordinates[0, -1]

            alignment.annotations = {}
            alignment.annotations["mafBlock"] = maf_block
            fixed_alignments.append(alignment)

        fixed_alignments.sort(
            key=lambda alignment: (
                alignment.target.id,
                alignment.coordinates[0, 0],
            )
        )

        record = alignments.targets[0]
        reference, chromosome = record.id.split(".", 1)

        targets = list(alignments.targets)
        targets[0] = SeqRecord(record.seq, id=chromosome)

        fixed_alignments.targets = targets

        bigbed.AlignmentWriter(
            stream,
            bedN=3,
            declaration=declaration,
            compress=self.compress,
        ).write(fixed_alignments)


class AlignmentIterator(bigbed.AlignmentIterator, maf.AlignmentIterator):
    """Alignment iterator for bigMaf files.

    The file may contain multiple alignments, which are loaded and returned
    incrementally.

    Alignment annotations are stored in the ``.annotations`` attribute of the
    ``Alignment`` object, except for the alignment score, which is stored as an
    attribute. Sequence information of empty parts in the alignment block
    (sequences that connect the previous alignment block to the next alignment
    block, but do not align to the current alignment block) is stored in the
    alignment annotations under the ``"empty"`` key. Annotations specific to
    each line in the alignment are stored in the ``.annotations`` attribute of
    the corresponding sequence record.
    """

    fmt = "bigMaf"
    mode = "b"

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
        - source - input file stream, or path to input file
        """
        self.reference = None
        super().__init__(source)

    def _read_reference(self, stream):
        # Supplemental Table 12: Binary BED-data format
        # chromId     4 bytes, unsigned
        # chromStart  4 bytes, unsigned
        # chromEnd    4 bytes, unsigned
        # rest        zero-terminated string in tab-separated format
        formatter = struct.Struct(self.byteorder + "III")
        size = formatter.size
        node = self.tree

        while True:
            try:
                children = node.children
            except AttributeError:
                break
            else:
                node = children[0]

        filepos = stream.tell()
        stream.seek(node.dataOffset)

        data_size = 256
        data = b""
        compressed_data = b""

        while True:
            chunk = stream.read(data_size)

            if self._compressed:
                compressed_data += chunk
                decompressor = zlib.decompressobj()
                data = decompressor.decompress(compressed_data)
            else:
                data += chunk

            try:
                i = data.index(b";s", size)
            except ValueError:
                continue

            words = data[i + 1 :].split()

            if len(words) > 2:
                break

        name = words[1]

        stream.seek(filepos)

        reference, chromosome = name.split(b".", 1)

        return reference.decode()

    def _read_header(self, stream):
        super()._read_header(stream)

        if self.reference is None:
            self.reference = self._read_reference(stream)
            self._index = 0

        self.targets[0].id = f"{self.reference}.{self.targets[0].id}"

    def _get_line_end(self, buffer, start, pattern):
        """Return the end position of a line matching the given pattern."""
        match = re.match(pattern, buffer[start:])

        if match is None:
            raise ValueError(
                "Error parsing alignment - invalid line format"
            )

        return start + match.span()[1]

    def _parse_annotation_line(
        self,
        buffer,
        start,
        annotations,
    ):
        """Parse an MAF annotation ('a') line."""
        end = self._get_line_end(
            buffer,
            start,
            b"^[^;]*",
        )

        line = buffer[start:end].tobytes()
        score = None

        for word in line[1:].split():
            key, value = word.split(b"=")

            if key == b"score":
                score = float(value)
            elif key == b"pass":
                value = int(value)

                if value <= 0:
                    raise ValueError(
                        "pass value must be positive (found %d)" % value
                    )

                annotations["pass"] = value
            else:
                raise ValueError(
                    "Unknown annotation variable '%s'" % key.decode()
                )

        return end, score

    def _parse_sequence_line(
        self,
        buffer,
        start,
        parser,
        data,
        data_start,
        records,
        strands,
    ):
        """Parse an MAF sequence ('s') line."""
        end = self._get_line_end(
            buffer,
            start,
            rb"^s\s*\S*\s*\d*\s*\d*\s*[+-]\s*\d*\s*",
        )

        line = buffer[start:end].tobytes()
        words = line.split(None, 5)

        if len(words) != 6:
            raise ValueError(
                "Error parsing alignment - 's' line must have 7 fields"
            )

        src = words[1].decode()
        start_position = int(words[2])
        size = int(words[3])
        strand = words[4]
        src_size = int(words[5])

        parser_position = end

        n, sequence = parser.feed(
            data,
            data_start + parser_position,
        )

        if len(sequence) != size:
            raise ValueError(
                "sequence size is incorrect (found %d, expected %d)"
                % (len(sequence), size)
            )

        seq = Seq(
            {start_position: sequence},
            length=src_size,
        )

        record = SeqRecord(
            seq,
            id=src,
            name="",
            description="",
        )

        records.append(record)
        strands.append(strand)

        return parser_position + n, src, record

    def _parse_insert_line(
        self,
        buffer,
        start,
        src,
        record,
    ):
        """Parse an MAF insertion ('i') line."""
        end = self._get_line_end(
            buffer,
            start,
            b"^[^;]*",
        )

        line = buffer[start:end].tobytes()
        words = line.split(None, 5)

        assert len(words) == 6
        assert words[1].decode() == src

        left_status = words[2].decode()
        left_count = int(words[3])
        right_status = words[4].decode()
        right_count = int(words[5])

        assert left_status in AlignmentIterator.status_characters
        assert right_status in AlignmentIterator.status_characters

        record.annotations["leftStatus"] = left_status
        record.annotations["leftCount"] = left_count
        record.annotations["rightStatus"] = right_status
        record.annotations["rightCount"] = right_count

        return end

    def _parse_empty_line(
        self,
        buffer,
        start,
        annotations,
    ):
        """Parse an MAF empty ('e') line."""
        end = self._get_line_end(
            buffer,
            start,
            b"^[^;]*",
        )

        line = buffer[start:end].tobytes()
        words = line.split(None, 6)

        assert len(words) == 7

        src = words[1].decode()
        start_position = int(words[2])
        size = int(words[3])
        strand = words[4]
        src_size = int(words[5])
        status = words[6].decode()

        assert status in AlignmentIterator.empty_status_characters

        sequence = Seq(None, length=src_size)

        record = SeqRecord(
            sequence,
            id=src,
            name="",
            description="",
        )

        end_position = start_position + size

        if strand == b"+":
            segment = (start_position, end_position)
        else:
            segment = (
                src_size - start_position,
                src_size - end_position,
            )

        empty = (record, segment, status)

        annotation = annotations.get("empty")

        if annotation is None:
            annotation = []
            annotations["empty"] = annotation

        annotation.append(empty)

        return end

    def _parse_quality_line(
        self,
        buffer,
        start,
        src,
        record,
    ):
        """Parse an MAF quality ('q') line."""
        end = self._get_line_end(
            buffer,
            start,
            b"^[^;]*",
        )

        line = buffer[start:end].tobytes()
        words = line.split(None, 2)

        assert len(words) == 3
        assert words[1].decode() == src

        value = words[2].replace(b"-", b"")

        record.annotations["quality"] = value.decode()

        return end

    def _build_coordinates(
        self,
        parser,
        records,
        strands,
    ):
        """Build alignment coordinates from parsed MAF records."""
        shape = parser.shape
        coordinates = np.empty(shape, np.intp)
        parser.fill(coordinates)

        for record, strand, row in zip(
            records,
            strands,
            coordinates,
        ):
            if strand == b"+":
                row += record.seq.defined_ranges[0][0]
            else:
                record.seq = record.seq.reverse_complement()
                row[:] = record.seq.defined_ranges[0][1] - row

        return coordinates

    def _create_alignment(
        self,
        chromId,
        chromStart,
        chromEnd,
        data,
        dataStart,
        dataEnd,
    ):
        assert data[dataEnd - 1] == 0

        buffer = memoryview(data)[dataStart:dataEnd]

        records = []
        strands = []
        annotations = {}
        score = None

        printed_alignment_parser = _aligncore.PrintedAlignmentParser(b";")

        j = -1
        src = None
        record = None

        while True:
            i = j + 1
            prefix = buffer[i : i + 1]

            if prefix == b"#":
                j = self._get_line_end(
                    buffer,
                    i,
                    b"^[^;]*",
                )

            elif prefix == b"a":
                j, score = self._parse_annotation_line(
                    buffer,
                    i,
                    annotations,
                )

            elif prefix == b"s":
                j, src, record = self._parse_sequence_line(
                    buffer,
                    i,
                    printed_alignment_parser,
                    data,
                    dataStart,
                    records,
                    strands,
                )

            elif prefix == b"i":
                j = self._parse_insert_line(
                    buffer,
                    i,
                    src,
                    record,
                )

            elif prefix == b"e":
                j = self._parse_empty_line(
                    buffer,
                    i,
                    annotations,
                )

            elif prefix == b"q":
                j = self._parse_quality_line(
                    buffer,
                    i,
                    src,
                    record,
                )

            elif prefix == b"\00":
                break

            else:
                raise ValueError(
                    f"Error parsing alignment - unexpected line:\n"
                    f"{buffer[i:j].tobytes()}"
                )

        coordinates = self._build_coordinates(
            printed_alignment_parser,
            records,
            strands,
        )

        alignment = Alignment(
            records,
            coordinates,
        )

        if annotations:
            alignment.annotations = annotations

        if score is not None:
            alignment.score = score

        return alignment