# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Bio.Align support for BED (Browser Extensible Data) files.

The Browser Extensible Data (BED) format stores a series of pairwise
alignments in a single file. Typically, they are used for transcript-to-genome
alignments. BED files store the alignment positions and alignment scores,
but not the aligned sequences.
"""

import sys

import numpy as np

from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Browser Extensible Data (BED) format."""

    def __init__(self, target, bed_n=12):
        """Create an AlignmentWriter object.

        Arguments:
         - target - output stream or file name
         - bed_n - number of columns in the BED file.
                   This must be between 3 and 12; default value is 12.
        """
        if not 3 <= bed_n <= 12:
            raise ValueError("bed_n must be between 3 and 12")

        super().__init__(target)
        self.bed_n = bed_n

    @staticmethod
    def _get_identifier(record, default):
        """Return the record identifier or a default value."""
        try:
            identifier = record.id
        except AttributeError:
            return default

        return default if identifier is None else identifier

    @staticmethod
    def _get_coordinates(alignment):
        """Return normalized alignment coordinates."""
        coordinates = alignment.coordinates

        if not coordinates.size:
            return None

        if coordinates[0, 0] > coordinates[0, -1]:
            return coordinates[:, ::-1]

        return coordinates

    @staticmethod
    def _get_strand(coordinates):
        """Return the BED strand from the alignment coordinates."""
        if coordinates[1, 0] > coordinates[1, -1]:
            return "-"

        return "+"

    @staticmethod
    def _get_blocks(coordinates):
        """Return BED block sizes and block starts."""
        block_sizes = []
        block_starts = []

        target_start, query_start = coordinates[:, 0]

        for target_end, query_end in coordinates[:, 1:].transpose():
            if target_start == target_end:
                query_start = query_end
                continue

            if query_start == query_end:
                target_start = target_end
                continue

            block_size = target_end - target_start
            block_starts.append(target_start)
            block_sizes.append(block_size)

            target_start = target_end
            query_start = query_end

        if not block_starts:
            return [], [], target_start, target_start

        chrom_start = block_starts[0]
        chrom_end = block_starts[-1] + block_sizes[-1]

        return block_sizes, block_starts, chrom_start, chrom_end

    def _build_fields(self, alignment, chrom, strand, blocks):
        """Build BED fields according to the requested BED column count."""
        block_sizes, block_starts, chrom_start, chrom_end = blocks
        fields = [chrom, str(chrom_start), str(chrom_end)]

        if self.bed_n == 3:
            return fields

        query = alignment.sequences[1]
        name = self._get_identifier(query, "query")
        fields.append(name)

        if self.bed_n == 4:
            return fields

        score = getattr(alignment, "score", 0)
        fields.append(format(score, "g"))

        if self.bed_n == 5:
            return fields

        fields.append(strand)

        if self.bed_n == 6:
            return fields

        thick_start = getattr(alignment, "thickStart", chrom_start)
        fields.append(str(thick_start))

        if self.bed_n == 7:
            return fields

        thick_end = getattr(alignment, "thickEnd", chrom_end)
        fields.append(str(thick_end))

        if self.bed_n == 8:
            return fields

        item_rgb = getattr(alignment, "itemRgb", "0")
        fields.append(str(item_rgb))

        if self.bed_n == 9:
            return fields

        fields.append(str(len(block_sizes)))

        if self.bed_n == 10:
            return fields

        fields.append(",".join(map(str, block_sizes)) + ",")

        if self.bed_n == 11:
            return fields

        block_offsets = [
            block_start - chrom_start for block_start in block_starts
        ]
        fields.append(",".join(map(str, block_offsets)) + ",")

        return fields

    def format_alignment(self, alignment):
        """Return a string with one alignment formatted as a BED line."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")

        coordinates = self._get_coordinates(alignment)

        if coordinates is None:
            return ""

        chrom = self._get_identifier(alignment.sequences[0], "target")
        strand = self._get_strand(coordinates)
        blocks = self._get_blocks(coordinates)

        fields = self._build_fields(
            alignment,
            chrom,
            strand,
            blocks,
        )

        return "\t".join(fields) + "\n"


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for Browser Extensible Data (BED) files.

    Each line in the file contains one pairwise alignment, which is loaded
    and returned incrementally. Additional alignment information is stored
    as attributes of each alignment.
    """

    fmt = "BED"

    @staticmethod
    def _parse_blocks(words):
        """Parse BED block sizes and starts."""
        block_count = int(words[9])

        block_sizes = [
            int(size) for size in words[10].rstrip(",").split(",")
        ]
        block_starts = [
            int(start) for start in words[11].rstrip(",").split(",")
        ]

        if len(block_sizes) != block_count:
            raise ValueError(
                "Inconsistent number of block sizes "
                f"({len(block_sizes)} found, expected {block_count})"
            )

        if len(block_starts) != block_count:
            raise ValueError(
                "Inconsistent number of block start positions "
                f"({len(block_starts)} found, expected {block_count})"
            )

        return np.array(block_sizes), np.array(block_starts)

    @staticmethod
    def _build_block_coordinates(block_sizes, block_starts):
        """Build alignment coordinates from BED blocks."""
        target_position = 0
        query_position = 0
        coordinates = [[target_position, query_position]]

        for block_size, block_start in zip(block_sizes, block_starts):
            if block_start != target_position:
                coordinates.append([block_start, query_position])
                target_position = block_start

            target_position += block_size
            query_position += block_size

            coordinates.append([target_position, query_position])

        return np.array(coordinates, np.intp).transpose()

    @staticmethod
    def _get_basic_coordinates(chrom_start, chrom_end):
        """Create coordinates for a BED record without blocks."""
        block_size = chrom_end - chrom_start
        coordinates = np.array(
            [[0, block_size], [0, block_size]],
            np.intp,
        )

        return coordinates, block_size

    @staticmethod
    def _validate_coordinates(
        coordinates,
        chrom_start,
        chrom_end,
    ):
        """Validate BED chromosome coordinates."""
        if chrom_start != coordinates[0, 0]:
            raise ValueError(
                "Inconsistent chromStart found "
                f"({chrom_start}, expected {coordinates[0, 0]})"
            )

        if chrom_end != coordinates[0, -1]:
            raise ValueError(
                "Inconsistent chromEnd found "
                f"({chrom_end}, expected {coordinates[0, -1]})"
            )

    @staticmethod
    def _create_records(chrom, name, query_size):
        """Create query and target sequence records."""
        query_sequence = Seq(None, length=query_size)
        query_record = SeqRecord(
            query_sequence,
            id=name,
            description="",
        )

        target_sequence = Seq(None, length=sys.maxsize)
        target_record = SeqRecord(
            target_sequence,
            id=chrom,
            description="",
        )

        return [target_record, query_record]

    @staticmethod
    def _set_optional_attributes(alignment, words, bed_n):
        """Set optional BED attributes on an alignment."""
        if bed_n <= 4:
            return alignment

        score = words[4]

        try:
            score = float(score)
        except ValueError:
            pass

        alignment.score = score

        if bed_n <= 6:
            return alignment

        alignment.thickStart = int(words[6])

        if bed_n <= 7:
            return alignment

        alignment.thickEnd = int(words[7])

        if bed_n <= 8:
            return alignment

        alignment.itemRgb = words[8]

        return alignment

    def _read_next_alignment(self, stream):
        """Read and return the next alignment from a BED stream."""
        for line in stream:
            words = line.split()
            bed_n = len(words)

            if not 3 <= bed_n <= 12:
                raise ValueError(
                    f"expected between 3 and 12 columns, found {bed_n}"
                )

            chrom = words[0]
            chrom_start = int(words[1])
            chrom_end = int(words[2])
            name = words[3] if bed_n > 3 else None
            strand = words[5] if bed_n > 5 else "+"

            if bed_n > 9:
                block_sizes, block_starts = self._parse_blocks(words)
                coordinates = self._build_block_coordinates(
                    block_sizes,
                    block_starts,
                )
                query_size = sum(block_sizes)
            else:
                coordinates, query_size = self._get_basic_coordinates(
                    chrom_start,
                    chrom_end,
                )

            coordinates[0, :] += chrom_start

            records = self._create_records(
                chrom,
                name,
                query_size,
            )

            if strand == "-":
                coordinates[1, :] = query_size - coordinates[1, :]

            self._validate_coordinates(
                coordinates,
                chrom_start,
                chrom_end,
            )

            alignment = Alignment(records, coordinates)

            return self._set_optional_attributes(
                alignment,
                words,
                bed_n,
            )