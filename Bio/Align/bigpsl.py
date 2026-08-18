# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for alignment files in the bigPsl format.

A bigPsl file is a bigBed file with a BED12+13 format consisting of the 12
predefined BED fields and 13 custom fields defined in the autoSql file
bigPsl.as. This module uses the Bio.Align.bigbed module to parse the file,
but stores the data in a PSL-consistent manner as defined in bigPsl.as. As the
bigPsl format is a special case of the bigBed format, bigPsl files are binary
and are indexed as bigBed files.

See http://genome.ucsc.edu/goldenPath/help/bigPsl.html for more information.

You are expected to use this module via the Bio.Align functions.
"""

import numpy as np

from Bio.Align import Alignment
from Bio.Align import Alignments
from Bio.Align import bigbed
from Bio.Align.bigbed import AutoSQLTable
from Bio.Align.bigbed import Field
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.Seq import UndefinedSequenceError
from Bio.SeqFeature import Location
from Bio.SeqFeature import SeqFeature
from Bio.SeqIO.InsdcIO import _insdc_location_string
from Bio.SeqRecord import SeqRecord

declaration = AutoSQLTable(
    "bigPsl",
    "bigPsl pairwise alignment",
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
            as_type="string",
            name="name",
            comment="Name or ID of item, ideally both human readable and unique",
        ),
        Field(
            as_type="uint",
            name="score",
            comment="Score (0-1000)",
        ),
        Field(
            as_type="char[1]",
            name="strand",
            comment="+ or - indicates whether the query aligns to the + or - strand on the reference",
        ),
        Field(
            as_type="uint",
            name="thickStart",
            comment="Start of where display should be thick (start codon)",
        ),
        Field(
            as_type="uint",
            name="thickEnd",
            comment="End of where display should be thick (stop codon)",
        ),
        Field(
            as_type="uint",
            name="reserved",
            comment="RGB value (use R,G,B string in input file)",
        ),
        Field(
            as_type="int",
            name="blockCount",
            comment="Number of blocks",
        ),
        Field(
            as_type="int[blockCount]",
            name="blockSizes",
            comment="Comma separated list of block sizes",
        ),
        Field(
            as_type="int[blockCount]",
            name="chromStarts",
            comment="Start positions relative to chromStart",
        ),
        Field(
            as_type="uint",
            name="oChromStart",
            comment="Start position in other chromosome",
        ),
        Field(
            as_type="uint",
            name="oChromEnd",
            comment="End position in other chromosome",
        ),
        Field(
            as_type="char[1]",
            name="oStrand",
            comment="+ or -, - means that psl was reversed into BED-compatible coordinates",
        ),
        Field(
            as_type="uint",
            name="oChromSize",
            comment="Size of other chromosome.",
        ),
        Field(
            as_type="int[blockCount]",
            name="oChromStarts",
            comment="Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand",
        ),
        Field(
            as_type="lstring",
            name="oSequence",
            comment="Sequence on other chrom (or edit list, or empty)",
        ),
        Field(
            as_type="string",
            name="oCDS",
            comment="CDS in NCBI format",
        ),
        Field(
            as_type="uint",
            name="chromSize",
            comment="Size of target chromosome",
        ),
        Field(
            as_type="uint",
            name="match",
            comment="Number of bases matched.",
        ),
        Field(
            as_type="uint",
            name="misMatch",
            comment="Number of bases that don't match",
        ),
        Field(
            as_type="uint",
            name="repMatch",
            comment="Number of bases that match but are part of repeats",
        ),
        Field(
            as_type="uint",
            name="nCount",
            comment="Number of 'N' bases",
        ),
        Field(
            as_type="uint",
            name="seqType",
            comment="0=empty, 1=nucleotide, 2=amino_acid",
        ),
    ],
)


class AlignmentWriter(bigbed.AlignmentWriter):
    """Alignment file writer for the bigPsl file format."""

    fmt = "bigPsl"

    def __init__(
        self,
        target,
        targets=None,
        compress=True,
        extraIndex=(),
        cds=False,
        fa=False,
        mask=None,
        wildcard="N",
    ):
        """Create an AlignmentWriter object.

        Arguments:
         - target      - output stream or file name.
         - targets     - A list of SeqRecord objects with the chromosomes in the
                         order as they appear in the alignments. The sequence
                         contents in each SeqRecord may be undefined, but the
                         sequence length must be defined.
         - compress    - If True (default), compress data using zlib.
         - extraIndex  - List of strings with the names of extra columns to be
                         indexed.
         - cds         - If True, look for a query feature of type CDS.
         - fa          - If True, include the query sequence in the PSL file.
         - mask        - Specify if repeat regions in the target sequence are
                         masked.
         - wildcard    - Report alignments to the wildcard character in the
                         nCount field.
        """
        super().__init__(
            target,
            bedN=12,
            declaration=declaration,
            targets=targets,
            compress=compress,
            extraIndex=extraIndex,
        )
        self.cds = cds
        self.fa = fa
        self.mask = mask
        self.wildcard = wildcard

    def _get_alignment_sequences(self, alignment):
        """Get target and query sequences from an alignment."""
        target, query = alignment.sequences

        try:
            query = query.seq
        except AttributeError:
            pass

        try:
            target = target.seq
        except AttributeError:
            pass

        return target, query

    def _normalize_strand(self, coordinates, target, query):
        """Normalize sequences and coordinates according to strand."""
        t_size = len(target)
        q_size = len(query)
        dnax = None

        if coordinates[1, 0] > coordinates[1, -1]:
            strand = "-"
            query = reverse_complement(query)
            coordinates = coordinates.copy()
            coordinates[1, :] = q_size - coordinates[1, :]

        elif coordinates[0, 0] > coordinates[0, -1]:
            strand = "-"
            target = reverse_complement(target)
            coordinates = coordinates.copy()
            coordinates[0, :] = t_size - coordinates[0, :]
            dnax = True

        else:
            strand = "+"

        return (
            coordinates,
            target,
            query,
            t_size,
            q_size,
            strand,
            dnax,
        )

    def _sequence_bytes(self, sequence):
        """Convert a sequence to bytes when possible."""
        try:
            return bytes(sequence)
        except TypeError:
            return bytes(sequence, "ASCII")
        except UndefinedSequenceError:
            return None

    def _count_sequence_matches(
        self,
        t_seq,
        q_seq,
        wildcard,
        mask,
        q_count,
    ):
        """Count matches, mismatches, repeats and wildcard bases."""
        matches = 0
        mis_matches = 0
        rep_matches = 0
        n_count = 0

        if t_seq is None or q_seq is None:
            return q_count, mis_matches, rep_matches, n_count

        if mask == "lower":
            values = zip(
                t_seq.upper(),
                q_seq.upper(),
                t_seq,
            )

            for u1, u2, c1 in values:
                if u1 == wildcard or u2 == wildcard:
                    n_count += 1
                elif u1 == u2:
                    if u1 == c1:
                        matches += 1
                    else:
                        rep_matches += 1
                else:
                    mis_matches += 1

        elif mask == "upper":
            values = zip(
                t_seq.lower(),
                q_seq.lower(),
                t_seq,
            )

            for u1, u2, c1 in values:
                if u1 == wildcard or u2 == wildcard:
                    n_count += 1
                elif u1 == u2:
                    if u1 == c1:
                        matches += 1
                    else:
                        rep_matches += 1
                else:
                    mis_matches += 1

        else:
            for u1, u2 in zip(t_seq.upper(), q_seq.upper()):
                if u1 == wildcard or u2 == wildcard:
                    n_count += 1
                elif u1 == u2:
                    matches += 1
                else:
                    mis_matches += 1

        return matches, mis_matches, rep_matches, n_count

    def _calculate_blocks(
        self,
        alignment,
        target,
        query,
        coordinates,
        dnax,
    ):
        """Calculate alignment blocks and sequence statistics."""
        wildcard = self.wildcard
        mask = self.mask

        matches = 0
        mis_matches = 0
        rep_matches = 0
        n_count = 0

        block_sizes = []
        q_starts = []
        t_starts = []

        t_start, q_start = coordinates[:, 0]

        for t_end, q_end in coordinates[:, 1:].transpose():
            if t_start == t_end:
                q_start = q_end

            elif q_start == q_end:
                t_start = t_end

            else:
                t_count = t_end - t_start
                q_count = q_end - q_start

                t_starts.append(t_start)
                q_starts.append(q_start)
                block_sizes.append(q_count)

                if t_count == q_count:
                    assert dnax is not True
                    dnax = False
                else:
                    assert t_count == 3 * q_count
                    assert dnax is not False
                    dnax = True

                t_seq = target[t_start:t_end]
                q_seq = query[q_start:q_end]

                t_seq = self._sequence_bytes(t_seq)
                q_seq = self._sequence_bytes(q_seq)

                (
                    block_matches,
                    block_mis_matches,
                    block_rep_matches,
                    block_n_count,
                ) = self._count_sequence_matches(
                    t_seq,
                    q_seq,
                    wildcard,
                    mask,
                    q_count,
                )

                matches += block_matches
                mis_matches += block_mis_matches
                rep_matches += block_rep_matches
                n_count += block_n_count

                t_start = t_end
                q_start = q_end

        return (
            np.array(t_starts),
            np.array(q_starts),
            np.array(block_sizes),
            matches,
            mis_matches,
            rep_matches,
            n_count,
            dnax,
        )

    def _get_alignment_statistics(self, alignment):
        """Get statistics already stored in the alignment."""
        statistics = {
            "matches": "matches",
            "mis_matches": "misMatches",
            "rep_matches": "repMatches",
            "n_count": "nCount",
        }

        values = {}

        for key, attribute in statistics.items():
            values[key] = getattr(
                alignment,
                attribute,
                None,
            )

        return values

    def _get_cds(self, alignment):
        """Get the CDS annotation when requested."""
        if not self.cds:
            return ""

        for feature in alignment.query.features:
            if feature.type == "CDS":
                return _insdc_location_string(
                    feature.location,
                    len(alignment.query),
                )

        return "n/a"

    def _get_sequence_type(self, alignment):
        """Determine the sequence type."""
        molecule_type = alignment.query.annotations.get(
            "molecule_type"
        )

        if molecule_type == "DNA":
            return "1"

        if molecule_type == "protein":
            return "2"

        return "0"

    def _set_alignment_annotations(
        self,
        alignment,
        t_size,
        q_size,
        q_start,
        q_end,
        o_strand,
        q_starts,
        o_sequence,
        o_cds,
        matches,
        mis_matches,
        rep_matches,
        n_count,
    ):
        """Set the bigPsl annotations."""
        alignment.annotations["oChromStart"] = str(q_start)
        alignment.annotations["oChromEnd"] = str(q_end)
        alignment.annotations["oStrand"] = o_strand
        alignment.annotations["oChromSize"] = str(q_size)
        alignment.annotations["oChromStarts"] = ",".join(
            map(str, q_starts)
        )
        alignment.annotations["oSequence"] = o_sequence
        alignment.annotations["oCDS"] = o_cds
        alignment.annotations["chromSize"] = str(t_size)
        alignment.annotations["match"] = str(matches)
        alignment.annotations["misMatch"] = str(mis_matches)
        alignment.annotations["repMatch"] = str(rep_matches)
        alignment.annotations["nCount"] = str(n_count)
        alignment.annotations["seqType"] = self._get_sequence_type(
            alignment
        )

    def _prepare_alignment(self, alignment):
        """Prepare a single alignment for writing."""
        coordinates = alignment.coordinates

        if not coordinates.size:
            return None

        target, query = self._get_alignment_sequences(
            alignment
        )

        (
            coordinates,
            target,
            query,
            t_size,
            q_size,
            strand,
            dnax,
        ) = self._normalize_strand(
            coordinates,
            target,
            query,
        )

        (
            t_starts,
            q_starts,
            block_sizes,
            matches,
            mis_matches,
            rep_matches,
            n_count,
            dnax,
        ) = self._calculate_blocks(
            alignment,
            target,
            query,
            coordinates,
            dnax,
        )

        q_start = q_starts[0]
        q_end = q_starts[-1] + block_sizes[-1]

        o_strand = "+"

        if strand == "-":
            if dnax is True:
                o_strand = "-"
                q_starts = q_size - (
                    q_starts + block_sizes
                )
                q_starts = q_starts[::-1]
                alignment.coordinates = (
                    alignment.coordinates[:, ::-1]
                )
            else:
                q_start = q_size - q_end
                q_end = q_size - q_starts[0]

        o_sequence = ""

        if self.fa:
            o_sequence = str(alignment.query.seq)

        o_cds = self._get_cds(alignment)

        statistics = self._get_alignment_statistics(
            alignment
        )

        if statistics["matches"] is not None:
            matches = statistics["matches"]

        if statistics["mis_matches"] is not None:
            mis_matches = statistics["mis_matches"]

        if statistics["rep_matches"] is not None:
            rep_matches = statistics["rep_matches"]

        if statistics["n_count"] is not None:
            n_count = statistics["n_count"]

        self._set_alignment_annotations(
            alignment,
            t_size,
            q_size,
            q_start,
            q_end,
            o_strand,
            q_starts,
            o_sequence,
            o_cds,
            matches,
            mis_matches,
            rep_matches,
            n_count,
        )

        return alignment

    def write_file(self, stream, alignments):
        """Write the file."""
        fixed_alignments = Alignments()

        for alignment in alignments:
            if not isinstance(alignment, Alignment):
                raise TypeError(
                    "Expected an Alignment object"
                )

            alignment = self._prepare_alignment(
                alignment
            )

            if alignment is not None:
                fixed_alignments.append(alignment)

        fixed_alignments.sort(
            key=lambda alignment: (
                alignment.target.id,
                alignment.coordinates[0, 0],
            )
        )

        fixed_alignments.targets = alignments.targets

        bigbed.AlignmentWriter(
            stream,
            bedN=12,
            declaration=declaration,
            compress=self.compress,
        ).write(fixed_alignments)

class AlignmentIterator(bigbed.AlignmentIterator):
    """Alignment iterator for bigPsl files.

    The pairwise alignments stored in the bigPsl file are loaded and returned
    incrementally.  Additional alignment information is stored as attributes
    of each alignment.
    """

    fmt = "bigPsl"

    def _analyze_fields(self, fields, fieldCount, definedFieldCount):
        names = (
            "chrom",
            "chromStart",
            "chromEnd",
            "name",  # 0
            "score",  # 1
            "strand",  # 2
            "thickStart",  # 3
            "thickEnd",  # 4
            "reserved",  # 5
            "blockCount",  # 6
            "blockSizes",  # 7
            "chromStarts",  # 8
            "oChromStart",  # 9
            "oChromEnd",  # 10
            "oStrand",  # 11
            "oChromSize",  # 12
            "oChromStarts",  # 13
            "oSequence",  # 14
            "oCDS",  # 15
            "chromSize",  # 16
            "match",  # 17
            "misMatch",  # 18
            "repMatch",  # 19
            "nCount",  # 20
            "seqType",  # 21
        )
        for i, name in enumerate(names):
            if name != fields[i].name:
                raise ValueError(
                    f"Expected field name '{name}'; found '{fields[i].name}'"
                )

    def _create_alignment(self, chromId, tStart, tEnd, rest, dataStart, dataEnd):
        assert rest[dataEnd - 1] == 0
        words = rest[dataStart : dataEnd - 1].decode().split("\t")
        if len(words) != 22:
            raise ValueError(
                "Unexpected number of fields (%d, expected 22)" % len(words)
            )
        target_record = self.targets[chromId]
        tSize = int(words[16])
        if len(target_record) != tSize:
            raise ValueError(
                "Unexpected chromosome size %d (expected %d)"
                % (tSize, len(target_record))
            )
        strand = words[2]
        qName = words[0]
        qSize = int(words[12])
        blockCount = int(words[6])
        blockSizes = [int(blockSize) for blockSize in words[7].rstrip(",").split(",")]
        tStarts = [int(start) for start in words[8].rstrip(",").split(",")]
        qStarts = [int(start) for start in words[13].rstrip(",").split(",")]
        if len(blockSizes) != blockCount:
            raise ValueError(
                "Inconsistent number of blocks (%d found, expected %d)"
                % (len(blockSizes), blockCount)
            )
        if len(qStarts) != blockCount:
            raise ValueError(
                "Inconsistent number of query start positions (%d found, expected %d)"
                % (len(qStarts), blockCount)
            )
        if len(tStarts) != blockCount:
            raise ValueError(
                "Inconsistent number of target start positions (%d found, expected %d)"
                % (len(qStarts), blockCount)
            )
        qStarts = np.array(qStarts)
        tStarts = np.array(tStarts)
        tBlockSizes = np.array(blockSizes)
        query_sequence = words[14]
        if query_sequence == "":
            query_sequence = Seq(None, length=qSize)
        else:
            query_sequence = Seq(query_sequence)
            if len(query_sequence) != qSize:
                raise ValueError(
                    "Inconsistent query sequence length (%d, expected %d)"
                    % (len(query_sequence), qSize)
                )
        query_record = SeqRecord(query_sequence, id=qName)
        cds = words[15]
        if cds and cds != "n/a":
            location = Location.fromstring(cds)
            feature = SeqFeature(location, type="CDS")
            query_record.features.append(feature)
        seqType = words[21]
        if seqType == "0":
            qBlockSizes = tBlockSizes
        elif seqType == "1":
            query_record.annotations["molecule_type"] = "DNA"
            qBlockSizes = tBlockSizes
        elif seqType == "2":
            query_record.annotations["molecule_type"] = "protein"
            qBlockSizes = tBlockSizes // 3
        else:
            raise ValueError("Unexpected sequence type '%s'" % seqType)
        tStarts += tStart
        qStrand = words[11]
        if qStrand == "-" and strand == "-":
            tStart, tEnd = tEnd, tStart
            qStarts = qSize - qStarts - qBlockSizes
            tStarts = tSize - tStarts - tBlockSizes
            qStarts = qStarts[::-1]
            tStarts = tStarts[::-1]
            qBlockSizes = qBlockSizes[::-1]
            tBlockSizes = tBlockSizes[::-1]
        qPosition = qStarts[0]
        tPosition = tStarts[0]
        coordinates = [[tPosition, qPosition]]
        for tB, qB, tS, qS in zip(tBlockSizes, qBlockSizes, tStarts, qStarts):
            if tS != tPosition:
                coordinates.append([tS, qPosition])
                tPosition = tS
            if qS != qPosition:
                coordinates.append([tPosition, qS])
                qPosition = qS
            tPosition += tB
            qPosition += qB
            coordinates.append([tPosition, qPosition])
        coordinates = np.array(coordinates, np.intp).transpose()
        qStart = int(words[9])
        qEnd = int(words[10])
        if strand == "-":
            if qStrand == "-":
                coordinates[0, :] = tSize - coordinates[0, :]
            else:
                qStart, qEnd = qEnd, qStart
                coordinates[1, :] = qSize - coordinates[1, :]
        if tStart != coordinates[0, 0]:
            raise ValueError(
                "Inconsistent tStart found (%d, expected %d)"
                % (tStart, coordinates[0, 0])
            )
        if tEnd != coordinates[0, -1]:
            raise ValueError(
                "Inconsistent tEnd found (%d, expected %d)" % (tEnd, coordinates[0, -1])
            )
        if qStart != coordinates[1, 0]:
            raise ValueError(
                "Inconsistent qStart found (%d, expected %d)"
                % (qStart, coordinates[1, 0])
            )
        if qEnd != coordinates[1, -1]:
            raise ValueError(
                "Inconsistent qEnd found (%d, expected %d)" % (qEnd, coordinates[1, -1])
            )
        records = [target_record, query_record]
        alignment = Alignment(records, coordinates)
        alignment.annotations = {}
        score = words[1]
        try:
            score = float(score)
        except ValueError:
            pass
        else:
            if score.is_integer():
                score = int(score)
        alignment.score = score
        alignment.thickStart = int(words[3])
        alignment.thickEnd = int(words[4])
        alignment.itemRgb = words[5]
        alignment.matches = int(words[17])
        alignment.misMatches = int(words[18])
        alignment.repMatches = int(words[19])
        alignment.nCount = int(words[20])
        return alignment