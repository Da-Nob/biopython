# Copyright 2012 by Eric Talevich.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for accessing sequences in PDB and mmCIF files."""

import collections
import warnings

from Bio import BiopythonParserWarning
from Bio.Data.PDBData import protein_letters_3to1
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import _TextIOSource
from .Interfaces import SequenceIterator

_aa3to1_dict = {}
_aa3to1_dict.update(protein_letters_3to1)
_aa3to1_dict.update(protein_letters_3to1_extended)


def _res2aacode(residue, undef_code="X"):
    """Return the one-letter amino acid code from the residue name.

    Non-amino acid are returned as "X".
    """
    if isinstance(residue, str):
        return _aa3to1_dict.get(residue, undef_code)

    return _aa3to1_dict.get(residue.resname, undef_code)


def _residues_for_chain(chain):
    """Return the chain's residues that map to a known amino acid (PRIVATE).

    HETATM mod. res. policy: keep a modified residue if it is recognised,
    else discard it.
    """
    return [
        res
        for res in chain.get_unpacked_list()
        if _res2aacode(res.get_resname().upper()) != "X"
    ]


def _find_gaps(rnumbers):
    """Return the list of gaps in a chain's residue numbering (PRIVATE).

    Each gap is a tuple ``(index, pre_gap_number, post_gap_number)``, found
    wherever two consecutive residue numbers are neither equal nor exactly
    one apart.
    """
    gaps = []
    for i, rnum in enumerate(rnumbers[:-1]):
        if rnumbers[i + 1] != rnum + 1 and rnumbers[i + 1] != rnum:
            gaps.append((i + 1, rnum, rnumbers[i + 1]))
    return gaps


def _fill_gaps(residues, gaps):
    """Return one-letter codes for ``residues``, padding gaps with 'X' (PRIVATE).

    Where amino acids are missing from the structure, as indicated by
    ``gaps``, the sequence is filled in with 'X' characters to match the
    size of the missing region.
    """
    if not gaps:
        return [_res2aacode(x) for x in residues]

    res_out = []
    prev_idx = 0
    for i, pregap, postgap in gaps:
        if postgap <= pregap:
            warnings.warn(
                "Ignoring out-of-order residues after a gap",
                BiopythonParserWarning,
            )
            # Keep the normal part, drop the out-of-order segment
            # (presumably modified or hetatm residues, e.g. 3BEG)
            res_out.extend(_res2aacode(x) for x in residues[prev_idx:i])
            return res_out
        gapsize = postgap - pregap - 1
        res_out.extend(_res2aacode(x) for x in residues[prev_idx:i])
        prev_idx = i
        res_out.append("X" * gapsize)

    # Last segment
    res_out.extend(_res2aacode(x) for x in residues[prev_idx:])
    return res_out


def AtomIterator(pdb_id, structure):
    """Return SeqRecords from Structure objects.

    Base function for sequence parsers that read structures Bio.PDB parsers.

    Once a parser from Bio.PDB has been used to load a structure into a
    Bio.PDB.Structure.Structure object, there is no difference in how the
    sequence parser interprets the residue sequence. The functions in this
    module may be used by SeqIO modules wishing to parse sequences from lists
    of residues.

    Calling functions must pass a Bio.PDB.Structure.Structure object.


    See Bio.SeqIO.PdbIO.PdbAtomIterator and Bio.SeqIO.PdbIO.CifAtomIterator for
    details.
    """
    model = structure[0]
    for chn_id, chain in sorted(model.child_dict.items()):
        residues = _residues_for_chain(chain)
        if not residues:
            continue

        # Identify missing residues in the structure
        # (fill the sequence with 'X' residues in these regions)
        rnumbers = [r.id[1] for r in residues]
        gaps = _find_gaps(rnumbers)
        res_out = _fill_gaps(residues, gaps)

        record_id = f"{pdb_id}:{chn_id}"
        # ENH - model number in SeqRecord id if multiple models?
        # id = "Chain%s" % str(chain.id)
        # if len(structure) > 1 :
        #     id = ("Model%s|" % str(model.id)) + id

        record = SeqRecord(Seq("".join(res_out)), id=record_id, description=record_id)
        # TODO: Test PDB files with DNA and RNA too:
        record.annotations["molecule_type"] = "protein"

        record.annotations["model"] = model.id
        record.annotations["chain"] = chain.id

        record.annotations["start"] = int(rnumbers[0])
        record.annotations["end"] = int(rnumbers[-1])
        yield record


def _build_seqres_record(chn_id, residues, metadata):
    """Build one SEQRES-derived SeqRecord for a chain (PRIVATE).

    ``metadata`` maps chain id to a list of dicts with the database
    cross-reference information for that chain (from DBREF/DBREF1/DBREF2 PDB
    records, or the equivalent mmCIF records); it may be empty for a chain
    with no such information.
    """
    record = SeqRecord(Seq("".join(residues)))
    record.annotations = {"chain": chn_id}
    # TODO: Test PDB files with DNA and RNA too:
    record.annotations["molecule_type"] = "protein"

    chain_metadata = metadata.get(chn_id)
    if not chain_metadata:
        record.id = chn_id
        return record

    m = chain_metadata[0]
    record.id = record.name = f"{m['pdb_id']}:{chn_id}"
    record.description = f"{m['database']}:{m['db_acc']} {m['db_id_code']}"
    for melem in chain_metadata:
        record.dbxrefs.extend(
            [
                f"{melem['database']}:{melem['db_acc']}",
                f"{melem['database']}:{melem['db_id_code']}",
            ]
        )
    return record


def _build_seqres_records(chains, metadata):
    """Build the sorted-by-chain-id list of SEQRES-derived SeqRecords (PRIVATE)."""
    return [
        _build_seqres_record(chn_id, residues, metadata)
        for chn_id, residues in sorted(chains.items())
    ]


def _parse_seqres_line(line, chains):
    """Extract the chain id and residues from a SEQRES line (PRIVATE)."""
    # NB: We only actually need chain ID and the residues here; commented
    # bits are placeholders from the wwPDB spec.
    # Serial number of the SEQRES record for the current chain.
    # Starts at 1 and increments by one each line.
    # Reset to 1 for each chain.
    # ser_num = int(line[8:10])
    # Chain identifier. This may be any single legal character,
    # including a blank which is used if there is only one
    # chain.
    chn_id = line[11]
    # Number of residues in the chain (repeated on every record)
    # num_res = int(line[13:17])
    residues = [_res2aacode(res) for res in line[19:].split()]
    chains[chn_id].extend(residues)


def _parse_dbref_line(line, metadata):
    """Extract a single-line DBREF record into ``metadata`` (PRIVATE)."""
    #  ID code of this entry (PDB ID)
    pdb_id = line[7:11]
    # Chain identifier.
    chn_id = line[12]
    # Initial sequence number of the PDB sequence segment.
    # seq_begin = int(line[14:18])
    # Initial insertion code of the PDB sequence segment.
    # icode_begin = line[18]
    # Ending sequence number of the PDB sequence segment.
    # seq_end = int(line[20:24])
    # Ending insertion code of the PDB sequence segment.
    # icode_end = line[24]
    # Sequence database name.
    database = line[26:32].strip()
    # Sequence database accession code.
    db_acc = line[33:41].strip()
    # Sequence database identification code.
    db_id_code = line[42:54].strip()
    # Initial sequence number of the database seqment.
    # db_seq_begin = int(line[55:60])
    # Insertion code of initial residue of the segment, if PDB
    # is the reference.
    # db_icode_begin = line[60]
    # Ending sequence number of the database segment.
    # db_seq_end = int(line[62:67])
    # Insertion code of the ending residue of the segment, if
    # PDB is the reference.
    # db_icode_end = line[67]
    metadata[chn_id].append(
        {
            "pdb_id": pdb_id,
            "database": database,
            "db_acc": db_acc,
            "db_id_code": db_id_code,
        }
    )


def _parse_dbref1_line(line):
    """Parse a DBREF1 line, returning its fields for pairing with DBREF2 (PRIVATE)."""
    return {
        # ID code of this entry (PDB ID)
        "pdb_id": line[7:11],
        # Chain identifier.
        "chn_id": line[12],
        # Sequence database name.
        "database": line[26:32].strip(),
        # Sequence database identification code.
        "db_id_code": line[47:67].strip(),
    }


def _parse_dbref2_line(line, dbref, metadata):
    """Parse a DBREF2 line, completing the preceding DBREF1 record (PRIVATE)."""
    pdb_id, chn_id = dbref["pdb_id"], dbref["chn_id"]
    # Ensure ID code and chain are consistent:
    if pdb_id != line[7:11] or chn_id != line[12]:
        raise ValueError("DBREF2 identifiers do not match")
    # Sequence database accession code.
    db_acc = line[18:40].strip()
    metadata[chn_id].append(
        {
            "pdb_id": pdb_id,
            "database": dbref["database"],
            "db_acc": db_acc,
            "db_id_code": dbref["db_id_code"],
        }
    )


def _parse_seqres_stream(stream):
    """Parse SEQRES/DBREF/DBREF1/DBREF2 records from a PDB file stream (PRIVATE).

    Returns ``(chains, metadata)`` where ``chains`` maps chain id to the
    list of one-letter residue codes found in SEQRES lines, and
    ``metadata`` maps chain id to a list of dicts holding the database
    cross-reference information found in DBREF/DBREF1/DBREF2 lines.
    """
    chains = collections.defaultdict(list)
    metadata = collections.defaultdict(list)
    dbref = None  # fields from the most recent DBREF1 line, if any

    rec_name = None
    for line in stream:
        rec_name = line[0:6].strip()
        if rec_name == "SEQRES":
            _parse_seqres_line(line, chains)
        elif rec_name == "DBREF":
            _parse_dbref_line(line, metadata)
        elif rec_name == "DBREF1":
            dbref = _parse_dbref1_line(line)
        elif rec_name == "DBREF2":
            _parse_dbref2_line(line, dbref, metadata)
        # ENH: 'SEQADV' 'MODRES'

    if rec_name is None:
        raise ValueError("Empty file.")

    return chains, metadata


class PdbSeqresIterator(SequenceIterator):
    """Parser for PDB files."""

    modes = "t"

    def __init__(self, source: _TextIOSource) -> None:
        """Iterate over chains in a PDB file as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file

        The sequences are derived from the SEQRES lines in the
        PDB file header, not the atoms of the 3D structure.

        Specifically, these PDB records are handled: DBREF, DBREF1, DBREF2, SEQADV, SEQRES, MODRES

        See: http://www.wwpdb.org/documentation/format23/sect3.html

        This gets called internally via Bio.SeqIO for the SEQRES based interpretation
        of the PDB file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("PDB/1A8O.pdb", "pdb-seqres"):
        ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...     print(record.dbxrefs)
        ...
        Record id 1A8O:A, chain A
        ['UNP:P12497', 'UNP:POL_HV1N5']

        Equivalently,

        >>> with open("PDB/1A8O.pdb") as handle:
        ...     for record in PdbSeqresIterator(handle):
        ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...         print(record.dbxrefs)
        ...
        Record id 1A8O:A, chain A
        ['UNP:P12497', 'UNP:POL_HV1N5']

        Note the chain is recorded in the annotations dictionary, and any PDB DBREF
        lines are recorded in the database cross-references list.
        """
        super().__init__(source, fmt="PDB")
        self.cache = None

    def __next__(self):
        """Iterate over the records in the PDB file."""
        if self.cache is None:
            chains, metadata = _parse_seqres_stream(self.stream)
            self.cache = _build_seqres_records(chains, metadata)
        if not self.cache:
            raise StopIteration
        return self.cache.pop(0)


class PdbAtomIterator(SequenceIterator):
    """Parser for structures in a PDB files."""

    modes = "t"

    def __init__(self, source: _TextIOSource) -> None:
        """Iterate over structures in a PDB file as SeqRecord objects.

        Argument source is a file-like object or a path to a file.

        The sequences are derived from the 3D structure (ATOM records), not the
        SEQRES lines in the PDB file header.

        Unrecognised three letter amino acid codes (e.g. "CSD") from HETATM
        entries are converted to "X" in the sequence.

        In addition to information from the PDB header (which is the same for
        all records), the following chain specific information is placed in the
        annotation:

        record.annotations["residues"] = List of residue ID strings
        record.annotations["chain"] = Chain ID (typically A, B ,...)
        record.annotations["model"] = Model ID (typically zero)

        Where amino acids are missing from the structure, as indicated by
        residue numbering, the sequence is filled in with 'X' characters to
        match the size of the missing region, and  None is included as the
        corresponding entry in the list record.annotations["residues"].

        This function uses the Bio.PDB module to do most of the hard work. The
        annotation information could be improved but this extra parsing should
        be done in parse_pdb_header, not this module.

        This gets called internally via Bio.SeqIO for the atom based
        interpretation of the PDB file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("PDB/1A8O.pdb", "pdb-atom"):
        ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...
        Record id 1A8O:A, chain A

        Equivalently,

        >>> with open("PDB/1A8O.pdb") as handle:
        ...     for record in PdbAtomIterator(handle):
        ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...
        Record id 1A8O:A, chain A

        """
        # NB: SequenceIterator.__init__ is intentionally not called here: this
        # iterator hands ``source`` straight to Bio.PDB, which manages its own
        # file handle, and it eagerly materialises every record below instead
        # of lazily reading from ``self.stream``.
        # pylint: disable=super-init-not-called
        # TODO - Add record.annotations to the doctest, esp the residues (not working?)

        # Only import PDB when needed, to avoid/delay NumPy dependency in SeqIO
        from Bio.PDB.PDBParser import PDBParser

        structure = PDBParser().get_structure(None, source)
        pdb_id = structure.header["idcode"]
        if not pdb_id:
            warnings.warn(
                "'HEADER' line not found; can't determine PDB ID.",
                BiopythonParserWarning,
            )
            pdb_id = "????"

        records = []
        for record in AtomIterator(pdb_id, structure):
            # The PDB header was loaded as a dictionary, so let's reuse it all
            record.annotations.update(structure.header)

            # ENH - add letter annotations -- per-residue info, e.g. numbers

            records.append(record)

        self.records = iter(records)

    def __next__(self):
        """Return the next SeqRecord from the PDB ATOM stream."""
        return next(self.records)


PDBX_POLY_SEQ_SCHEME_FIELDS = (
    "_pdbx_poly_seq_scheme.asym_id",  # Chain ID
    "_pdbx_poly_seq_scheme.mon_id",  # Residue type
)

STRUCT_REF_FIELDS = (
    "_struct_ref.id",  # ID of this reference
    "_struct_ref.db_name",  # Name of the database
    "_struct_ref.db_code",  # Code for this entity
    "_struct_ref.pdbx_db_accession",  # DB accession ID of ref
)

STRUCT_REF_SEQ_FIELDS = (
    "_struct_ref_seq.ref_id",  # Pointer to _struct_ref
    "_struct_ref_seq.pdbx_PDB_id_code",  # PDB ID of this structure
    "_struct_ref_seq.pdbx_strand_id",  # Chain ID of the reference
)


def _normalize_mmcif_fields(records):
    """Ensure every field of interest is present in ``records`` as a list (PRIVATE).

    Explicitly convert records to list (See #1533). If an item is not
    present, use an empty list.
    """
    for field in (
        PDBX_POLY_SEQ_SCHEME_FIELDS + STRUCT_REF_SEQ_FIELDS + STRUCT_REF_FIELDS
    ):
        if field not in records:
            records[field] = []
        elif not isinstance(records[field], list):
            records[field] = [records[field]]


def _mmcif_chains(records):
    """Return chain id -> list of one-letter residue codes (PRIVATE)."""
    chains = collections.defaultdict(list)
    for asym_id, mon_id in zip(
        records["_pdbx_poly_seq_scheme.asym_id"],
        records["_pdbx_poly_seq_scheme.mon_id"],
    ):
        chains[asym_id].append(_res2aacode(mon_id))
    return chains


def _mmcif_struct_refs(records):
    """Return _struct_ref id -> {database, db_id_code, db_acc} (PRIVATE)."""
    struct_refs = {}
    for ref_id, db_name, db_code, db_acc in zip(
        records["_struct_ref.id"],
        records["_struct_ref.db_name"],
        records["_struct_ref.db_code"],
        records["_struct_ref.pdbx_db_accession"],
    ):
        struct_refs[ref_id] = {
            "database": db_name,
            "db_id_code": db_code,
            "db_acc": db_acc,
        }
    return struct_refs


def _mmcif_metadata(records):
    """Return chain id -> list of DB cross-reference dicts (PRIVATE).

    Look through _struct_ref_seq records, look up the corresponding
    _struct_ref record, and add an entry to the metadata list for that
    chain. The dict keys mirror those used elsewhere in this module.
    """
    struct_refs = _mmcif_struct_refs(records)
    metadata = collections.defaultdict(list)
    for ref_id, pdb_id, chain_id in zip(
        records["_struct_ref_seq.ref_id"],
        records["_struct_ref_seq.pdbx_PDB_id_code"],
        records["_struct_ref_seq.pdbx_strand_id"],
    ):
        entry = {"pdb_id": pdb_id}
        entry.update(struct_refs[ref_id])
        metadata[chain_id].append(entry)
    return metadata


class CifSeqresIterator(SequenceIterator):
    """Parser for chains in an mmCIF files."""

    modes = "t"

    def __init__(self, source: _TextIOSource) -> None:
        """Iterate over chains in an mmCIF file as SeqRecord objects.

        Argument source is a file-like object or a path to a file.

        The sequences are derived from the _entity_poly_seq entries in the
        mmCIF file, not the atoms of the 3D structure.

        Specifically, these mmCIF records are handled: _pdbx_poly_seq_scheme
        and _struct_ref_seq. The _pdbx_poly_seq records contain sequence
        information, and the _struct_ref_seq records contain database
        cross-references.

        See:
        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/pdbx_poly_seq_scheme.html
        and
        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/struct_ref_seq.html

        This gets called internally via Bio.SeqIO for the sequence-based
        interpretation of the mmCIF file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("PDB/1A8O.cif", "cif-seqres"):
        ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...     print(record.dbxrefs)
        ...
        Record id 1A8O:A, chain A
        ['UNP:P12497', 'UNP:POL_HV1N5']

        Equivalently,

        >>> with open("PDB/1A8O.cif") as handle:
        ...     for record in CifSeqresIterator(handle):
        ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...         print(record.dbxrefs)
        ...
        Record id 1A8O:A, chain A
        ['UNP:P12497', 'UNP:POL_HV1N5']

        Note the chain is recorded in the annotations dictionary, and any mmCIF
        _struct_ref_seq entries are recorded in the database cross-references
        list.
        """
        # NB: SequenceIterator.__init__ is intentionally not called here: this
        # iterator hands ``source`` straight to Bio.PDB, which manages its own
        # file handle, and it eagerly materialises every record below instead
        # of lazily reading from ``self.stream``.
        # pylint: disable=super-init-not-called

        # Only import PDB when needed, to avoid/delay NumPy dependency in SeqIO
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict

        records = MMCIF2Dict(source)
        _normalize_mmcif_fields(records)

        chains = _mmcif_chains(records)
        metadata = _mmcif_metadata(records)

        self.records = iter(_build_seqres_records(chains, metadata))

    def __next__(self):
        """Return the next SeqRecord from the mmCIF chain stream."""
        return next(self.records)


class CifAtomIterator(SequenceIterator):
    """Parser for structures in an mmCIF files."""

    modes = "t"

    def __init__(self, source: _TextIOSource) -> None:
        """Iterate over structures in an mmCIF file as SeqRecord objects.

        Argument source is a file-like object or a path to a file.

        The sequences are derived from the 3D structure (_atom_site.* fields)
        in the mmCIF file.

        Unrecognised three letter amino acid codes (e.g. "CSD") from HETATM
        entries are converted to "X" in the sequence.

        In addition to information from the PDB header (which is the same for
        all records), the following chain specific information is placed in the
        annotation:

        record.annotations["residues"] = List of residue ID strings
        record.annotations["chain"] = Chain ID (typically A, B ,...)
        record.annotations["model"] = Model ID (typically zero)

        Where amino acids are missing from the structure, as indicated by
        residue numbering, the sequence is filled in with 'X' characters to
        match the size of the missing region, and  None is included as the
        corresponding entry in the list record.annotations["residues"].

        This function uses the Bio.PDB module to do most of the hard work. The
        annotation information could be improved but this extra parsing should
        be done in parse_pdb_header, not this module.

        This gets called internally via Bio.SeqIO for the atom based
        interpretation of the PDB file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("PDB/1A8O.cif", "cif-atom"):
        ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...
        Record id 1A8O:A, chain A

        Equivalently,

        >>> with open("PDB/1A8O.cif") as handle:
        ...     for record in CifAtomIterator(handle):
        ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...
        Record id 1A8O:A, chain A

        """
        # NB: SequenceIterator.__init__ is intentionally not called here: this
        # iterator hands ``source`` straight to Bio.PDB, which manages its own
        # file handle, and it lazily delegates to AtomIterator instead of
        # reading from ``self.stream``.
        # pylint: disable=super-init-not-called
        # TODO - Add record.annotations to the doctest, esp the residues (not working?)

        # Only import parser when needed, to avoid/delay NumPy dependency in SeqIO
        from Bio.PDB.MMCIFParser import MMCIFParser

        structure = MMCIFParser().get_structure(None, source)
        pdb_id = structure.header["idcode"]
        if not pdb_id:
            warnings.warn("Could not determine the PDB ID.", BiopythonParserWarning)
            pdb_id = "????"
        self.records = AtomIterator(pdb_id, structure)

    def __next__(self):
        """Return the next SeqRecord from the mmCIF ATOM stream."""
        return next(self.records)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
