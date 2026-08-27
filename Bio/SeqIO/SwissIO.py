# Copyright 2006-2013,2020 by Peter Cock.
# Revisions copyright 2008-2009 by Michiel de Hoon.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the "swiss" (aka SwissProt/UniProt) file format.

You are expected to use this module via the Bio.SeqIO functions.
See also the Bio.SwissProt module which offers more than just accessing
the sequences as SeqRecord objects.

See also Bio.SeqIO.UniprotIO.py which supports the "uniprot-xml" format.
"""

from Bio import SeqFeature
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import _TextIOSource
from .Interfaces import SequenceIterator

# Reference keys that map directly onto a SeqFeature.Reference attribute;
# any other key not explicitly ignored below is unexpected and raises.
_REFERENCE_KEY_ATTRS = {"PubMed": "pubmed_id", "MEDLINE": "medline_id"}
_IGNORED_REFERENCE_KEYS = {"DOI", "AGRICOLA"}


def _build_dbxrefs(cross_references):
    """Build the list of 'database:accession' cross-references (PRIVATE)."""
    dbxrefs = []
    for cross_reference in cross_references:
        if len(cross_reference) < 2:
            continue
        database, accession = cross_reference[:2]
        dbxref = f"{database}:{accession}"
        if dbxref not in dbxrefs:
            dbxrefs.append(dbxref)
    return dbxrefs


def _build_reference_feature(reference):
    """Convert a SwissProt reference into a SeqFeature.Reference (PRIVATE)."""
    feature = SeqFeature.Reference()
    feature.comment = " ".join("%s=%s;" % k_v for k_v in reference.comments)
    for key, value in reference.references:
        attr = _REFERENCE_KEY_ATTRS.get(key)
        if attr:
            setattr(feature, attr, value)
        elif key not in _IGNORED_REFERENCE_KEYS:
            raise ValueError(f"Unknown key {key} found in references")
    feature.authors = reference.authors
    feature.title = reference.title
    feature.journal = reference.location
    return feature


def _build_annotations(swiss_record):
    """Build the SeqRecord.annotations dict for a SwissProt record (PRIVATE)."""
    annotations = {
        "molecule_type": "protein",
        "accessions": swiss_record.accessions,
        "organism": swiss_record.organism.rstrip("."),
        "taxonomy": swiss_record.organism_classification,
        "ncbi_taxid": swiss_record.taxonomy_id,
    }
    if swiss_record.protein_existence:
        annotations["protein_existence"] = swiss_record.protein_existence
    if swiss_record.created:
        date, version = swiss_record.created
        annotations["date"] = date
        annotations["sequence_version"] = version
    if swiss_record.sequence_update:
        date, version = swiss_record.sequence_update
        annotations["date_last_sequence_update"] = date
        annotations["sequence_version"] = version
    if swiss_record.annotation_update:
        date, version = swiss_record.annotation_update
        annotations["date_last_annotation_update"] = date
        annotations["entry_version"] = version
    if swiss_record.gene_name:
        annotations["gene_name"] = swiss_record.gene_name
    if swiss_record.host_organism:
        annotations["organism_host"] = swiss_record.host_organism
    if swiss_record.host_taxonomy_id:
        annotations["host_ncbi_taxid"] = swiss_record.host_taxonomy_id
    if swiss_record.comments:
        annotations["comment"] = "\n".join(swiss_record.comments)
    if swiss_record.references:
        annotations["references"] = [
            _build_reference_feature(reference)
            for reference in swiss_record.references
        ]
    if swiss_record.keywords:
        annotations["keywords"] = swiss_record.keywords
    return annotations


class SwissIterator(SequenceIterator):
    """Parser to break up a Swiss-Prot/UniProt file into SeqRecord objects."""

    modes = "t"

    def __init__(self, source: _TextIOSource) -> None:
        """Iterate over a Swiss-Prot file and return SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file

        Every section from the ID line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        This parser is for the flat file "swiss" format as used by:
         - Swiss-Prot aka SwissProt
         - TrEMBL
         - UniProtKB aka UniProt Knowledgebase

        For consistency with BioPerl and EMBOSS we call this the "swiss"
        format. See also the SeqIO support for "uniprot-xml" format.

        Rather than calling it directly, you are expected to use this
        parser via Bio.SeqIO.parse(..., format="swiss") instead.
        """
        super().__init__(source, fmt="SwissProt")

    def __next__(self):
        """Return the next SeqRecord from the Swiss-Prot/UniProt stream."""
        swiss_record = SwissProt._read(self.stream)
        if swiss_record is None:
            raise StopIteration
        # Convert the SwissProt record to a SeqRecord
        record = SeqRecord(
            Seq(swiss_record.sequence),
            id=swiss_record.accessions[0],
            name=swiss_record.entry_name,
            description=swiss_record.description,
            features=swiss_record.features,
        )
        record.dbxrefs = _build_dbxrefs(swiss_record.cross_references)
        record.annotations = _build_annotations(swiss_record)
        return record
