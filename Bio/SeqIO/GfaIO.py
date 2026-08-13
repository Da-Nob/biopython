"""Bio.SeqIO support for the Graphical Fragment Assembly format.

This format is output by many assemblers and includes linkage information for
how the different sequences fit together, however, we just care about the
segment (sequence) information.

Documentation:
- Version 1.x: https://gfa-spec.github.io/GFA-spec/GFA1.html
- Version 2.0: https://gfa-spec.github.io/GFA-spec/GFA2.html
"""

import hashlib
import re
import warnings

from Bio import BiopythonWarning
from Bio.Seq import _UndefinedSequenceData
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


from .Interfaces import _TextIOSource
from .Interfaces import SequenceIterator

# Expected value format for each standard tag type (GFA 1.0 spec), used by
# _check_tag_type to validate a tag's value against its declared type.
_TAG_TYPE_PATTERNS = {
    "A": (r"[!-~]", "printable character"),
    "i": (r"[-+]?[0-9]+", "signed integer"),
    "f": (r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?", "float"),
    "Z": (r"[ !-~]+", "printable string"),
    "J": (r"[ !-~]+", "JSON excluding new-line and tab characters"),
    "H": (r"[0-9A-F]+", "byte array in hex format"),
    "B": (
        r"[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+",
        "array of integers or floats",
    ),
}


def _check_tags(seq, tags):
    """Check a segment line's tags for inconsistencies (PRIVATE)."""
    for tag in tags:
        if tag[:2] == "LN":
            # Sequence length
            if not seq:
                # No sequence data, set the sequence length. There is no
                # public Seq API for this, so a private attribute must be
                # set directly.
                seq._data = _UndefinedSequenceData(  # pylint: disable=protected-access
                    int(tag[5:])
                )
            elif int(tag[5:]) != len(seq):
                warnings.warn(
                    f"Segment line has incorrect length. Expected {tag[5:]} but got {len(seq)}.",
                    BiopythonWarning,
                )
        elif tag[:2] == "SH":
            # SHA256 checksum
            checksum = hashlib.sha256(str(seq).encode()).hexdigest()
            if checksum.upper() != tag[5:]:
                warnings.warn(
                    f"Segment line has incorrect checksum. Expected {tag[5:]} but got {checksum}.",
                    BiopythonWarning,
                )


def _check_tag_type(tag_type, value):
    """Warn if a tag's value does not match its declared type (PRIVATE).

    These RegExs are part of the 1.0 standard.
    """
    pattern = _TAG_TYPE_PATTERNS.get(tag_type)
    if pattern is None:
        warnings.warn(f"Tag has invalid type: {tag_type}", BiopythonWarning)
        return
    regex, description = pattern
    if re.fullmatch(regex, value) is None:
        warnings.warn(
            f"Tag has incorrect type. Expected {description}, got {value}.",
            BiopythonWarning,
        )


def _parse_tag(tag):
    """Split one raw tag into (name, type, value) (PRIVATE).

    Warns if the tag's name looks malformed.
    """
    parts = tag.split(":")
    if len(parts) < 3:
        raise ValueError(f"Segment line has invalid tag: {tag}.")
    name, tag_type = parts[0], parts[1]
    value = ":".join(parts[2:])  # tag value may contain : characters
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9]", name) is None:
        warnings.warn(
            f"Tag has invalid name: {name}. Are they tab delimited?",
            BiopythonWarning,
        )
    return name, tag_type, value


def _tags_to_annotations(tags):
    """Build an annotations dictionary from a list of tags (PRIVATE)."""
    annotations = {}
    for tag in tags:
        name, tag_type, value = _parse_tag(tag)
        annotations[name] = (tag_type, value)
        _check_tag_type(tag_type, value)
    return annotations


def _read_segment_line(stream):
    """Return the (line, fields) of the next segment ("S") line (PRIVATE).

    Blank lines are skipped with a warning. Raises StopIteration once the
    stream is exhausted without finding a segment line.
    """
    for line in stream:
        if line == "\n":
            warnings.warn("GFA data has a blank line.", BiopythonWarning)
            continue
        fields = line.strip("\n").split("\t")
        if fields[0] == "S":
            return line, fields
    raise StopIteration


class Gfa1Iterator(SequenceIterator):
    """Parser for GFA 1.x files.

    Documentation: https://gfa-spec.github.io/GFA-spec/GFA1.html
    """

    modes = "t"

    def __init__(
        self,
        source: _TextIOSource,
    ) -> None:
        """Iterate over a GFA file as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
        """
        super().__init__(source, fmt="GFA 1.0")

    def __next__(self):
        """Return the next SeqRecord from the GFA 1.0 stream."""
        line, fields = _read_segment_line(self.stream)
        if len(fields) < 3:
            raise ValueError(
                f"Segment line must have name and sequence fields: {line}."
            )

        seq = Seq(None, length=0) if fields[2] == "*" else Seq(fields[2])

        tags = fields[3:]
        _check_tags(seq, tags)
        annotations = _tags_to_annotations(tags)

        return SeqRecord(seq, id=fields[1], name=fields[1], annotations=annotations)


class Gfa2Iterator(SequenceIterator):
    """Parser for GFA 2.0 files.

    Documentation for version 2: https://gfa-spec.github.io/GFA-spec/GFA2.html
    """

    modes = "t"

    def __init__(
        self,
        source: _TextIOSource,
    ) -> None:
        """Iterate over a GFA file as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
        """
        super().__init__(source, fmt="GFA 2.0")

    def __next__(self):
        """Return the next SeqRecord from the GFA 2.0 stream."""
        line, fields = _read_segment_line(self.stream)
        if len(fields) < 4:
            raise ValueError(
                f"Segment line must have name, length, and sequence fields: {line}."
            )
        try:
            int(fields[2])
        except ValueError:
            raise ValueError(
                f"Segment line must have an integer length: {line}."
            ) from None

        seq = Seq(None, length=0) if fields[3] == "*" else Seq(fields[3])

        tags = fields[4:]
        _check_tags(seq, tags)
        annotations = _tags_to_annotations(tags)

        return SeqRecord(seq, id=fields[1], name=fields[1], annotations=annotations)
