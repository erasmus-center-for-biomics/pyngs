"""A module to interpret CIGAR strings and their operations."""
import collections


CIGAR_OPERATIONS = ("M", "I", "D", "N", "S", "H", "P", "=", "X")
CIGAR_OPERATIONS_ON_REFERENCE = ("M", "D", "N", "P", "=", "X")
CIGAR_OPERATIONS_ON_QUERY = ("M", "I", "S", "P", "=", "X")

"""
a named tuple to hold annotated cigar operations with the
code, reference coordinates, query coordinates, sequence, and
quality information
"""
CigarOperation = collections.namedtuple("CigarOperation", [
    "code", "length", "reference", "query",
    "sequence", "quality"])


def cigar_operations(cigar: str):
    """
    Get the next operation in the CIGAR string.
    :param cigar: the CIGAR string
    :yield: the next number of bases
            and its operation
    """
    cstart = 0
    cend = 0
    for idx, op in enumerate(cigar):
        if op in CIGAR_OPERATIONS:
            nop = int(cigar[cstart:(cend+1)])
            cstart = idx + 1
            cend = idx + 1
            yield nop, op
        else:
            cend = idx
    return


def to_locations(cigar: str, offset=0):
    """
    Couple the CIGAR operations to locations on the
    reference and query sequences.

    :param cigar: the CIGAR string
    :param offset: the offset on the reference sequence
    :yield: the next operation and its reference
            and query sequence locations
    """
    rbegin = offset
    rend = offset
    qbegin = 0
    qend = 0
    for nop, op in cigar_operations(cigar):
        if op in CIGAR_OPERATIONS_ON_QUERY:
            qbegin = qend
            qend += nop
        if op in CIGAR_OPERATIONS_ON_REFERENCE:
            rbegin = rend
            rend += nop
        yield op, (rbegin, rend), (qbegin, qend)


def operations(alignment):
    """Partition an alignment per cigar operation."""
    rbegin = alignment.position
    rend = alignment.position
    qbegin = 0
    qend = 0
    sequence = ""
    quality = ""
    for nop, op in cigar_operations(alignment.cigar):

        # if the operation is on the query increase the
        # query begin, query end and the sequence and
        # quality
        if op in CIGAR_OPERATIONS_ON_QUERY:
            qbegin = qend
            qend += nop
            sequence = alignment.sequence[qbegin:qend]
            quality = alignment.quality[qbegin:qend]
        else:
            qbegin = qend
            sequence = ""
            quality = ""

        # if the operation is on the reference sequence, increase
        # the reference coordinates
        if op in CIGAR_OPERATIONS_ON_REFERENCE:
            rbegin = rend
            rend += nop
        else:
            rbegin = rend

        # yield the operation
        yield CigarOperation(
            code=op,
            length=nop,
            reference=(alignment.chromosome, rbegin, rend),
            query=(qbegin, qend),
            sequence=sequence,
            quality=quality)

