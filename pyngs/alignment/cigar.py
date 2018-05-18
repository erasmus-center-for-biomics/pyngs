"""A module to interpret CIGAR strings and their operations."""

CIGAR_OPERATIONS = ("M", "I", "D", "N", "S", "H", "P", "=", "X")
CIGAR_OPERATIONS_ON_REFERENCE = ("M", "D", "N", "P", "=", "X")
CIGAR_OPERATIONS_ON_QUERY = ("M", "I", "S", "P", "=", "X")


def next_operation(cigar):
    """Get the next operation in the CIGAR string."""
    cstart = 0
    cend = 0
    for idx in xrange(len(cigar)):
        if cigar[idx] in CIGAR_OPERATIONS:
            nop = int(cigar[cstart:(cend+1)])
            cstart = idx + 1
            cend = idx + 1
            yield nop, cigar[idx]
        else:
            cend = idx
    raise StopIteration


class Cigar(object):
    """a class to represent a CIGAR."""

    # each CIGAR has an operations slot
    __slots__ = ["operations"]

    @classmethod
    def on_reference(cls, opc):
        """CIGAR operation is an operation on the reference sequence."""
        return opc in CIGAR_OPERATIONS_ON_REFERENCE

    @classmethod
    def on_query(cls, opc):
        """CIGAR operation is an operation on the query sequence."""
        return opc in CIGAR_OPERATIONS_ON_QUERY

    def __init__(self, cigar=None, operations=None):
        """Intialize a new CIGAR string."""
        self.operations = []
        if cigar is not None:
            self.operations = [x for x in next_operation(cigar)]
        if operations is not None:
            self.operations = operations

    def __repr__(self):
        """Return the CIGAR as a string."""
        return "".join([str(o) + str(n) for n, o in self.operations])

