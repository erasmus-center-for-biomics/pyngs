"""Classes to represent genomic intervals."""

from functools import total_ordering


@total_ordering
class Interval(object):
    """Represents an interval in genome space."""

    __slots__ = ["chromosome", "start", "end"]

    def __init__(self, chrid, start, end):
        """Create a new interval."""
        self.chromosome = chrid
        self.start = start
        self.end = end

    def overlaps_with(self, oth):
        """Self overlaps with oth."""
        if self.chromosome != oth.chromosome:
            return False
        if self.start > oth.end:
            return False
        if self.end < oth.start:
            return False
        return True

    def contained_in(self, oth):
        """Self is contained in oth."""
        if self.chromosome != oth.chromosome:
            return False
        if self.start >= oth.start and self.end <= oth.end:
            return True
        return False

    def __str__(self):
        return "%s:%d-%d" % (str(self.chromosome), self.start, self.end)

    def __repr__(self):
        return "%s:%d-%d" % (str(self.chromosome), self.start, self.end)

    def __eq__(self, other):
        return (self.chromosome, self.start, self.end) == (other.chromosome, other.start, other.end)

    def __lt__(self, other):
        return (self.chromosome, self.start, self.end) < (other.chromosome, other.start, other.end)


class DataInterval(Interval):
    """A specialization of interval to hold data."""

    def __init__(self, chrid, start, end, data):
        """Create a new DataInterval."""
        super(DataInterval, self).__init__(chrid, start, end)
        self.data = data
    
    def __str__(self):
        return super(DataInterval, self).__str__() + ":" + str(self.data)
