"""A parser for bedtools intersect."""

import sys


def bed5(fields):
    """Parse the BED5 portion of the fields."""
    return fields[:5], fields[5:]


def bed3(fields):
    """Parse the BED3 portion of the fields."""
    return fields[:2], fields[2:]


def gtf(fields):
    """Parse the BED3 portion of the fields."""
    return fields[:8], fields[8:]


class Parser(object):
    """Parser parses intersects."""

    def __init__(self, stream=sys.stdin, sep="\t", func_a=None, func_b=None):
        """Initialize a new intersect parser."""
        self.parser_a = func_a
        self.parser_b = func_b
        self.stream = stream
        self.sep = sep
        if not callable(self.parser_a):
            raise ValueError("Column parsing functions should be callable")
        if not callable(self.parser_b):
            raise ValueError("Column parsing functions should be callable")

    def __iter__(self):
        """Mark self as an iterator."""
        return self

    def __next__(self):
        """Get the next entry."""
        for line in self.stream:
            fields = line.rstrip().split(self.sep)
            resa, fields = self.parser_a(fields)
            resb, fields = self.parser_b(fields)
            return resa, resb, fields

        # raise a stop iteration to signify the end of the file
        raise StopIteration

    def next(self):
        """Get the next entry."""
        return self.__next__()
