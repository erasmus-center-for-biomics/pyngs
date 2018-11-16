"""
A module to provide a FastQ file reader.
"""

import sys


class FastQ(object):
    """A class to iterate over FastQ file streams."""

    def __init__(self, instream=sys.stdin):
        """Initialize the FastQ reader."""
        self.instream = instream
        self.total_reads = 0

    def __iter__(self: object):
        """Mark this class as an iterator."""
        return self

    def __next__(self: object) -> (tuple):
        """Get the next entry."""
        name = next(self.instream).rstrip()[1:]
        sequence = next(self.instream).rstrip()
        next(self.instream)
        quality = next(self.instream).rstrip()
        self.total_reads += 1
        return name, sequence, quality

    def next(self: object) -> (tuple):
        """Make this class compatible with python2."""
        return self.__next__()

    def nreads(self: object) -> (int):
        """Get the number of reads read."""
        return self.total_reads


def write_fastq(outstream=sys.stdout, name="", seq="", qual=""):
    """Write a FastQ entry to a stream."""
    outstream.write("@{name}\n{seq}\n+\n{qual}\n".format(
        name=name,
        seq=seq,
        qual=qual
    ))
