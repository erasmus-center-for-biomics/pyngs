"""
A module to provide a FastQ file reader.
"""

import sys


class FastQ(object):
    """A class to iterate over FastQ file streams."""

    def __init__(self, instream=sys.stdin, clean_readname=False):
        """Initialize the FastQ reader."""
        self.instream = instream
        self.total_reads = 0
        self.clean_readname = clean_readname

    def __iter__(self: object):
        """Mark this class as an iterator."""
        return self

    def readname(self, name: str="") -> (str):
        """Process the readname."""
        retval = name[1:]
        if self.clean_readname:
            # remove everything before the first space
            retval = retval.split(" ")[0]
        return retval

    def __next__(self: object) -> (tuple):
        """Get the next entry."""
        # get the readname
        name = self.readname(next(self.instream).rstrip())

        # get the sequence
        sequence = next(self.instream).rstrip()
        next(self.instream)

        # get the quality string
        quality = next(self.instream).rstrip()

        # increment the total reads
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
