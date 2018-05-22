#!/bin/env python3

import sys
import argparse
import gzip
import pyngs.sequence


class AddUMI(object):
    """An object to add a UMI sequence."""

    def __init__(self, umistream=None, astream=None, tagname="um"):
        """."""
        self.umis = pyngs.sequence.FastQ(umistream)
        self.reads = pyngs.sequence.FastQ(astream)
        self.tag = tagname

    def __iter__(self):
        """Mark self as an iterator."""
        return self

    def __next__(self):
        """Get the next reads."""
        _, useq, _ = next(self.umis)
        aname, aseq, aqual = next(self.reads)

        # add the UMI to the read name
        aname = aname.split(" ", 1)[0]
        aname += ":{tag}={seq}".format(tag=self.tag, seq=useq)
        return aname, aseq, aqual

    def __call__(self, outstream=None):
        """Convert all the reads and write to the output stream."""
        for name, seq, qual in iter(self):
            pyngs.sequence.write_fastq(outstream, name, seq, qual)


if __name__ == "__main__":

    def main():
        """Parse the commandline options and run the script."""
        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""
            A script to add UMI sequences to the
            readnames in FastQ files.
            """
        )
        parser.add_argument(
            "-f", "--fastq", dest="fastq",
            type=str, nargs="?", default="stdin",
            help="The FastQ file with the reads.")
        parser.add_argument(
            "-u", "--umi", dest="umi",
            type=str,
            help="The FastQ file with the UMI sequences.")
        parser.add_argument(
            "-o", "--output", dest="output",
            type=str, nargs="?", default="stdout",
            help="The output FastQ file.")
        parser.add_argument(
            "-t", "--tag", dest="tag",
            type=str, nargs="?", default="um",
            help="The tag-name to use in the header.")

        # parse the command line parameters
        args = parser.parse_args()

        # open the in and output streams for regular and compressed files
        instream = sys.stdin
        if args.fastq != "stdin":
            if args.fastq.endswith(".gz"):
                instream = gzip.open(args.fastq, "rt")
            else:
                instream = open(args.fastq, "r")

        outstream = sys.stdout
        if args.fastq != "stdout":
            if args.output.endswith(".gz"):
                outstream = gzip.open(args.output, "wt")
            else:
                outstream = open(args.output, "w")

        if args.umi.endswith(".gz"):
            umistream = gzip.open(args.umi, "rt")
        else:
            umistream = open(args.umi, "r")

        # add the umis
        addumi = AddUMI(umistream, instream, args.tag)
        addumi(outstream)

        # close the in and output files
        if not outstream.closed and outstream != sys.stdout:
            outstream.close()
        if not instream.closed and instream != sys.stdin:
            instream.close()
        if not umistream.closed and umistream != sys.stdin:
            umistream.close()
    
    main()
