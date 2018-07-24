#!/bin/env python3 

import sys
import argparse
import gzip
import pyngs.sequence


def determine_N_positions(instream=sys.stdin, outstream=sys.stdout):
    """Determine the N positions in reads."""
    readcounter = 0
    for _, seq, _ in pyngs.sequence.FastQ(instream):
        readcounter += 1
        for pos, base in enumerate(seq.upper()):
            if base != "N":
                continue
            outstream.write(
                "{read}\t{pos}\n".format(read=readcounter, pos=pos))


if __name__ == "__main__":

    def main():
        """Run the main program loop."""
        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""
            A script to filter alignments
            based on tags.""")
        parser.add_argument(
            "-i", "--input", dest="input",
            type=str, nargs="?", default="stdin",
            help="A FastQ file for which to determine N positions.")
        parser.add_argument(
            "-o", "--output", dest="output",
            type=str, nargs="?", default="stdout",
            help="A table with the N positions.")
    
        # parse the command line parameters
        args = parser.parse_args()

        # open the input streams for regular and compressed files
        instream = sys.stdin
        if args.input != "stdin":
            if args.input.endswith(".gz"):
                instream = gzip.open(args.input, "rt")
            else:
                instream = open(args.input, "rt")

        outstream = open(args.output, "wt") if args.output != "stdout" else sys.stdout

        # determine all the N positionsin the reads
        determine_N_positions(instream, outstream)

        # close file handles when we are done
        if not outstream.closed and outstream != sys.stdout:
            outstream.close()
        if not instream.closed and instream != sys.stdin:
            instream.close()
    
    # run the program
    main()
