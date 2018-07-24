#!/bin/env python3

import sys
import argparse

if __name__ == "__main__":

    def main():
        """Run the main program loop."""
        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""A script to demultiplex PE demux files.""")
        parser.add_argument(
            "-i", "--input", dest="input",
            type=str, nargs="?", default="stdin",
            help="The input table.")
        parser.add_argument(
            "-p", "--prefix", dest="prefix",
            type=str, default="prefix_",
            help="The prefix of the output files.")
        # parse the command line parameters
        args = parser.parse_args()
        prefix = args.prefix
        # open the in and output streams for regular and compressed files
        instream = open(args.input, "r") if args.input != "stdin" else sys.stdin

        handles = {}
        for line in instream:
            line = line.rstrip()
            tokens = line.split("\t") 

            if len(tokens) < 9:
                continue

            out = tokens[6]
            if out not in handles.keys():
                handles[out] = (
                    open(prefix + out + "_R1.fastq", "wt"),
                    open(prefix + out + "_R2.fastq", "wt"))
            data = handles[out]
            data[0].write("@{0}\n{1}\n+\n{2}\n".format(tokens[0], tokens[1], tokens[2]))
            data[1].write("@{0}\n{1}\n+\n{2}\n".format(tokens[3], tokens[4], tokens[5]))
        
        for key, value in handles.items():
            value[0].close()
            value[1].close()

    # run the program
    main()
