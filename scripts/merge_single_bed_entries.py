#!/bin/env python

import sys
import argparse


__version__ = "1.0"


def bed_iterator(handle=sys.stdin):
    """Get the next BED entry."""
    header = ["chromosome", "start", "end",
              "comment", "score", "strand"]
    for line in handle:
        line = line.rstrip()
        fields = line.split("\t")

        try:
            fields[1] = int(fields[1])
            fields[2] = int(fields[2])
        except ValueError:
            pass

        # return the zipped data
        yield zip(header, fields)


def process_entries(output=sys.stdout, pid=None, data=None):
    """Process the entries in the data buffer."""
    region = get_regions(data)
    output.write("%(chromosome)s\t%(start)d\t%(end)d\t%(comment)s\t1\t%(strand)s\n" % {
        "chromosome": region[0],
        "start": region[1],
        "end": region[2],
        "comment": pid,
        "strand": region[3]
    })


def get_regions(data=None):
    """Get the region(s) covered by the entries."""
    starts = []
    ends = []
    for elem in data:
        starts.append(elem[1][1])
        ends.append(elem[2][1])
    return data[0][0], min(starts), max(ends), data[0][5]
    


def merge_paired_bed_entries(istream=sys.stdin, ostream=sys.stdout):
    """Merge paired bed entries."""
    iterator = bed_iterator(istream)

    # define a buffer with data
    current = []
    previd = None

    # foreach entry in the bed file
    for entry in iterator:
        sfields = entry[3][1].split(";")
        pid = sfields[0]

        # process the entry
        if previd is not None and pid != previd:
            process_entries(ostream, previd, data=current)
            current = []

        # always add a new entry
        current.append(entry)
        previd = pid

    # process whats left in the buffer
    process_entries(ostream, previd, data=current)


if __name__ == "__main__":

    def main():
        """Parse the commandline options and run the script."""
        finput = sys.stdin
        foutput = sys.stdout

        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""
            A script to merge name-sorted BED entries and
            report the most 5- and 3-prime coordinates, if the
            read names are mapped to the same chromosomes.
            """
        )
        parser.add_argument(
            "-i", "--input", dest="input",
            type=str, nargs="?", default="stdin",
            help="The input BED file for which the entries will be merged.")
        parser.add_argument(
            "-o", "--output", dest="output",
            type=str, nargs="?", default="stdout",
            help="The output BED file with the merged entries.")
        args = parser.parse_args()

        # open files if necessary
        if args.input != "stdin":
            finput = open(args.input, "r")

        if args.output != "stdout":
            foutput = open(args.output, "w")

        # run the analysis
        merge_paired_bed_entries(istream=finput, ostream=foutput)

        # close the in and output streams
        if not finput.closed and finput is not sys.stdin:
            finput.close()
        if not foutput.closed and foutput is not sys.stdout:
            foutput.close()

    # run the main loop
    main()
