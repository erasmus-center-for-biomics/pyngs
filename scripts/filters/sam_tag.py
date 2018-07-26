#!/bin/env python3

import os
import os.path
import sys
import argparse
import pyngs.alignment


class TagFilter(object):
    """An object to filter on tag value."""

    def __init__(self, samin=sys.stdin, samout=sys.stdout):
        """Initialize a tag filter."""
        self.samin = samin
        self.samout = samout
        self.tag = None
        self.value = None
        self.compare = None
        self.discard = False

    def __call__(self):
        """Process the SAM file."""
        # prepare the parser and alignment writer
        parser = pyngs.alignment.SAMParser(self.samin)
        writer = pyngs.alignment.SAMWriter(self.samout)

        alncnt = 0
        for aln in parser:
            alncnt += 1
            # write header with the first alignment
            if alncnt == 1:
                writer.set_header(parser.header)

            # get the tag
            tag = aln.get_tag(self.tag)

            # if the tag is absent
            if tag is None:
                if not self.discard:
                    writer.write(aln)
                continue
            
            # otherwise compare the tag to the value
            # and write the alignment if the comparison
            # evaluates to true
            if self.compare(tag, self.value):
                writer.write(aln)
            

    @classmethod
    def equals(cls, tag, value):
        """Is the tag value equal to value."""
        try:
            if tag[1] == ("Z", "A", "H", "B"):
                return tag[2] == value
            elif tag[1] == "i":
                return int(tag[2]) == int(value)
            elif tag[1] == "f":
                return float(tag[2]) == float(value)
        except ValueError:
            return False

    @classmethod
    def greater(cls, tag, value):
        """Is the tag value greater than value."""
        try:
            if tag[1] == ("Z", "A", "H", "B"):
                return tag[2] > value
            elif tag[1] == "i":
                return int(tag[2]) > int(value)
            elif tag[1] == "f":
                return float(tag[2]) > float(value)
        except ValueError:
            return False

    @classmethod
    def less(cls, tag, value):
        """Is the tag value less than value."""
        try:
            if tag[1] == ("Z", "A", "H", "B"):
                return tag[2] > value
            elif tag[1] == "i":
                return int(tag[2]) < int(value)
            elif tag[1] == "f":
                return float(tag[2]) < float(value)
        except ValueError:
            return False

if __name__ == "__main__":

    def main():
        """Run the main program loop."""
        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""
            A script to filter alignments
            based on tags.""")
        parser.add_argument(
            "-s", "--sam", dest="sam",
            type=str, nargs="?", default="stdin",
            help="The SAM file with the reads.")
        parser.add_argument(
            "-o", "--output", dest="out",
            type=str, nargs="?", default="stdout",
            help="The output SAM file.")
        parser.add_argument(
            "-t", "--tag", dest="tag",
            type=str, help="The tag-name to filter on.")
        parser.add_argument(
            "-y", "--type", dest="type",
            choices=["equals", "greater", "less"], default="equals",
            type=str, help="The type of filter to use.")
        parser.add_argument(
            "-v", "--value", dest="value",
            type=str, help="The value to filter on.")
        parser.add_argument(
            "--discard-absent", dest="discard", default=False,
            type=bool,
            help="Discard entries for which the tag could not be found.")

        # parse the command line parameters
        args = parser.parse_args()

        # open the in and output streams for regular and compressed files
        instream = open(args.sam, "r") if args.sam != "stdin" else sys.stdin
        outstream = open(args.out, "w") if args.out != "stdout" else sys.stdout

        # initialize the filter and process the alignments
        filter = TagFilter(instream, outstream)

        # set the filter options
        filter.tag = args.tag
        filter.value = args.value
        if args.discard:
            filter.discard = True

        # set the compare function
        if args.type == "equals":
            filter.compare = filter.equals
        elif args.type == "greater":
            filter.compare = filter.greater
        elif args.type == "less":
            filter.compare = filter.less
        
        # run the filter
        filter()

        # close file handles when we are done
        if not outstream.closed and outstream != sys.stdout:
            outstream.close()
        if not instream.closed and instream != sys.stdin:
            instream.close()

    # run the main program loop
    main()
