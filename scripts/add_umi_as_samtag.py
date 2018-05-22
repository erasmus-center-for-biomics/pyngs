"""
A script to add umis as a tag in a SAM file
"""

import sys
import argparse
import pyngs.alignment


def add_umi_to_sam(instream=sys.stdin, outstream=sys.stdout, tag="um"):
    """."""
    # prepare the parser and alignment writer
    parser = pyngs.alignment.SAMParser(instream)
    writer = pyngs.alignment.SAMWriter(outstream)

    alncnt = 0
    for aln in parser:
        alncnt += 1
        # write header with the first alignment
        if alncnt == 1:
            writer.write_header(parser.header)

        # get the UMI field from the read name
        fields = aln.name.split(":")
        umifld = [fld for fld in fields if fld.startswith(tag)]

        # if the UMI field is found, remove
        # this from the name and add it as
        # a tag in the alignment
        if umifld:
            umi = umifld.split("=", 1)[1]
            aln.name = ":".join(fields.remove(umifld))
            aln.tags.append((tag, "Z", umi))

        # write the (modified) alignment to
        # the output file
        writer.write(aln)


if __name__ == "__main__":

    def main():
        """Run the main program loop."""
        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""
            A script to add UMI sequence to a SAM file from the read name.
            """
        )
        parser.add_argument(
            "-s", "--sam", dest="sam",
            type=str, nargs="?", default="stdin",
            help="The SAM file with the reads.")
        parser.add_argument(
            "-o", "--output", dest="out",
            type=str, nargs="?", default="stdout",
            help="The output FastQ file.")
        parser.add_argument(
            "-t", "--tag", dest="tag",
            type=str, nargs="?", default="um",
            help="The tag-name to retrieve from the header.")

        # parse the command line parameters
        args = parser.parse_args()

        # open the in and output streams for regular and compressed files
        instream = open(args.sam, "r") if args.sam != "stdin" else sys.stdin
        outstream = open(args.out, "r") if args.out != "stdout" else sys.stdout

        # add the UMI as a SAM tag
        add_umi_to_sam(instream, outstream, args.tag)

        # close file handles when we are done
        if not outstream.closed and outstream != sys.stdout:
            outstream.close()
        if not instream.closed and instream != sys.stdin:
            instream.close()
        
    main()
