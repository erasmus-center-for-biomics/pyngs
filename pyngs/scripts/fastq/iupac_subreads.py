
import argparse
import re
from pyngs import open_stream
from pyngs.fastq import fastq
from pyngs.bases import iupac


def iupac_subreads(args):
    """Run the script."""
    regstr = iupac.to_regexp(args.iupac)
    regexp = re.compile(regstr)
    fastqfmt = "@{0}\n{1}\n+\n{2}\n"

    # open the input files
    with open_stream(args.input, "rt") as instream:
        with open_stream(args.output, "wt") as outstream:

            # foreach read
            for name, seq, qual in fastq(instream):

                # match the sequence
                mtch = regexp.search(seq)
                if mtch is None:
                    continue

                # print the read
                start = mtch.span()[0]
                outstream.write(
                    fastqfmt.format(
                        name, seq[start:], qual[start:]))


if __name__ == "__main__":

    sparser = argparse.ArgumentParser(
        prog="iupac_subread",
        description="""Extract the subread after an IUPAC sequence.""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, default="stdin",
        help="""The input FastQ file.""")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, default="stdout",
        help="The output FastQ file.")
    sparser.add_argument(
        "-s", "--iupac", dest="iupac",
        type=str,
        help="""The iupac sequence to search for.""")
    sparser.set_defaults(func=iupac_subreads)
    args = sparser.parse_args()
    args.func(args)

