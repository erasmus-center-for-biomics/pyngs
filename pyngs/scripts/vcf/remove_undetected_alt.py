import sys
import gzip
import argparse
from typing import TextIO, Dict, AnyStr, Generator, List
import pyngs.vcf as vcf

"""
A script to remove unknown alternate allele from mpileup output, such as the <*> allele from the mpileup VCF output.
"""

class FilterAlt:

    def __init__(self, ploidy: int=2):
        """Initialize the object."""
        self.ploidy = ploidy

    def __call__(self, instream: TextIO, outstream: TextIO) -> None:
        """Filter the other alt alleles."""
        reader = vcf.Reader(instream)
        index = vcf.HeaderIndex(reader.header)
        writer = vcf.Writer(outstream, reader.header, reader.samples)

        # pass through the variants
        for var in reader:
            variant = vcf.Variant.from_row(var, index)
            try:
                altidx = variant.alternate.index("<*>")
            except ValueError:
                # just write variants without a <*> allele
                writer.write(variant)
                continue

            # remove the alternate allele
            variant = variant.filter_alternate(altidx, ploidy=self.ploidy)

            # write the modified variant
            writer.write(variant)


def run(args):
    """Run the script."""
    instream = sys.stdin
    if args.input != "stdin":
        if args.input.endswith(".gz"):
            instream = gzip.open(args.input, "rt")
        else:
            instream = open(args.input, "rt")

    outstream = sys.stdout
    if args.output != "stdout":
        if args.output.endswith(".gz"):
            outstream = gzip.open(args.output, "wt")
        else:
            outstream = open(args.output, "wt")

    # filter the alternate alleles
    fltalt = FilterAlt(ploidy=args.ploidy)
    fltalt(instream, outstream)

    # close the in and output files
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()


if __name__ == "__main__":
    sparser = argparse.ArgumentParser(
        prog="remove_undetected_alt",
        description="""Remove the <*> alternate allele from mpileup input.""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, default="stdin",
        help="""The input VCF file.""")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, default="stdout",
        help="The output VCF file.")
    sparser.add_argument(
        "-p", "--ploidy", dest="ploidy",
        type=int, default=2,
        help="The ploidy of the organism")
    sparser.set_defaults(func=run)
    args = sparser.parse_args()
    args.func(args)
