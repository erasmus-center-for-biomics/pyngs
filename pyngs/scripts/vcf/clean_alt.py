import sys
import gzip
import argparse
import pyngs.vcf as vcf
from pyngs import open_stream
from pyngs.vcf.tools import FilterAlt


"""
A script to remove unknown alternate allele from mpileup output, such as the <*> allele from the mpileup VCF output.
"""

def clean_alt(args):
    """Run the script."""
    fltalt = FilterAlt(ploidy=args.ploidy)
    with open_stream(args.input, "rt") as instream:
        with open_stream(args.output, "wt") as outstream:
            reader = vcf.Reader(instream)
            writer = vcf.Writer(outstream, reader.meta, reader.samples)
            for variant in reader:
                nvariant = fltalt(variant, remove=["<*>"])
                writer.write(nvariant)


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
    sparser.set_defaults(func=clean_alt)
    args = sparser.parse_args()
    args.func(args)
