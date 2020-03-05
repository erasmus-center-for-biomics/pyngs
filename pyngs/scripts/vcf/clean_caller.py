import sys
import gzip
import math
import argparse
from typing import TextIO, Dict, AnyStr
import pyngs.vcf


class CallerOptions:

    def __init__(self):
        self.input_tag = "XAF"

def clean_variants(opt: CallerOptions, instream: TextIO, outstream: TextIO):
    """Call the variants in the input stream."""
    reader = pyngs.vcf.Reader(instream)
    parser = reader.field("FORMAT", opt.input_tag)
    if parser.number != "A":
        raise ValueError(
            "{tag} value does not have number A".format(opt.input_tag))
    
    # get the field parser for the allelic depth tag
    writer = pyngs.vcf.Writer(outstream, reader.header, reader.samples)

    # foreach variant in the reader
    for variant in reader:
        try:
            formatidx = variant.format.index(opt.input_tag)
        except ValueError:
            # could not find the index of tag
            # so we will just continue
            writer.write(variant)
            continue
        
        # determine the number of alternate alleles
        altalleles = 0
        for alt in pyngs.vcf.quote_tokenizer(variant.alternate, sep=","):
            if alt != ".":
                altalleles += 1
        
        # determine the frequencies for each sample
        for sidx, _ in enumerate(variant.samples):
            frqstr = variant.samples[sidx][formatidx]
            if frqstr == ".":
                continue
            frqparts = [f for f in pyngs.vcf.quote_tokenizer(frqstr, sep=",")]
            if altalleles < len(frqparts):
                frqparts = frqparts[0:altalleles]
                variant.samples[sidx][formatidx] = ",".join(frqparts)
    
        # write the variant with the modified 
        # alternate allele frequencies
        writer.write(variant)


def clean_caller(args):
    """Clean the format fields of the variants after bcftools call."""
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

    # prepare the options
    opt = CallerOptions()
    opt.input_tag = args.input_tag

    # call_variants
    clean_variants(opt, instream, outstream)

    # close the in and output files
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()



if __name__ == '__main__':
    sparser = argparse.ArgumentParser(
        prog="clean_caller",
        description="""Clean the variant after BCFtools call. """)
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, default="stdin",
        help="""The input VCF file.""")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, default="stdout",
        help="The output VCF file.")    
    sparser.add_argument(
        "-t", "--input-tag", dest="input_tag",
        type=str, default="XAF",
        help="""The format tag to parse for the reads per genotype.""")

    sparser.set_defaults(func=clean_caller)
    args = sparser.parse_args()
    args.func(args)
