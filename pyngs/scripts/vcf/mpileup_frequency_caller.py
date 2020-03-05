import sys
import gzip
import math
import argparse
from typing import TextIO, Dict, AnyStr
import pyngs.vcf


class CallerOptions:

    def __init__(self):
        self.input_tag = "AD"
        self.minimum_frequency = 0.01
        self.minimum_alternate = 0
        self.output_tag = "XAF"
        self.output_header = """FORMAT=<ID={tag},Number=A,Type=Float,Description="Allele Frequency (AD/total)">"""
        self.digits = 5


def call_variants(opt: CallerOptions, instream: TextIO, outstream: TextIO):
    """Call the variants in the input stream."""
    reader = pyngs.vcf.Reader(instream)

    # get the field parser for the allelic depth tag
    parser = reader.field("FORMAT", opt.input_tag)
    header = reader.header
    header.append(
        pyngs.vcf.Header(opt.output_header.format(tag=opt.output_tag)))
    writer = pyngs.vcf.Writer(outstream, header, reader.samples)

    # foreach variant in the reader
    for variant in reader:
        try:
            formatidx = variant.format.index(opt.input_tag)
        except ValueError:
            # could not find the index of tag
            # so we will just continue
            continue
        
        # prepare the output data
        freqstr = ["."] * len(variant.samples)
        keep = False

        # determine the frequencies for each sample
        for sampleidx, sample in enumerate(variant.samples):
            values = list(parser.interpret(sample[formatidx]))
            totaldepth = sum([v for v in values if v is not None])

            # determine the variant frequencies
            freqs = [0.0] * (len(values) - 1)
            for idx, val in enumerate(values):
                if val is None or idx == 0:
                    continue
                freqs[idx - 1] = val / totaldepth
                if freqs[idx - 1] >= opt.minimum_frequency and val >= opt.minimum_alternate:
                    keep = True

            # print the frequencies to a string
            freqstr[sampleidx] = ",".join([str(round(f, opt.digits)) for f in freqs])
        
        # if we keep the variant, add the alternate 
        # allele frequencies and write the variant 
        if keep:
            variant.add_data(opt.output_tag, freqstr)
            writer.write(variant)


def mpileup_frequency_caller(args):
    """Call variants based on the frequency from the mpileup output."""
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
    opt.minimum_alternate = args.minimum_alternate
    opt.minimum_frequency = args.minimum_frequency
    opt.digits = args.digits

    # call_variants
    call_variants(opt, instream, outstream)

    # close the in and output files
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()



if __name__ == '__main__':
    sparser = argparse.ArgumentParser(
        prog="mpileup_to_vcf",
        description="""Call variants based on the frequency of the 
        alternate allele. 
        
        Typically this script should be used between bcftools 
        mpileup and bcftools call --keep-alts to create a 
        proper VCF file with the lowly present variants included""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, default="stdin",
        help="""The input VCF file.""")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, default="stdout",
        help="The output VCF file.")
    sparser.add_argument(
        "-f", "--minimum-frequency", dest="minimum_frequency",
        type=float, default=0.01,
        help="""The minimum alternate allele frequency.""")
    sparser.add_argument(
        "-r", "--minimum-alternate-reads", dest="minimum_alternate",
        type=int, default=0,
        help="""The minimum number of reads supporting an alternate allele.""")
    sparser.add_argument(
        "-d", "--output-digits", dest="digits",
        type=int, default=5,
        help="""The number of decimals to report in the alternate allele frequency.""")
    sparser.add_argument(
        "-t", "--input-tag", dest="input_tag",
        type=str, default="AD",
        help="""The format tag to parse for the reads per genotype.""")

    sparser.set_defaults(func=mpileup_frequency_caller)
    args = sparser.parse_args()
    args.func(args)
