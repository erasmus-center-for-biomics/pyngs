import sys
import gzip
import argparse
import logging
from typing import TextIO
import pyngs.vcf as vcf
from pyngs import open_stream


def add_frequency_stream(instream: TextIO, outstream: TextIO, digits=5):
    """Add the frequencies to the variants in the stream."""
    # prepare the reader
    reader = vcf.Reader(instream)
    metas = [m for m in reader.meta]

    # errors
    if "AD" not in reader.format.keys():
        raise ValueError("FORMAT/AD not found in the VCF file")

    if "AF" not in reader.format.keys():
        afstr = """##FORMAT=<ID=AF,Type=Float,Number=R,Description="The allele frequencies of the alternate alleles">"""
        afmeta = vcf.Meta.from_str(afstr)
        afformat = vcf.Format.from_meta(afmeta)
        metas.append(afmeta)
    else:
        logging.warning("FORMAT/AF tag is already present in the header only filling in missing values")
    if "DP" not in reader.format.keys():
        dpstr = """##FORMAT=<ID=DP,Type=Integer,Number=1,Description="The total depth at this position">"""
        dpmeta = vcf.Meta.from_str(dpstr)
        dpformat = vcf.Format.from_meta(dpmeta)
        metas.append(dpmeta)
    else:
        logging.warning("FORMAT/DP tag is already present, using previous DP values")
    #
    writer = vcf.Writer(outstream, metas, reader.samples)
    for variant in reader:
        if not variant.fstore.has_tag("AD"):
            # write the new variant
            logging.warning(
                "Variant %s has no FORMAT/AD tag", variant.to_simple_repr())
            writer.write(variant)
            continue

        # add to the format fields
        format_fields = []
        todo = []
        if not variant.fstore.has_tag("AF"):
            format_fields.append(afformat)
            todo.append("AF")

        if not variant.fstore.has_tag("DP"):
            format_fields.append(dpformat)
            todo.append("DP")

        for fmt in format_fields:
            variant.fstore.add_format(fmt)

        # add the data
        for store in variant.fstore.stores:

            # nothing todo skip
            if not todo:
                continue

            # add the fields that are to be calculated
            adval = store["AD"]
            if "DP" in todo:
                dp = 0
                if adval is not None:
                    dp = sum(adval)
                store["DP"] = [dp]

            if "AF" in todo:
                afval = None
                dpval = store["DP"]
                if dpval is not None and dpval[0] > 0:
                    afval = [round(v / dpval[0], digits) for v in adval]
                store["AF"] = afval

        # write the new variant
        writer.write(variant)


def add_frequency(args):
    """Add allele frequencies to the VCF."""
    with open_stream(args.input, "rt") as instream:
        with open_stream(args.output, "wt") as outstream:
            add_frequency_stream(
                instream,
                outstream,
                digits=args.digits)

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
        "-d", "--digits", dest="digits",
        type=int, default=5,
        help="""The number of digits in the frequency.""")


    sparser.set_defaults(func=add_frequency)
    args = sparser.parse_args()
    args.func(args)
