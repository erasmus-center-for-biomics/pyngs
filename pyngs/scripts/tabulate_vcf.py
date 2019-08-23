import sys
import gzip
import math
import pyngs.vcf


def tabulate(instream, outstream, tags=None):
    """Tabulate the tags in a VCF entry."""
    reader = pyngs.vcf.Reader(instream)
    outstream.write(
        "{chrom}\t{position}\t{reference}\t{alternate}\t{sample}\t{values}".format(
            chrom="chrom",
            position="position",
            reference="reference",
            alternate="alternate",
            sample="samples",
            values="\t".join(tags)))

    # foreach variant in the reader
    for variant in reader:
        for sidx, sample in enumerate(variant.samples):
            values = []
            for tag in tags:
                try:
                    values.append(sample[variant.format.index(tag)])
                except ValueError:
                    values.append(".")

            outstream.write(
                "{chrom}\t{position}\t{reference}\t{alternate}\t{sample}\t{values}".format(
                    chrom=variant.chrom,
                    position=variant.position,
                    reference=variant.reference,
                    alternate=variant.alternate,
                    sample=reader.samples[sidx],
                    values="\t".join(values)))


def tabulate_vcf(args):
    """Convert a VCF file to a tab delimited file."""
    instream = sys.stdin
    if args.vcf != "stdin":
        if args.vcf.endswith(".gz"):
            instream = gzip.open(args.vcf, "rt")
        else:
            instream = open(args.vcf, "rt")

    outstream = sys.stdout
    if args.output != "stdout":
        if args.output.endswith(".gz"):
            outstream = gzip.open(args.output, "wt")
        else:
            outstream = open(args.output, "wt")
    
    # add the BAF value
    tabulate(instream, outstream, args.tags)