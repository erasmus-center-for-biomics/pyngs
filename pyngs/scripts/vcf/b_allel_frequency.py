import sys
import gzip
import math
import pyngs.vcf


def add_baf(instream, outstream, tag="AD"):
    """Adds the B-allel frequency to the VCF entry."""
    reader = pyngs.vcf.Reader(instream)

    # get the field parser for the allelic depth tag
    parser = reader.field("FORMAT", tag)
    
    # prepare the header
    header = reader.header
    header.append(pyngs.vcf.Header("""FORMAT=<ID=BAF,Number=1,Type=Float,Description="B-Allele Frequency">"""))
    writer = pyngs.vcf.Writer(outstream, header, reader.samples)

    # foreach variant in the reader
    for variant in reader:
        try:
            fidx = variant.format.index(tag)
        except ValueError:
            # could not find the index of tag
            # so just add the data
            writer.write(variant)
            continue
        
        # prepare the data
        bafs = [None] * len(variant.samples)

        # determine the BAF and LRR per sample
        for sidx, sample in enumerate(variant.samples):
            values = list(parser.interpret(sample[fidx]))
            mtot = sum([v for v in values if v is not None])

            # get the BAF value
            bafs[sidx] = values[0] / mtot if mtot != 0 else "."

        # add the BAF and LRR data
        variant.add_data("BAF", bafs)
        writer.write(variant)


def add_b_allel_frequency(args):
    """Add the B-allel frequency and the Log R ratios to a VCF."""
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
    add_baf(instream, outstream, args.tag)