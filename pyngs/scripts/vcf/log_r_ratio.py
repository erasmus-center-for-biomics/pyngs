import sys
import gzip
import math
import pyngs.vcf


def add_lrr(instream, outstream, tag="AD", expected=None):
    """Adds the Log-R ratio to a VCF file."""
    reader = pyngs.vcf.Reader(instream)
    assert len(reader.samples) == len(expected)

    # get the field parser for the allelic depth tag
    parser = reader.field("FORMAT", tag)
    
    # prepare the header
    header = reader.header
    header.append(pyngs.vcf.Header("""FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R Ratio values">"""))
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
        lrrs = [None] * len(variant.samples)

        # determine the BAF and LRR per sample
        for sidx, sample in enumerate(variant.samples):
            values = list(parser.interpret(sample[fidx]))
            mtot = sum([v for v in values if v is not None])

            # get the LRR value
            lrrs[sidx] = math.log(mtot / expected[sidx]) / math.log(2)

        # add the BAF and LRR data
        variant.add_data("LRR", lrrs)
        writer.write(variant)


def add_log_r_ratio(args):
    """Add the log R ratio to a VCF."""
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
    
    # add the LRR value
    add_lrr(instream, outstream, args.tag, args.expected)