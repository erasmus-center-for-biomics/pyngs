import sys
import gzip
from pyngs.fastq import fastq


def trim_reads(provider, length=50):
    """Shorten the reads from the provider."""
    for name, sequence, quality in provider:
        if len(sequence) > length:
            sequence = sequence[0:length]
            quality = quality[0:length]
        yield name, sequence, quality

def trim_fastq(args):
    """Shorten the reads in a fastq file to length."""
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
    
    # parse the length
    length = args.length
    if length <= 0:
        raise ValueError("The read length cannot be smaller than 1")

    # shorten the reads
    provider = fastq(instream)
    for name, sequence, quality in trim_reads(provider, length):
        outstream.write("@{0}\n{1}\n+\n{2}\n".format(
            name, sequence, quality))
    
    # close the in and output files
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()