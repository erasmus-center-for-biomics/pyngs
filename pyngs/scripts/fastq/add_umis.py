import sys
import gzip
import argparse
from typing import Generator
from pyngs.fastq import fastq, encode_in_readname


def add_umis_to_readname(args: argparse.Namespace) -> None:
    """Add UMIs to the read names."""
    instream = sys.stdin
    if args.input != "stdin":
        if args.input.endswith(".gz"):
            instream = gzip.open(args.input, "rt")
        else:
            instream = open(args.input, "rt")

    if args.umis.endswith(".gz"):
        umistream = gzip.open(args.umis, "rt")
    else:
        umistream = open(args.umis, "rt")

    outstream = sys.stdout
    if args.output != "stdout":
        if args.output.endswith(".gz"):
            outstream = gzip.open(args.output, "wt")
        else:
            outstream = open(args.output, "wt")

    # shorten the reads
    umiprovider = fastq(umistream)
    for name, sequence, quality in fastq(instream):
        _, umi, _ = next(umiprovider)
        eumi = "{tag}={umi}".format(args.tag, umi)
        name = encode_in_readname(name, [eumi])       
        outstream.write("@{0}\n{1}\n+\n{2}\n".format(
            name, sequence, quality))
    
    # close the in and output files
    umistream.close()
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()


if __name__ == "__main__":

    sparser = argparse.ArgumentParser(
        prog="extract_subread_from_fastq",
        description="""Extract a subread from a FastQ file.""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, default="stdin",
        help="""The input FastQ file.""")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, default="stdout",
        help="The output FastQ file.")
    sparser.add_argument(
        "-u", "--umis", dest="umis",
        type=str,
        help="""The FastQ with the UMIs.""")
    sparser.add_argument(
        "-t", "--tag", dest="tag",
        type=str, default="umi",
        help="""The tag to use for the UMI sequence""")

    sparser.set_defaults(func=add_umis_to_readname)
    args = sparser.parse_args()
    args.func(args)
