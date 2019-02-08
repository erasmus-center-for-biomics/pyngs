#!/bin/env python

import argparse
import gzip
import sys

from pyngs import sam


def to_bed(instream, outstream, operations, tags):
    """Convert the CIGAR strings to a BED file."""
    bedline = "{chromosome}\t{start}\t{end}\t{comment}\t{score}\t{strand}\n"
    reader = sam.Reader(instream)
    readgroups = reader.readgroups()

    for alignment in reader:

        # get the strand
        strand = "-" if alignment.is_reverse() else "+"

        taglst = []
        if tags:
            tagfnd = alignment.get_tags(tags)
            for idx, tagname in enumerate(tags):
                if tagname == "__sample__":
                    rgrp = alignment.get_tag("RG")
                    if rgrp:
                        taglst.append("__sample__={0}".format(readgroups[rgrp[2]]))
                    else:
                        taglst.append("__sample__=None")
                else:
                    tagres = tagfnd[idx] if tagfnd[idx] else "None"
                    taglst.append("{tag}={result}".format(
                        tag=tagname,
                        result=tagres[2]))
        tagstr = ";".join(taglst)

        # iterate over the cigar operations
        for cigar in alignment.cigar_regions():

            # skip cigar operations that we are not interested in
            if cigar[3] not in operations:
                continue
            comment = "READNAME={name};COP={op};{tags}".format(
                name=alignment.name,
                op=cigar[3],
                tags=tagstr)
            outstream.write(
                bedline.format(
                    chromosome=cigar[0],
                    start=cigar[1],
                    end=cigar[2],
                    comment=comment,
                    score=alignment.mapping_quality,
                    strand=strand))


def main():
    """Parse the commandline parameters and run to_bed."""
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description="""
        A script to generate BED entries for bases
        covered by the specified CIGAR operation in a SAM
        file.

        Each BED entry will be annotated with the
        read-name and CIGAR operation type in the COMMENT
        column of the BED file. This script is particularly
        useful for determining the bases covered by ATAC and
        RNA-seq data.

        Currently only SAM format is supported, but BAM files
        can be piped in via samtools view.
        """
    )
    parser.add_argument(
        "-s", "--sam", dest="sam",
        type=str, nargs="?", default="stdin",
        help="""The SAM file of which the CIGAR entries will
                be converted to BED entries.""")
    parser.add_argument(
        "-b", "--bed", dest="bed",
        type=str, nargs="?", default="stdout",
        help="The path to the output BED file.")
    parser.add_argument(
        "-c", "--cigar-operations", dest="operations",
        type=str, nargs="+",
        help="The CIGAR operations to convert to BED entries.")
    parser.add_argument(
        "-t", "--tags", dest="tags",
        type=str, nargs="*",
        help="""The BAM tags to add to the comment column in the
                BED entries. Use the keyword sample to add the
                sample as a comment.""")
    args = parser.parse_args()

    # open the in and output files
    instream = sys.stdin
    if args.sam != "stdin":
        if args.sam.endswith(".gz"):
            instream = gzip.open(args.sam, "rt")
        else:
            instream = open(args.sam, "rt")

    outstream = sys.stdout
    if args.bed != "stdout":
        if args.bed.endswith(".gz"):
            outstream = gzip.open(args.bed, "wt")
        else:
            outstream = open(args.bed, "wt")

    # write the BED entries
    to_bed(instream, outstream, args.operations, args.tags)

    # close the streams
    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":

    # run the main program entry point
    main()
