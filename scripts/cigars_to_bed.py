#!/bin/env python

import sys
import os
import os.path
import re
import logging
import argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../"))
import pyngs.alignment as alignment
import pyngs.interval as interval

__version__ = "1.0"


def alignment_batch(samiterator, batchsize=100000):
    """
    Get a batch of reads from a iterator.

    Alignments are read from a alignment..AlignmentFactory
    iterator (samiterator) and returned as a list. The
    number of alignments to read can be controlled with the
    batchsize argument.

    Arguments:
        samiterator - alignment.AlignmentFactory iterator
                      that provides alignments
        batchsize - the number of reads to read on each call

    Returns:
        a list with reads.
    """
    batch = []
    chridx = 0
    while True:
        aln = samiterator.next()
        if len(batch) == batchsize or aln.chridx != chridx:
            yield batch
            batch = []
            chridx = aln.chridx
        batch.append(aln)


def alignments_to_bed(samstream=sys.stdin, outstream=sys.stdout,
                      tags_to_add=None, cigar_operations=None):
    """
    Map the alignments to BED regions.

    Arguments:
        samstream - the input stream of the SAM file.
        outstream - the stream to which to write the BED file to
        tags_to_add - the tags to add to the name
        cigar_operations - the cigar operations to keep
    Returns:
        None?
    """
    # a list with the chromosome Ids
    chromosomes = []

    # prepare the SAM iterator
    logging.info("Preparing SAM parsing")
    factory = alignment.AlignmentFactory(chromosome_list=chromosomes)
    cigarcount = 0
    alignmentcount = 0
    samples = {}

    # process the alignments
    logging.info("Converting CIGARs to BED entries")
    for aln in factory.iterator(samstream):

        # parse the samples from the header
        if alignmentcount == 0:
            matcher = re.compile("^.*ID:([^\t]+).*SM:([^\t]+).*$")
            for line in factory.header:
                if line.startswith("@RG"):
                    regm = matcher.search(line)
                    if regm:
                        samples[regm.group(1)] = regm.group(2)

        # increase the aligment counter
        alignmentcount += 1

        # get the strand
        strand = "-" if aln.is_reverse() else "+"

        # prepare the name from the taglist
        tagslist = ["" for _ in tags_to_add]
        if len(tags_to_add) > 0:
            for tag in aln.tags:
                if tag[0] == "RG":
                    aln.tags.append(("__sample__", ":", samples[tag[2]]))
                try:
                    idx = tags_to_add.index(tag[0])
                    tagslist[idx] = tag[2]
                except ValueError:
                    pass

        # record data from the tags in the BED file name
        name = ';'.join(tagslist)

        # get the cigar intervals and map these to the alignments
        for ival_cigar in aln.cigar_regions():
            if ival_cigar.data in cigar_operations:
                cigarcount += 1
                outstream.write("%s\t%d\t%d\t%s;%s;%s\t%d\t%s\n" % (
                    chromosomes[ival_cigar.chromosome],
                    ival_cigar.start,
                    ival_cigar.end,
                    aln.name,
                    ival_cigar.data,
                    name,
                    aln.mapping_quality,
                    strand))
    logging.info(
        "Finished conversion: wrote %d CIGAR entries from %d alignments",
        cigarcount,
        alignmentcount)


if __name__ == "__main__":

    def main():
        """Parse the commandline options and run the script."""
        # setup the logger
        logginglevel = logging.DEBUG
        logging.basicConfig(
            level=logginglevel,
            format='[%(asctime)s %(name)s %(process)d %(levelname)-6s]: %(message)s',
            stream=sys.stderr)

        # default variables
        samin = sys.stdin
        bedout = sys.stdout

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
            help="The SAM file of which the CIGAR entries will be converted to BED entries.")
        parser.add_argument(
            "-o", "--output", dest="output",
            type=str, nargs="?", default="stdout",
            help="The path to the output BED file.")
        parser.add_argument(
            "-c", "--cigar-operation", dest="operations",
            type=str, nargs="+", default=[],
            help="The CIGAR operations to convert to BED entries.")
        parser.add_argument(
            "-t", "--tag", dest="tags",
            type=str, nargs="*", default=[],
            help="The BAM tags to add to the comment column in the BED entries.")
        args = parser.parse_args()

        # check that we will report CIGAR operations
        if len(args.operations) == 0:
            logging.fatal("No CIGAR operations to keep were specified.")

        # open files if necessary
        if args.sam != "stdin":
            samin = open(args.sam, "r")

        if args.output != "stdout":
            bedout = open(args.output, "w")

        # run the main program loop
        alignments_to_bed(
            samin, bedout,
            tags_to_add=list(args.tags),
            cigar_operations=list(args.operations))

        # close the in and output files
        if not bedout.closed and bedout != sys.stdout:
            bedout.close()
        if not samin.closed and samin != sys.stdin:
            samin.close()

    # run the main loop
    main()
