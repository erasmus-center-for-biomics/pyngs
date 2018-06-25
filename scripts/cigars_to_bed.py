#!/bin/env python

import sys
import logging
import argparse
import pyngs.alignment
import pyngs.interval


def cigar_operations_to_bed(samstream=sys.stdin,
                            outstream=sys.stdout,
                            tags_to_add=None,
                            cigar_operations=None):
    """
    Map cigar operations to BED regions.

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
    samfile = pyngs.alignment.SAMParser(samstream)
    cigarcount = 0
    samples = {}

    # process the alignments
    logging.info("Converting CIGARs to BED entries")
    for aln in samfile:

        # parse the samples from the header
        if samfile.records == 1:
            samples = samfile.sample_map()

        # get the strand
        strand = "-" if aln.is_reverse() else "+"

        # prepare the name from the taglist
        tagslist = ["{0}=".format(l) for l in tags_to_add]
        if len(tags_to_add) > 0:
            for tag in aln.tags:
                if tag[0] == "RG":
                    aln.tags.append(("sample", "Z", samples[tag[2]]))
                try:
                    idx = tags_to_add.index(tag[0])
                    tagslist[idx] = "{label}={val}".format(
                        label=tags_to_add[idx], val=tag[2])
                except ValueError:
                    pass

        # record data from the tags in the BED file name
        name = ';'.join(tagslist)

        # get the cigar intervals and map these to the alignments
        for ival_cigar in aln.cigar_regions():
            if ival_cigar.data in cigar_operations:
                cigarcount += 1
                outstream.write("%s\t%d\t%d\tREADNAME=%s;COP=%s;%s\t%d\t%s\n" % (
                    ival_cigar.chromosome,
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
        samfile.count)


if __name__ == "__main__":

    def main():
        """Parse the commandline options and run the script."""
        # setup the logger
        logginglevel = logging.DEBUG
        logging.basicConfig(
            level=logginglevel,
            format="[%(asctime)s %(name)s %(process)d %(levelname)-6s]: %(message)s",
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
            help="""The SAM file of which the CIGAR entries will
                    be converted to BED entries.""")
        parser.add_argument(
            "-o", "--output", dest="output",
            type=str, nargs="?", default="stdout",
            help="The path to the output BED file.")
        parser.add_argument(
            "-c", "--cigar-operation", dest="operations",
            type=str, nargs="+",
            help="The CIGAR operations to convert to BED entries.")
        parser.add_argument(
            "-t", "--tag", dest="tags",
            type=str, nargs="*",
            help="""The BAM tags to add to the comment column in the
                    BED entries. Use the keyword sample to add the 
                    sample as a comment.""")
        args = parser.parse_args()

        # check that we will report CIGAR operations
        if len(args.operations) == 0:
            logging.fatal("No CIGAR operations to keep were specified.")

        # open files if necessary
        if args.sam != "stdin":
            samin = open(args.sam, "r")

        if args.output != "stdout":
            bedout = open(args.output, "w")

        # no tags added
        tags = []
        if args.tags is None:
            tags = args.tags

        # run the main program loop
        cigar_operations_to_bed(
            samin, bedout,
            tags_to_add=tags,
            cigar_operations=args.operations)

        # close the in and output files
        if not bedout.closed and bedout != sys.stdout:
            bedout.close()
        if not samin.closed and samin != sys.stdin:
            samin.close()

    # run the main loop
    main()
