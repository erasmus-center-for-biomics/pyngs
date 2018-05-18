#!/bin/env python

import sys
import argparse
import logging


__version__ = "1.0"


class TableWriter(object):
    """The TableWriter writes the output table."""

    def __init__(self, handle):
        """Initialize the TableWriter."""
        self.handle = handle
        self.handle.write("#chromosome\tstart\tend\tdepth\n")

    def __call__(self, data):
        """Write the next batch of results."""
        for target, depth in data.items():
            self.handle.write("%s\t%s\t%s\t%d\n" % (
                str(target[0]),
                str(target[1]),
                str(target[2]),
                depth
            ))

    def close(self):
        """Close the output file handle."""
        if not self.handle.closed:
            if self.handle is not sys.stdout or self.handle is not sys.stderr:
                self.handle.close()

    def __del__(self):
        """Close the output file handle on delete."""
        self.close()


def overlap_processor(handle=sys.stdin):
    """Yield each line in the overlap file."""
    while True:
        line = next(handle).rstrip()
        if len(line) == 0:
            continue
        if line[0] == "#":
            continue
        fields = line.split("\t")
        part_a = (fields[0], int(fields[1]), int(fields[2]),
                  fields[3], fields[4], fields[5])
        part_b = (fields[6], int(fields[7]), int(fields[8]))
        overlap = int(fields[-1])
        yield part_a, part_b, overlap


def score_buffer(data):
    """
    Process the reads.

    Arguments:
        alignments - a list with alignments to process

    Returns:
        statistics per read group/wbc

    """
    score = 0

    # per data element
    data.sort()
    previous = None
    for current in data:
        # Starting case
        # when pval is None: pval is not the current value,
        # so we can increase the score
        #
        # Follow up cases
        # don't count double values
        #
        if previous != current:
            score += 1
        previous = current

    # return the stats
    return score


def process_transcripts(t_buffer=None, toprocess=None):
    """Process the transcripts in the buffer."""
    retval = {}
    for name in toprocess:
        retval[name] = score_buffer(t_buffer[name])
        del t_buffer[name]
    return retval


def count_per_target(hin=sys.stdin, hout=sys.stdout,
                     grace_distance=1000000, checkinterval=100000):
    """Count the number of reads per transcript."""
    target_buffer = {}
    lasttarget = {}
    entrycount = 0
    writer = TableWriter(hout)
    unmapped = 0
    targets = 0

    for read, target, overlap in overlap_processor(hin):
        # increase the processed read counter
        entrycount += 1

        # add the next read to the buffer if it overlaps with a target
        if overlap != 0:

            if target not in target_buffer.keys():
                target_buffer[target] = []

            # otherwise count the read only if the strand is the same
            target_buffer[target].append(read[3])
            
            # update the last encountered exon for the gene
            lasttarget[target] = target
        else:
            unmapped += 1
        # process the buffer each 100000 reads or if the chromosome changed
        if entrycount % checkinterval == 0:
            # get the genes to process
            toprocess = []
            for tgt, last in lasttarget.items():
                if read[2] - last[2] > grace_distance or read[0] != last[0]:
                    toprocess.append(tgt)
            logging.debug(
                "%d entries processed and resolving %d targets",
                entrycount, len(toprocess))

            # remove the exons that will be processed
            for tgt in toprocess:
                del lasttarget[tgt]

            # process the transcripts to be finalized
            values = process_transcripts(target_buffer, toprocess)
            writer(values)
            targets += len(toprocess)

    # process the remaining entries in the buffer
    values = process_transcripts(target_buffer, target_buffer.keys())
    writer(values)
    targets += len(toprocess)

    # log the last entries
    logging.info("%d entries processed and resolving %d targets with %d unmapped alignments",
                 entrycount, targets, unmapped)


if __name__ == "__main__":

    def main():
        """Parse the commandline options and run the script."""
        # preset the essential variables
        overlapsin = sys.stdin
        countsout = sys.stdout
        grace = 1000000
        interval = 1000000
        debug = False

        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""
            A script to count the number of overlaps
            generated by bedtools intersect -wao -sorted

            """
        )
        parser.add_argument(
            "-o", "--overlaps", dest="overlaps",
            type=str, nargs="?", default="stdin",
            help="The overlaps file from bedtools.")
        parser.add_argument(
            "-c", "--counts", dest="counts",
            type=str, nargs="?", default="stdout",
            help="The path to the output count file.")
        parser.add_argument(
            "-b", "--grace", dest="grace",
            type=int, nargs="?", default=5000000,
            help="the number of bases after which we can safely process a transcript.")
        parser.add_argument(
            "-i", "--interval", dest="interval",
            type=int, nargs="?", default=1000000,
            help="the number of overlaps which are processed before the transcripts are counted.")
        parser.add_argument(
            "-s", "--stranded", dest="stranded",
            type=bool, nargs="?", default=False,
            help="should we perform a stranded analysis.")
        parser.add_argument(
            "-d", "--dedug", dest="debug",
            type=bool, nargs="?", default=False,
            help="report debug messages.")
        args = parser.parse_args()

        # open files if necessary
        if args.overlaps != "stdin":
            overlapsin = open(args.overlaps, "r")

        if args.counts != "stdout":
            countsout = open(args.counts, "w")

        # setup the logger
        logginglevel = logging.INFO
        if args.debug:
            logginglevel = logging.DEBUG
        logging.basicConfig(
            level=logginglevel,
            format='[%(asctime)s %(name)s %(process)d %(levelname)-6s]: %(message)s',
            stream=sys.stderr)

        # count the alignments per transcript
        count_per_target(
            hin=overlapsin,
            hout=countsout,
            grace_distance=args.grace,
            checkinterval=args.interval)

        # close the in and output streams
        if not overlapsin.closed and overlapsin is not sys.stdin:
            overlapsin.close()
        if not countsout.closed and countsout is not sys.stdout:
            countsout.close()

    # run the main loop
    main()