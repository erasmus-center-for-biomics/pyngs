#!/bin/env python3

import sys
import argparse
import itertools
import logging
import pyngs.intersect


def parse_gtf_comment(comment, labels=None):
    """Parse the GTF comment."""
    retval = ["" for _ in labels]
    fields = [f.lstrip().rstrip() for f in comment.split(";")]

    for fld in fields:
        fld = fld.replace('"', "")
        try:
            label, value = fld.split(" ", 1)
            idx = labels.index(label)
            retval[idx] = value
        except ValueError:
            pass
    return retval


def strand_any(a, b):
    """Any direction of the strand."""
    return True


def strand_same(a, b):
    """Entries are on the same strand."""
    return a == b


def strand_opposite(a, b):
    """Entries are on the opposite strands."""
    return a != b


class Writer(object):
    """A class to write the output to."""

    def __init__(self, outstream=sys.stdout):
        """Initialize the object."""
        self.outstream = outstream

    def write_header(self, genes=None, samples=None, scores=None):
        """Write the header."""
        self.outstream.write(
            "{type}\t{gene}\t{sample}\t{scores}\n".format(
                type="Type",
                gene="\t".join(genes),
                sample="\t".join(samples),
                scores="\t".join(scores)))

    def __call__(self, *args, **kwargs):
        """Make self callable."""
        gdata, sdata = args[0]
        gene = "\t".join(gdata[1])
        sample = "\t".join(sdata)
        scores = "\t".join([str(x) for x in args[1]])
        self.outstream.write(
            "{type}\t{gene}\t{sample}\t{scores}\n".format(
                type=gdata[0],
                gene=gene,
                sample=sample,
                scores=scores))


class Quantify(object):
    """An object to quantify genes."""

    def __init__(self, stream=sys.stdin, ncheck=10000, types=None,
                 tocount=None, partition=None, strand_check=strand_any,
                 grace=10000000, resultwriter=None):
        """Initialize a new quantification object."""
        self.ncheck = ncheck
        self.types = types
        self.instream = stream
        self.tocount = tocount
        self.partition = partition
        self.strand_check = strand_check
        self.grace = grace
        self.result = resultwriter
        # internals
        self.buffer = []
        self.records = 0

    def parse_bed_comment(self, comment):
        """Parse the comment in the BED file."""
        tocount = [None for _ in self.tocount]
        topart = [None for _ in self.partition]
        fields = comment.split(";")

        for fld in fields:
            try:
                # try to find the label
                label, value = fld.split("=", 1)
                idx = self.tocount.index(label)
                tocount[idx] = value
            except ValueError:
                pass
            # find the label
            try:
                idx = self.partition.index(label)
                topart[idx] = value
            except ValueError:
                pass
        return topart, tocount

    def __call__(self):
        """Make the object callable."""
        self.run()

    def run(self):
        """Process the input stream."""
        intersect = pyngs.intersect.Parser(
            self.instream,
            func_a=pyngs.intersect.bed5,
            func_b=pyngs.intersect.gtf)
        #
        prevchr = None

        # process all the entries in the intersect
        for bedentry, gtfentry, _ in intersect:

            # skip entries we don't need to process
            if not gtfentry[2] in self.types:
                continue

            # skip entries that do not conform the strand criteria
            if not self.strand_check(bedentry[5], gtfentry[6]):
                continue

            # process the buffer every n entries
            if self.records % self.ncheck == 0:
                curpos = int(bedentry[1])
                data = []
                keep = []
                nitem = 0
                for item in self.buffer:
                    if curpos - item[0] > self.grace:
                        data.append(item)
                        nitem += 1
                    else:
                        keep.append(item)
                logging.info(
                    "Processing %d input records (%d records in buffer; %d total; interval check)",
                    nitem, len(keep), self.records)
                self.__process__(data)
                self.buffer = keep

            # clear the buffer on chromosome changes
            if bedentry[0] != prevchr:
                logging.info(
                    "Processing %d input records (%d total; chromosome change)",
                    len(self.buffer), self.records)
                self.__process__(self.buffer)
                self.buffer.clear()

            # parse the bed comment
            topart, tocount = self.parse_bed_comment(bedentry[3])
            ebed = (topart, tocount)

            # parse the gtf comment
            pos = int(gtfentry[4])
            gtftype = gtfentry[2]
            gtfids = parse_gtf_comment(
                gtfentry[8], labels=["gene_id", "gene_name"])
            egtf = (gtftype, gtfids)

            # add the current entry to the buffer
            self.buffer.append((pos, egtf, ebed))

            # increase the numbers records
            self.records += 1
            prevchr = bedentry[0]


        # process the rest of the buffer
        logging.info("Processing %d input records (%d total; final)", len(self.buffer), self.records)
        self.__process__(self.buffer)

    def __process__(self, data):
        """Process the buffer."""
        def sorter(item):
            """Sort data."""
            return (item[1], item[2][0])

        # sort the data per entry
        nrecords = 0
        data.sort(key=sorter)
        for key, values in itertools.groupby(data, key=sorter):
            tocnt = [b[2][1] for b in values]
            scores = []
            for idx in range(len(self.tocount)):
                val = [l[idx] for l in tocnt]
                scores.append(len(set(val)))

            if self.result is not None:
                self.result(key, scores)
            nrecords += 1
        logging.info("%d Output records sent to output", nrecords)


if __name__ == "__main__":

    def main():
        # setup the logger
        logginglevel = logging.DEBUG
        logging.basicConfig(
            level=logginglevel,
            format="[%(asctime)s %(name)s %(process)d %(levelname)-6s]: %(message)s",
            stream=sys.stderr)

        """Run the main program loop."""
        # prepare the argument parser
        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""
            A script to count the number of overlaps
            generated by an intersect of a BED and a
            GTF file. The intersects sould be generated
            by bedtools with the following command
            bedtools intersect -wao -sorted
            -a [BED5 file] -b [GTF file]""")
        parser.add_argument(
            "-i", "--input", dest="input",
            type=str, default="stdin",
            help="The path of the intersect data.")
        parser.add_argument(
            "-o", "--output", dest="output",
            type=str, default="stdout",
            help="The path of the output data.")
        parser.add_argument(
            "-d", "--direction", dest="direction",
            type=str, default="any",
            choices=["any", "same", "opposite"],
            help="How to treat the strand in the quantification.")
        parser.add_argument(
            "-c", "--count", dest="count",
            type=str, nargs="+", default=["READNAME"],
            help="Readname encoded tags to count.")
        parser.add_argument(
            "-p", "--partition", dest="partition",
            type=str, nargs="+",
            help="The fields by which to partition the counts (eg. sample).")
        parser.add_argument(
            "-t", "--type", dest="type",
            type=str, nargs="+", default="exon",
            help="The type of GTF to quantify.")
        parser.add_argument(
            "--check-interval", dest="ncheck",
            type=int, default=1000000,
            help="Check for transcripts to quantify every n entries.")
        parser.add_argument(
            "--grace", dest="grace",
            type=int, default=10000000,
            help="The number of bases after which entries can be processed.")

        # parse the arguments
        args = parser.parse_args()
        logging.info("Reading from %s", args.input)
        logging.info("Writing to %s", args.output)

        # set the input and output streams
        instream = sys.stdin
        if args.input != "stdin":
            instream = open(args.input, "rt")
        outstream = sys.stdout
        if args.output != "stdout":
            outstream = open(args.output, "wt")

        func = strand_any
        if args.direction == "same":
            func = strand_same
        if args.direction == "opposite":
            func = strand_opposite

        # Initialize the writer and write the header
        writer = Writer(outstream)
        writer.write_header(
            ["gene_id", "gene_name"],
            args.partition,
            args.count)

        # run the quantification
        quant = Quantify(
            stream=instream,
            types=args.type,
            ncheck=args.ncheck,
            tocount=args.count,
            partition=args.partition,
            strand_check=func,
            grace=args.grace,
            resultwriter=writer)
        quant.run()

        logging.info("Finished quantification")
        
        # close the streams
        if instream is not sys.stdin and not instream.closed:
            instream.close()
        if outstream is not sys.stdout and not outstream.closed:
            outstream.close()

    # run the main program
    main()
