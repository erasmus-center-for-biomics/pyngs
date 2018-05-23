#!/bin/env python3

import sys
import argparse
import itertools
import pyngs.alignment


"""
A script to determine the consensus entries from a SAM
file annotated with unique molecular indexes (UMIs).
"""


class CreateConsensus(object):
    """An object to create consensus alignments."""

    def __init__(self, instream=sys.stdin, outstream=sys.stdout, tag="um"):
        """Initialize the consensus caller."""
        self.tag = tag
        self.buffer = []
        self.umi = None
        self.parser = pyngs.alignment.SAMParser(instream)
        self.writer = pyngs.alignment.SAMWriter(outstream)

    def run(self):
        """Create the consensus sequences."""
        # for each alignment
        for aln in self.parser:

            # write the header to the output SAM file
            if self.parser.count == 1:
                self.writer.write_header(self.parser.header)

            # get the UMI
            umitag = aln.get_tag(self.tag)

            # if the UMI is absent write the entry as is
            if umitag is None:
                self.writer.write(aln)
                continue

            # process the buffer if the umi switches
            if umitag[2] != self.umi:
                self.__paired_consensus__()
                self.umi = umitag[2]

            # add the alignment for the current UMI
            self.buffer.append(aln)

        # process the last buffer entry
        self.__paired_consensus__()

    def __paired_consensus__(self):
        """
        Determine the consensus alignments.

        the buffer is first sorted based on the readname to obtain the
        pairs, followed by sorting on the first_in_pair.

        Then, the alignment-pairs are sorted based on chromosome and
        position of both.

        Groups are formed based on the position before the alignments
        are fed into the consensus calling.
        """
        def alignment_sorter(aln):
            return (
                aln.name,
                aln.is_last_in_pair(),
                not aln.is_secondary_alignment(),
                not aln.is_supplementary_alignment())

        def pair_sorter(prs):
            return (
                prs[0].chromosome,
                prs[0].position,
                prs[1].chromosome,
                prs[1].position)

        # only the forward and reverse read can be handled quickly
        if len(self.buffer) == 2:
            for aln in self.buffer:
                aln.tags.append(("nr", "i", 1))
                self.writer.write(aln)
            self.buffer.clear()
            return

        # sort the buffer first on the name
        self.buffer.sort(key=alignment_sorter)
        pairs = []

        # group the read pairs
        for _, grp in itertools.groupby(self.buffer, key=lambda aln: aln.name):
            # get the paired alignments together
            prs = list(grp)

            # write single alignments directly to the output
            if len(prs) < 2:
                for aln in prs:
                    self.writer.write(aln)
                continue
            # add targets to be considered
            pairs.append(prs)

        # only a single pair of reads are not
        # sufficient for consensus calling
        if len(pairs) == 1:
            for aln in pairs[0]:
                aln.tags.append(("nr", "i", 1))
                self.writer.write(aln)
                self.buffer.clear()
                return

        # sort the pairs on position
        pairs.sort(key=pair_sorter)
        
        # determine groups of reads at the same approximate locations

        # printer
        sys.stderr.write("--{umi}--\n".format(umi=self.umi))
        for prs in pairs:
            sys.stderr.write("{n}\t{chroma}\t{posa}\t{chromb}\t{posb}\t{namea}\t{nameb}\n".format(
                n=len(prs),
                chroma=prs[0].chromosome,
                chromb=prs[1].chromosome,
                posa=prs[0].position,
                posb=prs[1].position,
                namea=prs[0].name,
                nameb=prs[1].name
            ))
        # at the end empty the buffer for the next round
        # (NOTE list clear is python 3 specific)
        self.buffer.clear()


if __name__ == "__main__":

    def main():
        """Parse commandline arguments and run the script."""
        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""
                A script to determine the consensus entries from a SAM
                file annotated with unique molecular indexes (UMIs).
                """
        )
        parser.add_argument(
            "-s", "--sam", dest="sam",
            type=str, nargs="?", default="stdin",
            help="The input SAM file.")
        parser.add_argument(
            "-t", "--tag", dest="tag",
            type=str, nargs="?", default="um",
            help="The tag-name for the UMI tag.")
        parser.add_argument(
            "-o", "--output", dest="out",
            type=str, nargs="?", default="stdout",
            help="The output SAM file with the consensus sequences.")
        args = parser.parse_args()

        # open the files
        instream = open(args.sam, "r") if args.sam != "stdin" else sys.stdin
        outstream = open(args.out, "w") if args.out != "stdout" else sys.stdout

        # create the consensus sequences per UMI
        consobj = CreateConsensus(instream, outstream, args.tag)
        consobj.run()

        # close opened files
        if not outstream.closed and outstream != sys.stdout:
            outstream.close()
        if not instream.closed and instream != sys.stdin:
            instream.close()
    main()        
