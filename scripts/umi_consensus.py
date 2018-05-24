#!/bin/env python3

import sys
import argparse
import itertools
from functools import total_ordering
import pyngs.alignment


"""
A script to determine the consensus entries from a SAM
file annotated with unique molecular indexes (UMIs).
"""


def consensus_fragments(fragments, distance=20):
    """Generate groups of fragments for the consensus calling."""
    def in_group(frga, frgb):
        """Determine whether a and b are of the same group."""
        if len(frga.content) != len(frgb.content):
            return False
        for alna, alnb in zip(frga.content, frgb.content):
            if alna.chromosome != alnb.chromosome:
                return False
            if abs(alna.position - alnb.position) > distance:
                return False
        return True

    # sort the pairs on position
    fragments.sort()

    confrg = []
    for frg in fragments:

        # set the previous alignment for the first time
        if not confrg:
            confrg.append(frg)
            continue

        # check whether the alignment is part
        # of a new group
        if not in_group(confrg[0], frg):
            yield confrg
            confrg.clear()

        # always add the current alignment
        confrg.append(frg)

    # yield the last alignments in the group as well
    yield confrg

    # final return to signal the end the generator
    return


def sort_alignments(aln):
    """Sort alignments in a fragment."""
    return (
        not aln.is_supplementary_alignment(),
        not aln.is_secondary_alignment(),
        aln.chromosome,
        aln.position
    )


@total_ordering
class DNAFragment(object):
    """A class to represent DNA fragments."""

    def __init__(self, content=None):
        """Initialize a new DNAfragment object."""
        self.content = sorted(content, key=sort_alignments)

    def __eq__(self, other):
        """2 DNA fragments are the same."""
        if len(self.content) != len(other.content):
            return False
        for aln, bln in zip(self.content, other.content):
            if aln.chromosome != bln.chromosome:
                return False
            if aln.position != bln.position:
                return False
            if aln.flag != bln.flag:
                return False
            if aln.sequence != bln.sequence:
                return False
        # yes
        return True

    def __lt__(self, other):
        """Self is less than the other."""
        if len(self.content) < len(other.content):
            return True
        for aln, bln in zip(self.content, other.content):
            if aln.chromosome < bln.chromosome:
                return True
            if aln.position < bln.position:
                return True
            if aln.flag < bln.flag:
                return True
            if aln.sequence < bln.sequence:
                return True
        return False

    def __repr__(self):
        """Represent self as a string."""
        # len(self.content)
        data = []
        for aln in self.content:
            part = "{chrom}\t{pos}\t{sup}\t{sec}".format(
                chrom=aln.chromosome,
                pos=aln.position,
                sup=aln.is_supplementary_alignment(),
                sec=aln.is_secondary_alignment())
            data.append(part)
        return "{n}\t{parts}".format(
            n=len(self.content),
            parts="\t".join(data))


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

    def __group_reads_per_fragment__(self):
        """Group alignments per DNA fragment."""
        retval = []

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
            retval.append(DNAFragment(prs))
        # return the reads per DNA fragment
        return retval

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
            """Sort alignments to based on their DNA fragment."""
            return (
                aln.name,
                aln.is_last_in_pair(),
                not aln.is_secondary_alignment(),
                not aln.is_supplementary_alignment())

        # only the forward and reverse read can be handled quickly
        if len(self.buffer) == 2:
            for aln in self.buffer:
                aln.tags.append(("nr", "i", 1))
                self.writer.write(aln)
            self.buffer.clear()
            return

        # sort the buffer first on the name
        self.buffer.sort(key=alignment_sorter)

        # get the alignments for the same DNA fragments (read
        # pairs with supplementary/secondary alignments)
        fragments = self.__group_reads_per_fragment__()

        # only a single pair of reads are not
        # sufficient for consensus calling
        if len(fragments) == 1:
            for aln in fragments[0].content:
                aln.tags.append(("nr", "i", 1))
                self.writer.write(aln)
                self.buffer.clear()
                return

        # make groups of DNA fragments which are close enough 
        # on the genome for consensus calling.
        sys.stderr.write("--{umi}--\n".format(umi=self.umi))

        concnt = 0
        for consen_frg in consensus_fragments(fragments):
            concnt += 1

            # if only 1 fragment is in a group just write it
            if len(consen_frg) == 1:
                for aln in consen_frg[0].content:
                    aln.tags.append(("nr", "i", 1))
                    aln.tags.append(("gr", "i", concnt))
                    self.writer.write(aln)
                continue

            # otherwise do consensus calling
            sys.stderr.write("{g} ({n}):\n".format(g=concnt, n=len(consen_frg)))
            for frag in consen_frg:
                sys.stderr.write("\t{frag}\n".format(frag=repr(frag)))

            for aln in zip(*[frg.content for frg in consen_frg]):
                if "S" in aln[0].cigar:
                    consensus = pyngs.alignment.consensus(aln)
            
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
