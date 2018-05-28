#!/bin/env python3

import sys
import argparse
import itertools
from functools import total_ordering
import pyngs.alignment
import pyngs.alignment.consensus


"""
A script to determine the consensus entries from a SAM
file annotated with unique molecular indexes (UMIs).
"""


def group_dna_fragments(fragments, distance=20):
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
        aln.is_supplementary_alignment(),
        aln.is_secondary_alignment(),
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
        data = []
        for aln in self.content:
            part = "{chrom}\t{pos}\t{is_first}\t{sup}\t{sec}\t{name}".format(
                chrom=aln.chromosome,
                pos=aln.position,
                is_first=aln.is_first_in_pair(),
                sup=aln.is_supplementary_alignment(),
                sec=aln.is_secondary_alignment(),
                name=aln.name)
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

    def __call__(self):
        """Create the consensus sequences."""
        # for each alignment
        for aln in self.parser:

            # write the header to the output SAM file
            if self.writer.header is None:
                self.writer.set_header(self.parser.header)

            # get the UMI
            umitag = aln.get_tag(self.tag)
            if self.umi is None:
                self.umi = umitag[2]

            # if the UMI is absent write the entry as is
            if umitag is None or aln.is_unmapped():
                self.writer.write(aln)
                continue

            # process the buffer if the umi switches
            if umitag[2] != self.umi:
                self.paired_consensus()
                self.buffer.clear()
                self.umi = umitag[2]

            # add the alignment for the current UMI
            self.buffer.append(aln)

        # process the last buffer entry
        self.paired_consensus()

    def get_fragments(self):
        """Group alignments in DNA fragment."""
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

    def write_alignment(self, aln):
        """Write a single alignment to the output."""
        self.writer.write(aln)

    def write_alignments(self, alignments):
        """Write alignments to the output."""
        for aln in alignments:
            self.write_alignment(aln)

    def paired_consensus(self):
        """
        Determine the consensus alignments.

        the buffer is first sorted based on the readname to obtain the
        pairs, followed by sorting on the first_in_pair.

        Then, the alignment-pairs are sorted based on chromosome and
        position of both.

        Groups are formed based on the position before the alignments
        are fed into the consensus calling.

        If the alignments in a group are all the same, we can still not
        simply report the alignment, as the base quality scores still need
        to be summed*.

        * The sum of the quality scores will likely exceed the maximum score
        that can be encoded in ASCII based encoding. This needs to be taken
        into account either in a tag or in the score itself.
        """
        def alignment_sorter(aln):
            """Sort alignments to based on their DNA fragment."""
            return (
                aln.name,
                aln.is_last_in_pair(),
                aln.is_secondary_alignment(),
                aln.is_supplementary_alignment())

        # only the forward and reverse read can be handled quickly
        if len(self.buffer) == 2:
            self.write_alignments(self.buffer)
            return

        # sort the buffer first on the name
        self.buffer.sort(key=alignment_sorter)

        # get the alignments for the same DNA fragments (read
        # pairs with supplementary/secondary alignments)
        fragments = self.get_fragments()

        # only a single pair of reads are not
        # sufficient for consensus calling
        if len(fragments) == 1:
            self.write_alignments(fragments[0].content)
            return

        # make groups of DNA fragments which are close enough
        # on the genome for consensus calling.
        fragment_count = 0
        for dna_fragment in group_dna_fragments(fragments):
            fragment_count += 1

            # if only 1 fragment is in a group just write it
            if len(dna_fragment) == 1:
                self.write_alignments(dna_fragment[0].content)
                continue

            # get the consensusses of the first, second, etc
            # alignments from the fragments in order.
            pairs = []
            alignment_sets = zip(*[frg.content for frg in dna_fragment])
            for alignments in alignment_sets:

                # determine the consensus
                consensus = pyngs.alignment.consensus.to_alignment_info(
                    pyngs.alignment.consensus.consensus(alignments))

                # generate a consensus alignment
                consname = "{0}:{1}".format(self.umi, fragment_count)
                tags = [
                    ("um", "Z", self.umi),
                    ("rd", "i", len(list(alignments))),
                    ("qu", "Z", ",".join([repr(q) for q in consensus[5]])),
                    ("rp", "Z", consensus[6])]
                consaln = pyngs.alignment.SAMAlignment(
                    consname, alignments[0].flag,
                    alignments[0].chromosome, consensus[0],
                    alignments[0].mapping_quality, consensus[2],
                    "*", 0, 0,
                    consensus[3], consensus[4],
                    tags)
                pairs.append((consensus, consaln))

            # correct the pair information
            if len(pairs) >= 2:
                if not pairs[0][1].is_unmapped() and not pairs[1][1].is_unmapped():
                    if pairs[0][1].chromosome != pairs[1][1].chromosome:
                        pairs[0][1].mate_chromosome = pairs[1][1].chromosome
                        pairs[0][1].mate_position = pairs[1][1].position

                        pairs[1][1].mate_chromosome = pairs[0][1].chromosome
                        pairs[1][1].mate_position = pairs[0][1].position
                    else:
                        pairs[0][1].mate_chromosome = "="
                        pairs[1][1].mate_chromosome = "="
                        pairs[0][1].mate_position = pairs[1][1].position
                        pairs[1][1].mate_position = pairs[0][1].position
                        tlen = min(pairs[0][0][0], pairs[1][0][0]) - max(pairs[0][0][1], pairs[1][0][1])
                        if pairs[0][1].position < pairs[1][1].position:
                            pairs[0][1].tlen = tlen
                            pairs[1][1].tlen = tlen * -1
                        else:
                            pairs[0][1].tlen = tlen * -1
                            pairs[1][1].tlen = tlen
            # write the paired consensus alignments
            for _, aln in pairs:
                self.write_alignment(aln)
            
            # sys.stderr.write("{0}\n".format(pairs))
        

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
        consobj()

        # close opened files
        if not outstream.closed and outstream != sys.stdout:
            outstream.close()
        if not instream.closed and instream != sys.stdin:
            instream.close()

    # run the main loop
    main()        
