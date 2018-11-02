import sys
import argparse
import itertools
import pyngs.alignment
import pyngs.alignment.consensus

"""
A script to determine the consensus entries from a SAM
file annotated with unique molecular indexes (UMIs).
"""


def group_per_umi(samparser, tag="um"):
    """Get alignments per UMI from a UMI sorted samparser."""
    retval = []
    umi = None
    for alignment in samparser:

        # get the UMI for the alignment
        curumi = alignment.get_tag(tag)
        if curumi is None:
            yield "", [alignment]

        # if we don't have an UMI, set it
        if umi is None:
            umi = curumi[2]

        # if the UMI switches more
        if umi != curumi[2]:
            yield umi, retval
            umi = curumi[2]
            retval = []

        # add the alignment to the buffer
        retval.append(alignment)

    # yield the last umi with its alignments
    yield umi, retval


def fragment_sorter(aln: pyngs.alignment.SAMAlignment) -> tuple:
    """Get the alignment name and attributes for sorting."""
    return (
        aln.name,
        not aln.is_first_in_pair(),
        aln.is_unmapped(),
        aln.is_supplementary_alignment(),
        aln.is_secondary_alignment(),
        aln.chromosome,
        aln.position
    )


def position_sorter(fragment: list) -> tuple:
    """Get the fragmnent position."""
    return (
        [aln.chromosome for aln in fragment],
        [aln.position for aln in fragment],
    )


def distance_grouper(fragments: list, distance: int = 20) -> list:
    """Group the alignments by distance."""
    retval = []
    bait = []
    for fragment in fragments:

        # only check after the first fragment is added
        if retval:

            # get the values to compare
            nbait = position_sorter(fragment)
            dstcheck = True

            # check if the chromosomes (and number of entries) are the same
            if bait[0] != nbait[0]:
                yield retval
                dstcheck = False
                # reset the loop
                retval = []
                bait = None

            # do the distance check
            if dstcheck:
                # get the distance between all the dna fragments
                cdst = sum([abs(pos[0] - pos[1]) for pos in zip(bait[1], nbait[1])])
                if cdst > distance:
                    yield retval

                    # reset the loop
                    retval = []
                    bait = None

        # always add the current alignment
        retval.append(fragment)
        if not bait:
            bait = position_sorter(fragment)

    # yield the last buffer entry
    yield retval


def rle(lst: list, sep=":") -> str:
    """Encode the values in a list via run-length encoding."""
    retval = []
    for val in lst:
        if retval and val != retval[0]:
            yield "{0}{1}{2}".format(len(retval), sep, retval[0])
            retval = []
        retval.append(val)
    yield "{0}{1}{2}".format(len(retval), sep, retval[0])


def make_consensus(name: str, aln: list, umi: str) -> pyngs.alignment.SAMAlignment:
    """Make a consensus alignment."""
    cons = pyngs.alignment.consensus.to_alignment_info(
        pyngs.alignment.consensus.consensus(aln))

    # prepare the tags
    qstr = ",".join([str(q) for q in cons[5]])
    tags = [
        ("um", "Z", umi),
        ("rd", "i", len(aln)),
        ("qu", "Z", qstr),
        ("rp", "Z", cons[6])]
    rgtag = aln[0].get_tag("RG")
    if rgtag is not None:
        tags.append(rgtag)

    # get the mapping quality
    mapq = int(sum([a.mapping_quality / len(aln) for a in aln]))

    # create the alignment
    caln = pyngs.alignment.SAMAlignment(
        name, aln[0].flag, aln[0].chromosome, cons[0],
        mapq, cons[2], "*", 0, 0,
        cons[3], cons[4], tags)

    return caln


def mate_sorter(aln: pyngs.alignment.SAMAlignment) -> tuple:
    """Sort the alignments in order of mate."""
    return (
        aln.is_supplementary_alignment(),
        aln.is_secondary_alignment(),
        not aln.is_first_in_pair()
    )


def set_mates(aset: list) -> list:
    """Set the mates in a set of alignments."""
    aset.sort(key=mate_sorter)
    if len(aset) == 1:
        return aset
    # consider paired alignments differently
    if aset[0].is_paired():
        if not aset[1].is_last_in_pair():
            raise ValueError("second alignment should be last in pair")
        if aset[1].is_supplementary_alignment() or aset[1].is_secondary_alignment():
            raise ValueError("second alignment should be primary alignment")

        aset[0].set_mate(aset[1])
        aset[1].set_mate(aset[0])

        if len(aset) > 2:
            for idx in range(2, len(aset)):
                if aset[idx].is_first_in_pair():
                    aset[idx].set_mate(aset[0])
                else:
                    aset[idx].set_mate(aset[1])
    # consider single read alignments
    else:
        if len(aset) > 1:
            for idx in range(1, len(aset)):
                aset[idx].set_mate(aset[0])

    return aset


class ConsensusController(object):
    """An object to create consensus alignments."""

    def __init__(self, instream=sys.stdin, tag="um", distance=20):
        """Initialize the consensus controller."""
        self.tag = tag
        self.distance = distance
        self.buffer = []
        self.umi = None
        self.parser = pyngs.alignment.SAMParser(instream)

    def run(self):
        """Run the alignments."""
        for umi, alignments in group_per_umi(self.parser, self.tag):

            # sort the alignments based on the dna fragment
            alignments.sort(key=fragment_sorter)
            # sys.stdout.write("{umi}\t{naln}\n".format(
            #     umi=umi, naln=len(alignments)))

            # group the individual alignments into (DNA) fragments
            # (read-pairs + secondary + supplementary alignemnts)
            fragments = []
            for _, fragment in itertools.groupby(
                    alignments, key=lambda a: a.name):
                fragments.append(list(fragment))

            # sort the fragments on position
            fragments.sort(key=position_sorter)

            # get the fragments that start within `distance` bases
            # of each other in total. All alignments in the fragments
            # should occur on the same chromosomes in the same order.
            # This also implies that each DNA fragment must have the same
            # number of reads (first of pair, last of pair, supplementary
            # and secondary)
            for idx, fragset in enumerate(distance_grouper(
                    fragments, distance=self.distance)):

                # # just some code to print the data
                # sys.stdout.write("{umi}\t{idx}\t{key}\t{nfrag}\n".format(
                #     umi=umi,
                #     idx=idx,
                #     key=str(position_sorter(fragset[0])),
                #     nfrag=len(fragset)
                # ))

                # single dna fragment, so no merging required
                if len(fragset) == 1:
                    # write alignment here
                    continue

                # when more than 1 dna fragments are present merge these
                name = "{umi}:{idx}".format(umi=umi, idx=idx)
                cfragment = []

                # consider the alignments for the dna fragment in order
                for tidx in range(len(fragset[0])):

                    # get the consensus alignment
                    aln = [frg[tidx] for frg in fragset]
                    caln = make_consensus(name, aln, umi)

                    # add the consensus alignment to the fragment

                    cfragment.append(caln)

                # set the mate information in the fragment
                cfragment = set_mates(cfragment)
                if len(cfragment) > 2:
                    for caln in cfragment:
                        print(repr(caln))


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
            "-d", "--distance", dest="distance",
            type=int, nargs="?", default=20,
            help="The allowed distance between the start postion of alignments with the same UMI.")
        parser.add_argument(
            "-o", "--output", dest="out",
            type=str, nargs="?", default="stdout",
            help="The output SAM file with the consensus sequences.")
        args = parser.parse_args()

        # open the files
        instream = open(args.sam, "r") if args.sam != "stdin" else sys.stdin
        outstream = open(args.out, "w") if args.out != "stdout" else sys.stdout

        #
        cobj = ConsensusController(instream, args.tag, args.distance)
        cobj.run()

        # close opened files
        if not outstream.closed and outstream != sys.stdout:
            outstream.close()
        if not instream.closed and instream != sys.stdin:
            instream.close()


    # run the main program loop
    main()
