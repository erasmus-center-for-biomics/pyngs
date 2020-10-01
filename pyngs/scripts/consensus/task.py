from itertools import groupby
from pyngs import sam
from pyngs.sam import consensus



def pair_sorter(aln):
    """Get the alignment name and attributes for sorting."""
    return (
        aln.name,
        not aln.first_in_pair,
        aln.unmapped,
        aln.supplementary_alignment,
        aln.secondary_alignment)


def fragment_positions(fragment):
    """Get the fragmnent position."""
    return (
        len(fragment),
        [aln.chromosome for aln in fragment],
        [aln.position for aln in fragment],
        [aln.flag for aln in fragment])


def group_on_distance(fragments, maxd=20):
    """Group fragments based on distance."""

    def pairwise_compare(frg_a, frg_b, maxd=20):
        """Compare 2 fragments pairwise."""
        if len(frg_a) != len(frg_b):
            return False
        dst = 0
        for idx in range(len(frg_a)):
            if frg_a[idx].chromosome != frg_b[idx].chromosome:
                return False
            if frg_a[idx].flag != frg_b[idx].flag:
                return False
            dst += abs(frg_a[idx].position - frg_b[idx].position)
            if dst > maxd:
                return False
        return True

    batch = []
    for frg in fragments:
        if batch and not pairwise_compare(frg, batch[-1], maxd):
            yield batch
            batch = []
        batch.append(frg)
    if batch:
        yield batch


def set_mates(aset):
    """Set the mates in a set of alignments."""

    def mate_sorter(aln):
        """Sort the alignments in order of mate."""
        return (
            aln.supplementary_alignment,
            aln.secondary_alignment,
            not aln.first_in_pair)

    aset.sort(key=mate_sorter)
    if len(aset) == 1:
        return aset
    # consider paired alignments differently
    if aset[0].paired:
        if not aset[1].last_in_pair:
            raise ValueError("second alignment should be last in pair")
        if aset[1].supplementary_alignment or aset[1].secondary_alignment:
            raise ValueError("second alignment should be primary alignment")

        aset[0].mate = aset[1]
        aset[1].mate = aset[0]

        if len(aset) > 2:
            for idx in range(2, len(aset)):
                if aset[idx].first_in_pair:
                    aset[idx].mate = aset[0]
                else:
                    aset[idx].mate = aset[1]
    # consider single read alignments
    else:
        if len(aset) > 1:
            for idx in range(1, len(aset)):
                aset[idx].mate = aset[0]
    #
    return aset


class PairedConsensus:

    def __init__(self):
        self.max_distance = 20
        self.discard = False
        self.consensus_id = "x"
        self.consensus_factory = consensus.Consensus()

    def __call__(self, umi, alignments):
        # sort the alignments based on the dna fragment
        alignments.sort(key=pair_sorter)

        # group the individual alignments into (DNA) fragments
        # (read-pairs + secondary + supplementary alignemnts)
        fragments = []
        for _, fragment in groupby(alignments, key=lambda a: a.name):
            fragments.append(list(fragment))

        # sort the fragments on position
        fragments.sort(key=fragment_positions)

        # get the fragments that start within `distance` bases
        # of each other in total. All alignments in the fragments
        # should occur on the same chromosomes in the same order.
        # This also implies that each DNA fragment must have the same
        # number of reads (first of pair, last of pair, supplementary
        # and secondary)
        for idx, fragset in enumerate(group_on_distance(
                fragments, self.max_distance)):

            # single dna fragment, so no merging required
            if len(fragset) == 1:
                # write alignment here
                if not self.discard:
                    for alnset in fragset:
                        yield alnset
                continue

            # when more than 1 dna fragments are present merge these
            name = "{umi}:{idx}".format(umi=umi, idx=idx)
            cfragments = []

            # consider the alignments for the dna fragment in order
            for aln in zip(*fragset):

                # skip unmapped fragments
                if aln[0].unmapped:
                    caln = aln[0]
                    caln.name = name
                    caln.tags.append(("rd", "i", len(aln)))
                    cfragments.append(caln)
                    continue

                # get the consensus alignment
                caln = self.consensus_factory.sam_alignment(aln)
                caln.name = name
                caln.flag = aln[0].flag
                mapq = int(sum([a.mapping_quality for a in aln]) / len(aln))
                caln.mapping_quality = mapq
                caln.tags.append(("rd", "i", len(aln)))
                caln.tags.append(("RG", "Z", self.consensus_id))

                # add the consensus alignment to the fragment
                cfragments.append(caln)

            # set the mate information in the fragment
            mfragments = set_mates(cfragments)
            yield mfragments
