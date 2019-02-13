import sys
import gzip
import itertools

from pyngs import sam
from pyngs.sam import consensus


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
    if retval:
        yield umi, retval


def fragment_sorter(aln):
    """Get the alignment name and attributes for sorting."""
    return (
        aln.name,
        not aln.is_first_in_pair(),
        aln.is_unmapped(),
        aln.is_supplementary_alignment(),
        aln.is_secondary_alignment(),
        aln.chromosome,
        aln.position)


def fragment_positions(fragment):
    """Get the fragmnent position."""
    return (
        [aln.chromosome for aln in fragment],
        [aln.position for aln in fragment])

def same_fragment_position(pos_a, pos_b, allowed_distance):
    """Do 2 fragments have the same position."""
    if len(pos_a[0]) != len(pos_b[0]):
        return False
    if pos_a[0] != pos_b[0]:
        return False
    if pos_a[1] != pos_b:
        distance = 0
        for pos in zip(pos_a[1], pos_b[1]):
            distance += abs(pos[0] - pos[1])
            if distance > allowed_distance:
                return False

    return True


def within_distance(fragments, max_distance=20):
    """Group the alignments by distance."""
    batch = []
    for fragment in fragments:

        positions = fragment_positions(fragment)

        # only check after the first fragment is added
        if batch:
            if same_fragment_position(batch[0][0], positions, allowed_distance=max_distance):
                yield batch
                batch = []

        # always add the current alignment
        batch.append((positions, fragment))

    # yield the last entries
    if batch:
        yield batch


def mate_sorter(aln):
    """Sort the alignments in order of mate."""
    return (
        aln.is_supplementary_alignment(),
        aln.is_secondary_alignment(),
        not aln.is_first_in_pair())


def set_mates(aset):
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


def make_consensus(self, reader, writer, tag="um", max_distance=20):
    """Run the alignments."""
    consensus_factory = consensus.Consensus()
    for umi, alignments in group_per_umi(reader, tag):

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
        fragments.sort(key=fragment_positions)

        # get the fragments that start within `distance` bases
        # of each other in total. All alignments in the fragments
        # should occur on the same chromosomes in the same order.
        # This also implies that each DNA fragment must have the same
        # number of reads (first of pair, last of pair, supplementary
        # and secondary)
        for idx, (_, fragset) in enumerate(within_distance(
                fragments, max_distance)):

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
                for aln in fragset:
                    writer.write(aln)
                continue

            # when more than 1 dna fragments are present merge these
            name = "{umi}:{idx}".format(umi=umi, idx=idx)
            cfragments = []

            # consider the alignments for the dna fragment in order
            for aln in zip(fragset):

                # get the consensus alignment
                caln = consensus_factory.sam_alignment(aln)
                caln.name = name
                caln.flag = aln[0].flag
                mapq = int(sum([a.mapping_quality for a in aln]) / len(aln))
                caln.mapping_quality = mapq

                # add the consensus alignment to the fragment
                cfragments.append(caln)

            # set the mate information in the fragment
            mfragments = set_mates(cfragments)
            for aln in mfragments:
                writer.write(aln)


def run_consensus(args):
    """Run the consensus calling."""
    # open the files
    instream = sys.stdin
    if args.sam != "stdin":
        if args.sam.endswith(".gz"):
            instream = gzip.open(args.sam, "rt")
        else:
            instream = open(args.sam, "rt")

    outstream = sys.stdout
    if args.out != "stdout":
        if args.out.endswith(".gz"):
            outstream = gzip.open(args.out, "wt")
        else:
            outstream = open(args.out, "wt")
    #
    reader = sam.Reader(instream)
    writer = sam.Writer(outstream, reader.header)
    make_consensus(reader, writer, args.tag, args.distance)

    # close opened files
    if not outstream.closed and outstream != sys.stdout:
        outstream.close()
    if not instream.closed and instream != sys.stdin:
        instream.close()
