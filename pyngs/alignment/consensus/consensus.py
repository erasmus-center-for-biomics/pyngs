"""A module to calculate the consensus alignments from multiple alignments."""
from operator import attrgetter
from itertools import groupby
from .segment import Segment, segmenter


def consensus(alignments, quality_offset=32):
    """Get the consensus from a set of alignments.

    The consensus is determined by first segmenting the
    alignments based on their CIGAR operations. Each individual
    operation becomes a segment. These segments are then first
    solved per reference position. Afterwards, soft-clipped,
    insertions and deletions are solved.
    """
    def grouper(aln):
        """Group segments."""
        anchor = []
        for bln in aln.anchored:
            anchor.append((bln.operation, bln.sequence))
        return (
            aln.refpos,
            aln.operation,
            aln.sequence,
            len(aln.anchored),
            anchor)

    # segment the alignments
    segments = segment_alignments(alignments, quality_offset=quality_offset)

    # aggregate the hits
    aggregate = []
    segments.sort(key=grouper)
    for _, atpos in groupby(segments, key=grouper):
        aggregate.append(summarize_segments(atpos))

    # remove internal S anchors from the aggregate
    for idx in range(len(aggregate)):
        if idx > 0 and idx < len(aggregate) - 1:
            for idb in reversed(range(len(aggregate[idx].anchored))):
                if aggregate[idx].anchored[idb].operation == "S":
                    del aggregate[idx].anchored[idb]

    # get the consensus sequence and operations
    consensus_segments = solve_consensus(aggregate)

    # return the consensus segments
    return consensus_segments


def solve_consensus(segments):
    """Solve the consensus from the aggregated segments."""
    retval = []
    # for each segment per segment
    for _, choices in groupby(segments, key=attrgetter("refpos")):

        # get the segments from which to choose
        choices = list(choices)
        choices.sort(key=attrgetter("quality", "content"), reverse=True)

        # get the winning segment
        hit = choices[0]

        # if there are more segments, decrease
        # the score accordingly.
        if len(choices) > 1:
            for cho in choices[1:]:
                hit.quality -= cho.quality
                # hit.content -= cho.content
        retval.append(hit)

    # return the winning (modified) segments
    return retval


def summarize_segments(segments):
    """Summarize segments with the same sequence and anchors."""
    retval = None
    for segment in segments:
        # first hit
        if retval is None:

            # prepare the return value
            retval = Segment(
                segment.refpos, segment.qpos,
                segment.operation, segment.sequence,
                0.0, 0)
            # copy the anchors over
            for anch in segment.anchored:
                retval.anchored.append(
                    Segment(
                        anch.refpos, anch.qpos,
                        anch.operation, anch.sequence,
                        0.0, 0))

        # increase the scores
        retval.content += 1
        if segment.quality is not None:
            retval.quality += segment.quality

        for idx in range(len(segment.anchored)):
            if segment.anchored[idx].quality is not None:
                retval.anchored[idx].quality += segment.anchored[idx].quality
            retval.anchored[idx].content += 1

    # return segment
    return retval


def segment_alignments(alignments, quality_offset=32):
    """Create segments from a set of alignments."""
    segments = []
    for idx in range(len(alignments)):
        anchorbuffer = []
        current = []
        for segment in segmenter(alignments[idx],
                                 quality_offset=quality_offset,
                                 content=idx):

            # insertions and soft clips are added to the last entry
            if segment.operation in ("S", "I"):
                if current:
                    current[-1].anchored.append(segment)
                else:
                    anchorbuffer.append(segment)

            # "normal" case with (mis)-matches and deletions
            elif segment.operation in ("M", "D"):
                # add the anchor buffer to the next normal entry
                if anchorbuffer:
                    segment.anchored.extend(anchorbuffer)
                    anchorbuffer.clear()
                current.append(segment)

            # set extended CIGAR to M
            elif segment.operation in ("=", "X"):
                if anchorbuffer:
                    segment.anchored.extend(anchorbuffer)
                    anchorbuffer.clear()
                segment.operation = "M"
                current.append(segment)

        # if the alignment contained valid CIGAR
        # operations extend the segments
        if current:
            segments.extend(current)

    # return the segments
    return segments
