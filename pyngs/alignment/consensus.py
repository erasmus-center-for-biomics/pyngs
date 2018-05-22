"""A module to calculate the consensus from multiple alignments."""
import sys
import operator
import itertools
from .cigar import Cigar


def segment(alignment=None, aid=0):
    """
    Segment an alignment per base.

    Args:
        alignment - an alignment to segments
    Yields:
        a tuple with the
            [0] reference position
            [1] query position
            [2] operator
            [3] sequence
            [4] quality
    """
    cigar = Cigar(alignment.cigar)
    rpos = alignment.position
    qpos = 0
    for nbase, cop in cigar.operations:
        seq = None
        qual = None
        for _ in range(nbase):
            if Cigar.on_query(cop):
                seq = alignment.sequence[qpos]
                qual = alignment.quality[qpos]

            yield(aid, rpos, qpos, cop, seq, qual)

            if Cigar.on_query(cop):
                qpos += 1
            if Cigar.on_reference(cop):
                rpos += 1


class Consensus(object):
    """
    Merges mutltiple alignments to a single alignment.

    In various protocols unique molecular identifiers (UMI) are used to
    identify PCR duplicates from a single DNA fragment. By merging the
    alignments for a single DNA fragment more information on the DNA
    fragment can be obtained than by simply removing alignments with
    duplicate UMIs. In merging the alignments, one must take into
    account base quality scores, CIGAR strings and the position of
    alignments.

    We can visualize the alignments as follows:

    alignment 1 SSSMMMIMMMMMMMDMMMM_______________
    alignment 2 _____MIMMMMMMM_MMMM_______________
    alignment 3 ____________MM_MMMMMMMM___________
    alignment 4 ________________MMMMMMMMMMMSSSS___
    consensus   SSSMMMIMMMMMMM_MMMMMMMMMMMMSSSS___

    To calculate the consensus, the alignments are split per base in the
    CIGAR operations. Each base segment is represented by a tuple with the
    alignment_id, reference position, query position, CIGAR operation, base,
    and quality character. By sorting the segments on reference, position and
    alignment_id, the consensus can be determined in a single subsequent pass.

    Quality values will be determined as follows:

    qual = sum(qual|base_max) - sum(qual|base_other)

    WARNING:
        Due to the method used, (left-most) soft-clipped bases can have a
        lower quality than the actual alignment content.
    """

    def __init__(self, qualityoffset=32, remove_internal_s=True):
        """
        Initialize the MergeAlignment object.

        Args:
            self - reference to self
            qualityoffset - the offset of the quality score
            remove_internal_s - remove the softclipped bases
                                that fall within an alignment.
        """
        self.qualityoffset = qualityoffset
        self.remove_internal_s = remove_internal_s

    @classmethod
    def __merge_insertions__(cls, group=None):
        """
        Merge the insertions prior to further analysis.

        Args:
            self

        Returns:
            reference position [0], query position* [1], operator [2],
            sequence** [3], and quality** [4]

            * First encountered base
            ** Can be multiple

        """
        retval = []
        for _, grp in itertools.groupby(group, key=operator.itemgetter(0)):

            # get the segments per read at this position
            segments = [aln for aln in grp]
            if len(segments) == 1:
                retval.append(segments[0])
                continue

            # merge insertions with the subsequent real character
            temp = [segments[0][0], segments[0][1], segments[0][2], "", "", ""]
            for seg in segments:
                temp[3] += seg[3]
                temp[4] += seg[4]
                temp[5] += seg[5]
            retval.append(tuple(temp))

        # return the cleaned segment
        return retval

    def __score_segments__(self, segments=None, qualdefault=32):
        """
        Score the segments.

        Args:
            self - self
            segments    - a list of segments sorted by cigar and sequence
            qualoffset  - the quality code offset
            qualdefault - the default quality of a cigar without its own
                          quality

        Return:
            a list with scored segment tuples that are formed as follows
            [0] reference coordinate
            [1] cigar
            [2] sequence
            [3] depth
            [4] quality
            [5] ids

        """
        retval = []

        for k, grp in itertools.groupby(segments, key=operator.itemgetter(3, 4)):
            depth = 0
            qual = 0
            ids = []
            base = None
            for seg in grp:
                if base is None:
                    base = [seg[1], seg[3], seg[4]]
                ids.append(seg[0])
                depth += 1
                if seg[5] is None:
                    qual += qualdefault
                else:
                    score = sum([ord(char) - self.qualityoffset for char in seg[5]])
                    qual += operator.truediv(
                        score,
                        len(seg[5]))
            scored_segment = (base[0], base[1], base[2], depth, qual, ids)
            retval.append(scored_segment)

        # return the score
        return retval

    def __rle__(self, noncompressed=""):
        """Run length encode the provided string."""
        parts = []
        pchar = None
        cnt = 0
        for char in noncompressed:
            if pchar is not None and char != pchar:
                parts.append((cnt, pchar))
                cnt = 0
            cnt += 1
            pchar = char
        parts.append((cnt, pchar))

        return "".join(["%d%s" % p for p in parts])

    def __top_score__(self, scored_segments=None):
        """Select the top scoring segments."""
        scored_segments.sort(key=operator.itemgetter(3, 4))

        # select the top scoring hit
        hsegment = list(scored_segments.pop())
        for segm in scored_segments:
            hsegment[4] -= segm[4]
        return hsegment

    def __consensus__(self, segments=None):
        """
        Calculate the consensus from a set of singular scored segments.

        Args:
            self - reference to self
            segments - the top segments

        Returns:
            a continuous consensus sequence in the form
            [0]start
            [1] end
            [2] sequence
            [3] quality
            [4] cigar
        """
        retval = []
        seq, qual, cigar = "", "", ""
        rqual = []
        prev = segments[0][0]
        start = prev
        ascii_threshold = 126 - self.qualityoffset

        # iterate over the segments
        for seg in segments:

            # report non connected part separately
            if prev - seg[0] > 1:
                retval.append([start, prev, seq, qual, self.__rle__(cigar)])
                start = seg[0]
                seq, qual, cigar = "", "", ""

            # add the next base
            cigar += seg[1]
            seq += seg[2] if seg[2] is not None else ""
            qtmp = "~"
            if seg[4] <= ascii_threshold:
                qtmp = chr(int(seg[4]) + self.qualityoffset)
            qual += qtmp * len(seg[2]) if seg[2] is not None else ""

            # add the real qualities
            rqual.append(seg[4])

            # set the previous for the next iteration
            prev = seg[0]

        # add the current entry
        retval.append([start, prev, seq, qual, self.__rle__(cigar), rqual])
        if len(retval) == 1:
            return retval[0]
        else:
            pidx = 0
            base = retval[0]
            for idx in range(1, len(retval)):
                # determine the number of n's to add
                couple = "%dN" % (retval[idx][1] - retval[pidx][1])
                base[1] = retval[idx][1]
                base[2] += retval[idx][2]
                base[3] += couple + retval[idx][3]
                base[4] += retval[idx][4]
                pidx = idx
            return base

    def __call__(self, alignments=None):
        """
        Calculate the consensus alignment.

        Args:
            self - a reference to self
            alignments -
        Returns:
            a consensus sequence represented as a list with
            the following fields:
                [0] reference start
                [1] reference end
                [2] sequence
                [3] quality string
                [4] cigar
                [5] real qualities*

            * the real qualities are numbers and include values
            over the ascii encoding. The maximum quality in Sanger
            ascii encoding is 94 (126 - 32).

        """
        # create segments for each reference/query base in the alignments
        segments = []
        for idx in range(len(alignments)):
            segments.extend([x for x in segment(alignments[idx], aid=idx)])
        segments.sort(key=operator.itemgetter(1, 0, 2))

        # remove H and N CIGAR operations prior to the procedure
        segments = [x for x in segments if x[3] not in ("H", "N")]

        # solve the consensus parts
        parts = []
        for _, group in itertools.groupby(segments, key=operator.itemgetter(1)):
            osegments = [x for x in group]

            # merge the segments for insertions and softclips
            msegments = self.__merge_insertions__(osegments)
            msegments.sort(key=operator.itemgetter(3, 4))

            # score the segments
            ssegments = self.__score_segments__(msegments)
            parts.append(ssegments)

        # Filter S CIGARS that are not at the end
        fltsegment = parts
        if self.remove_internal_s and len(parts) > 2:
            fltsegment = [parts[0]]
            for idx in range(1, len(parts)-1):
                entry = [x for x in parts[idx] if "S" not in x[1]]
                fltsegment.append(entry)
            fltsegment.append(parts[-1])

        # get the top hits
        tophits = [self.__top_score__(x) for x in fltsegment]

        # determine the consensus sequence
        consensus = self.__consensus__(tophits)
        return consensus
