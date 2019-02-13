"""
Methods to generate consensus alignments from
multiple other alignments.
"""
import collections
from operator import itemgetter

from . import quality_to_score, score_to_quality
from . import encode_rle
from . import Alignment
from .cigar import operations, CigarOperation
from .cigar import CIGAR_OPERATIONS_ON_QUERY
from .cigar import CIGAR_OPERATIONS_ON_REFERENCE


SplitOp = collections.namedtuple(
    "SplitOp",
    ["code", "length", "reference", "sequence", "quality"])


def split_operations(alignment):
    """Split operations from an alignment."""

    # for each operation
    for cigarop in operations(alignment):

        # return non reference operations as is
        if cigarop.code not in CIGAR_OPERATIONS_ON_REFERENCE:
            yield SplitOp(
                code=cigarop.code,
                length=cigarop.length,
                reference=list(cigarop.reference),
                sequence=cigarop.sequence,
                quality=cigarop.quality)
        else:
            # split reference operations per base
            for offset in range(0, cigarop.length):
                seq = ""
                qual = ""
                if cigarop.code in CIGAR_OPERATIONS_ON_QUERY:
                    seq = cigarop.sequence[offset]
                    qual = cigarop.quality[offset]
                yield SplitOp(
                    code=cigarop.code,
                    length=1,
                    reference=[
                        cigarop.reference[0],
                        cigarop.reference[1] + offset,
                        cigarop.reference[1] + offset + 1],
                    sequence=seq,
                    quality=qual)


def splitobs_sort(obj):
    """Sort cigar operations for the consensus creation."""
    return (
        obj.reference[1],
        obj.reference[2],
        obj.code,
        obj.length,
        obj.sequence)


def differs(splitop_a, splitop_b):
    """Check whether 2 consensus alignments differ."""
    if splitop_a.reference != splitop_b.reference:
        return True
    if splitop_a.code != splitop_b.code:
        return True
    if splitop_a.length != splitop_b.length:
        return True
    if splitop_a.sequence != splitop_b.sequence:
        return True
    return False


def same_operation(splitops):
    """Batch the same operations in the consensus alignment."""
    batch = []
    for splitop in splitops:
        if batch:
            if differs(batch[0], splitop):
                yield batch
                batch = []
        batch.append(splitop)
    if batch:
        yield batch


def by_position(batches):
    """Group the consensus alignments by position."""
    sets = []
    for batch in batches:
        if sets:
            if sets[0][0].reference != batch[0].reference:
                yield sets
                sets = []
        sets.append(batch)

    if sets:
        yield sets


class Consensus:
    """A class to generate consensus alignments."""

    def __init__(self, quality_offset: int=32, quality_maxvalue=126):
        """Initialize a new Consensus object."""
        self.quality_offset = quality_offset
        self.quality_maxvalue = quality_maxvalue

    def performance(self, batch):
        """Determine the performance."""
        quals = [0.0] * len(batch[0].quality)

        # if we have a sequence
        if batch[0].code in CIGAR_OPERATIONS_ON_QUERY:

            # for each base in the query sequence
            for idx in range(len(batch[0].quality)):

                # for each split operation
                for sop in batch:
                    quals[idx] += quality_to_score(
                        sop.quality[idx],
                        self.quality_offset)
        # return a tuple with the quality
        # and performance of the batch
        return len(batch), quals

    def operations(self, alignments: list):
        """Generate a new consensus alignment."""
        parts = []
        for alignment in alignments:
            for cigarop in split_operations(alignment):
                parts.append(cigarop)
        parts.sort(key=splitobs_sort)

        # determine the preliminary consensus
        preliminary = []
        for batches in by_position(same_operation(parts)):

            # get the performance of the possible entries per
            # position in the consensus
            to_choose = []
            for batch in batches:
                meas = self.performance(batch)
                to_choose.append((meas, batch[0]))

            # append the top hit to the preliminary consensus
            to_choose.sort(key=itemgetter(0), reverse=True)
            preliminary.append(to_choose[0])

        # remove insertions with fewer than half the reads
        # of the surrounding bases, internal soft-clipped and
        # hard-clipped bases.
        cleaned = []
        for idx, (meas, splitop) in enumerate(preliminary):

            # decide to keep or skip insertions
            if splitop.code in "I":
                if idx > 0 and meas[0] < preliminary[idx-1][0][0] / 2:
                    continue
                if idx < len(preliminary) - 1 and meas[0] < preliminary[idx+1][0][0] / 2:
                    continue

            # remove internal clipped bases
            if splitop.code in "SH":
                if idx > 0 and idx < len(preliminary) - 1:
                    continue
            cleaned.append((meas, splitop))

        # add N CIGAR operations for bases covered by the alignment
        # but absent in the reference.
        consensus = []
        for idx, (meas, consop) in enumerate(cleaned):

            # only check after the first
            if idx > 0:
                prevop = cleaned[idx-1][1]

                # insert an N stretch if the previous operation
                # does not end at the current operations
                if prevop.reference[2] != consop.reference[1]:
                    size = consop.reference[1] - prevop.reference[2]
                    insert = SplitOp(
                        code="N", length=size,
                        reference=[
                            consop.reference[0],
                            prevop.reference[2],
                            consop.reference[1]],
                        sequence="", quality="")
                    consensus.append((0, 0.0), insert)
            consensus.append((meas, consop))

        # return the consensus operations
        return consensus

    def to_sam(self, consensus):
        """Convert a list of consensus operations to a SAM alignment."""
        # get the cigar and sequence as a list
        sequence = []
        quality_scores = []
        cigarops = []
        nalntag = []
        for perf, cons in consensus:
            cigarops.extend(cons.code * cons.length)
            sequence.append(cons.sequence)
            quality_scores.extend(perf[1])
            nalntag.append(perf[0])

        # convert the quality to a string
        quality_strings = []
        quality = []
        for score in quality_scores:
            quality_strings.append("{:.{prec}f}".format(score, prec=0))
            quality.append(
                score_to_quality(
                    score, self.quality_offset, self.quality_maxvalue))

        # get the tags
        tags = [
            ("za", "Z", ",".join([str(v) for v in nalntag])),
            ("zq", "Z", ",".join(quality_strings))]

        # get the sequence
        seq = "".join(sequence)
        qual = "".join(quality)
        cigar = "".join(
            ["{0}{1}".format(op[0], op[1]) for op in encode_rle(cigarops)])

        # return a base alignment
        return Alignment(
            "", 0,
            consensus[0][1].reference[0], consensus[0][1].reference[1],
            0, cigar,
            "*", 0, 0,
            seq, qual, tags)

    def sam_alignment(self, alignments):
        """Get the consensus sequence as a SAM alignment."""
        # get the consensus
        consensus = self.operations(alignments)
        if not consensus:
            return None
        return self.to_sam(consensus)
