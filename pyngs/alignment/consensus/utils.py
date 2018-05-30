def rle(noncompressed=""):
    """Run-lenght encode a list."""
    retval = []
    pchar = None
    cnt = 0
    for char in noncompressed:
        if pchar is not None and char != pchar:
            retval.append((cnt, pchar))
            cnt = 0
        cnt += 1
        pchar = char
    retval.append((cnt, pchar))
    return retval


def encode_quality(score, max_value=126, offset=33):
    """Encode a quality score."""
    val = int(score) + offset
    if val > max_value:
        return chr(max_value)
    if val < offset:
        return chr(offset)
    return chr(val)


def decode_quality(char, offset=33):
    """Decode a quality character."""
    return ord(char) - offset


def aln_generator(cons):
    """Get the segments in order to the alignment."""
    for segment in cons:
        # if the anchored segments are earlier in the sequence than the
        # segment itself, report these first
        if segment.anchored and segment.anchored[0].qpos < segment.qpos:
            for aseg in segment.anchored:
                yield aseg
            yield segment
        else:
            yield segment
            for aseg in segment.anchored:
                yield aseg
    return


def to_alignment_info(cons, offset=32, maxvalue=126):
    """Convert the consensus to an alignment."""
    cigar = []
    sequence = []
    qualities = []
    reads = []
    startpos = None
    endpos = None
    # for each segment
    for segment in aln_generator(cons):
        # get the starting position of the alignment
        if startpos is None:
            startpos = segment.refpos
        endpos = segment.refpos

        # append the list with the values per segment
        cigar.append(segment.operation)
        if segment.sequence is not None:
            sequence.append(segment.sequence)
        if segment.operation != "D":
            qualities.append(int(segment.quality))
        reads.append(segment.content)

    # prepare string representations
    qual = "".join([encode_quality(s, maxvalue, offset) for s in qualities])
    reads = ";".join(["{0},{1}".format(x, y) for x, y in rle(reads)])
    seq = "".join(sequence)
    cig = "".join(["{0}{1}".format(x, y) for x, y in rle(cigar)])

    # return the alignment data.
    return (startpos, endpos, cig, seq, qual, qualities, reads)
