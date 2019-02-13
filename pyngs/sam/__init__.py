"""Read, write and otherwise process alignments in SAM format."""
import typing
import re
from .cigar import CIGAR_OPERATIONS, \
    CIGAR_OPERATIONS_ON_REFERENCE, \
    cigar_operations


# SAM flag bits
IS_PAIRED = 0x0001
PROPER_PAIR = 0x0002
QUERY_UNMAPPED = 0x0004
MATE_UNMAPPED = 0x0008
IS_REVERSE = 0x0010
MATE_IS_REVERSE = 0x0020
FIRST_IN_PAIR = 0x0040
LAST_IN_PAIR = 0x0080
NOT_PRIMARY = 0x0100
QC_FAILED = 0x0200
DUPLICATE = 0x0400
SUPPLEMENTARY = 0x0800


class Alignment:
    """A class representing a SAM alignment."""

    __slots__ = [
        "name", "flag", "chromosome", "position",
        "mapping_quality", "cigar", "mate_chromosome",
        "mate_position", "tlen", "sequence", "quality",
        "tags"]

    def __init__(self, *args):
        """Create a new alignment."""
        self.name = args[0]
        self.flag = args[1]
        self.chromosome = args[2]
        self.position = args[3]
        self.mapping_quality = args[4]
        self.cigar = args[5]
        self.mate_chromosome = args[6]
        self.mate_position = args[7]
        self.tlen = args[8]
        self.sequence = args[9]
        self.quality = args[10]
        self.tags = args[11]

    def is_set(self, mask):
        """
        Check whether the bits in the mask are set in the flag.

        :param mask: the bits to check
        :return: True if the bits are set, False if not
        """
        return self.flag | mask == self.flag

    def set_bits(self, mask):
        """
        Set the bit in the mask.

        :param mask: the bits to set
        """
        self.flag |= mask

    def clear_bits(self, mask):
        """
        Clear the bit in the mask.

        :param mask: the bits to clear
        """
        self.flag &= ~mask

    @property
    def paired(self):
        return self.is_set(IS_PAIRED)

    @paired.setter
    def paired(self, value: bool):
        if value:
            self.set_bits(IS_PAIRED)
        else:
            self.clear_bits(IS_PAIRED)

    @property
    def proper_pair(self):
        return self.is_set(PROPER_PAIR)

    @proper_pair.setter
    def proper_pair(self, value: bool):
        if value:
            self.set_bits(PROPER_PAIR)
        else:
            self.clear_bits(PROPER_PAIR)

    @property
    def unmapped(self):
        return self.is_set(QUERY_UNMAPPED)

    @unmapped.setter
    def unmapped(self, value: bool):
        if value:
            self.set_bits(QUERY_UNMAPPED)
        else:
            self.clear_bits(QUERY_UNMAPPED)

    @property
    def mate_unmapped(self):
        return self.is_set(MATE_UNMAPPED)

    @mate_unmapped.setter
    def mate_unmapped(self, value: bool):
        if value:
            self.set_bits(MATE_UNMAPPED)
        else:
            self.clear_bits(MATE_UNMAPPED)

    @property
    def reverse(self):
        return self.is_set(IS_REVERSE)

    @reverse.setter
    def reverse(self, value: bool):
        if value:
            self.set_bits(IS_REVERSE)
        else:
            self.clear_bits(IS_REVERSE)

    @property
    def mate_reverse(self):
        return self.is_set(MATE_IS_REVERSE)

    @mate_reverse.setter
    def mate_reverse(self, value: bool):
        if value:
            self.set_bits(MATE_IS_REVERSE)
        else:
            self.clear_bits(MATE_IS_REVERSE)

    @property
    def first_in_pair(self):
        return self.is_set(FIRST_IN_PAIR)

    @first_in_pair.setter
    def first_in_pair(self, value: bool):
        if value:
            self.set_bits(FIRST_IN_PAIR)
        else:
            self.clear_bits(FIRST_IN_PAIR)

    @property
    def last_in_pair(self):
        return self.is_set(LAST_IN_PAIR)

    @last_in_pair.setter
    def last_in_pair(self, value: bool):
        if value:
            self.set_bits(LAST_IN_PAIR)
        else:
            self.clear_bits(LAST_IN_PAIR)

    @property
    def secondary_alignment(self):
        return self.is_set(NOT_PRIMARY)

    @secondary_alignment.setter
    def secondary_alignment(self, value: bool):
        if value:
            self.set_bits(NOT_PRIMARY)
        else:
            self.clear_bits(NOT_PRIMARY)

    @property
    def does_not_pass_filters(self):
        return self.is_set(QC_FAILED)

    @does_not_pass_filters.setter
    def does_not_pass_filters(self, value: bool):
        if value:
            self.set_bits(QC_FAILED)
        else:
            self.clear_bits(QC_FAILED)

    @property
    def duplicate(self):
        return self.is_set(DUPLICATE)

    @duplicate.setter
    def duplicate(self, value: bool):
        if value:
            self.set_bits(DUPLICATE)
        else:
            self.clear_bits(DUPLICATE)

    @property
    def supplementary_alignment(self):
        return self.is_set(SUPPLEMENTARY)

    @supplementary_alignment.setter
    def supplementary_alignment(self, value: bool):
        if value:
            self.set_bits(SUPPLEMENTARY)
        else:
            self.clear_bits(SUPPLEMENTARY)

    @property
    def mate(self) -> (tuple):
        """
        Get the mate chromosome and postion.

        :return: the chromosome, position of the mate.
        """
        rchrom = self.mate_chromosome
        if self.mate_chromosome == "=":
            rchrom = self.chromosome
        return rchrom, self.mate_chromosome

    @mate.setter
    def mate(self, other):
        """
        Set the mate of the current alignment.

        :param other: the alignment to set as the mate
        """
        # set the position
        chrom = "="
        position = other.position
        tlen = 0

        # if the other read is unmapped
        if other.unmapped:
            chrom = "*"
            position = 0
            tlen = 0
        elif self.unmapped:
            chrom = other.chromosome
            position = other.position
            tlen = 0

        # if other read maps to a different chromosome
        elif self.chromosome != other.chromosome:
            chrom = other.chromosome
            tlen = 0
        else:
            # get in the insert size
            if other.position > self.position:
                lastpos = 0
                for cig in other.cigar_regions():
                    lastpos = cig[2]
                tlen = lastpos - self.position
            else:
                lastpos = 0
                for cig in self.cigar_regions():
                    lastpos = cig[2]
                tlen = lastpos - other.position

            # correct the tlen
            if self.position > other.position:
                tlen *= -1
            elif self.position == other.position and self.last_in_pair:
                tlen *= -1

        # set the mate pair chromosome, position and tlen fields
        self.mate_chromosome = chrom
        self.mate_position = position
        self.tlen = tlen

    def cigar_regions(self):
        """
        Interpret the cigar regions as separate regions.

        :yield: a tuple with the chromosome, start, end and operation
        """
        delta = 0
        cdelta = 0
        for nop, operation in cigar_operations(self.cigar):
            if operation in CIGAR_OPERATIONS_ON_REFERENCE:
                cdelta = delta
                delta += nop
            yield (
                self.chromosome,
                self.position + cdelta,
                self.position + delta,
                operation)
        return

    def get_tag(self, code: str="XX") -> (tuple):
        """
        Get a tag with a specific code.

        :param code: the 2 letter code of the tag
        :return: the tag with the specified code
        """
        for tag in self.tags:
            if tag[0] == code:
                return tag
        return None

    def get_tags(self, codes: list) -> (list):
        """Get a set of tags."""
        retval = [None] * len(codes)
        for tag in self.tags:
            if tag[0] in codes:
                retval[codes.index(tag[0])] = tag
        return retval

    @classmethod
    def format_tag(cls, tag: list) -> (str):
        """Format a tag."""
        return tag[0] + ":" + tag[1] + ":" + str(tag[2])

    def __repr__(self) -> (str):
        """
        Print the alignment in SAM format.

        :return: A SAM formatted string
        """
        tags = "\t".join([self.format_tag(tag) for tag in self.tags])
        line = "{name}\t{flag}\t{chromosome}\t{position}" + \
                  "\t{mapq}\t{cigar}\t{mate_chr}\t{mate_pos}" + \
                  "\t{tlen}\t{sequence}\t{quality}\t{tags}"
        return line.format(
            name=self.name,
            flag=self.flag,
            chromosome=self.chromosome,
            position=self.position,
            mapq=self.mapping_quality,
            cigar=self.cigar,
            mate_chr=self.mate_chromosome,
            mate_pos=self.mate_position,
            tlen=self.tlen,
            sequence=self.sequence,
            quality=self.quality,
            tags=tags)


def next_tag(tags: list) -> (list):
    """
    Get the next tag.

    :param tags: a list of strings representing tags
    :yield: the next tag
    """
    for tag in tags:
        tmp = tag.split(":")
        if len(tmp) == 3:
            yield (tmp[0], tmp[1], tmp[2])


def from_tokens(tokens: list):
    """
    Parse a SAM alignment from a list of tokens.

    :param tokens: a list with tokens
    :param sep: the field separator
    :return: a SAM alignment
    """
    tokens[1] = int(tokens[1])
    tokens[3] = int(tokens[3])
    tokens[4] = int(tokens[4])
    tokens[7] = int(tokens[7])
    tokens[8] = int(tokens[8])
    tags = [tag for tag in next_tag(tokens[11:])]
    return Alignment(*tokens[:11], tags)


def from_string(line: str, sep: str="\t"):
    """
    Parse a SAM alignment from a string.

    :param line: the line to parse
    :param sep: the field separator
    :return: a SAM alignment
    """
    return from_tokens(line.rstrip().split(sep))


class Reader:
    """Parse a SAM file."""

    def __init__(self, stream: typing.TextIO):
        """Initialize a new SAM reader."""
        self.header = []
        self.stream = stream
        self.lastline = None
        for line in self.stream:
            if not line:
                continue
            if line.startswith("@"):
                self.header.append(line.rstrip())
            else:
                self.lastline = line
                break

    def readgroups(self):
        """Get the readgroup to sample map."""
        matcher = re.compile("^.*ID:([^\t]+).*SM:([^\t]+).*$")
        samples = {}
        for line in self.header:
            if line.startswith("@RG"):
                regm = matcher.search(line)
                if regm:
                    samples[regm.group(1)] = regm.group(2)
        return samples

    def samples(self):
        """Get the samples."""
        return set(self.readgroups().values())

    def __iter__(self):
        """Treat this object as a generator."""
        if self.lastline:
            yield from_string(self.lastline)
            self.lastline = ""
        for line in self.stream:
            if not line:
                continue
            yield from_string(line)


class Writer:
    """Write a SAM file."""

    def __init__(self, stream, header):
        """Initialize the SAM writer."""
        self.stream = stream
        self.header = header
        self.write_header()

    def write_header(self):
        """Write the header."""
        if self.header is not None:
            for line in self.header:
                self.stream.write("{line}\n".format(line=line))

    def write(self, aln):
        """Write a new alignment to the SAM file."""
        self.stream.write("{aln}\n".format(aln=repr(aln)))


def quality_to_score(qchar, offset=32):
    """Convert the quality to a score."""
    return ord(qchar) - offset


def score_to_quality(qval, offset=32, maxval=126):
    """Convert a score to quality value."""
    cval = int(qval) + offset
    if cval > maxval:
        cval = maxval
    return chr(cval)


def encode_rle(lst):
    """Run length encode a list."""
    cur = []
    for val in lst:
        if cur:
            if cur[0] != val:
                yield len(cur), cur[0]
                cur = []
        cur.append(val)
    if cur:
        yield len(cur), cur[0]
