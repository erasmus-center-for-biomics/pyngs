"""This module contains functions to parse and compare alignments."""

import sys
import re
import pyngs.interval
from .cigar import CIGAR_OPERATIONS, CIGAR_OPERATIONS_ON_REFERENCE

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


class SAMAlignment(object):
    """A class to represent an alignment."""

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

    def is_set(self, bit):
        """Check whetehr the bit set in the flag."""
        return self.flag | bit == self.flag

    def is_paired(self):
        return self.is_set(IS_PAIRED)

    def is_proper_pair(self):
        return self.is_set(PROPER_PAIR)

    def is_unmapped(self):
        return self.is_set(QUERY_UNMAPPED)

    def is_mate_unmapped(self):
        return self.is_set(MATE_UNMAPPED)

    def is_reverse(self):
        return self.is_set(IS_REVERSE)

    def is_mate_reverse(self):
        return self.is_set(MATE_IS_REVERSE)

    def is_first_in_pair(self):
        return self.is_set(FIRST_IN_PAIR)

    def is_last_in_pair(self):
        return self.is_set(LAST_IN_PAIR)

    def is_secondary_alignment(self):
        return self.is_set(NOT_PRIMARY)

    def does_not_pass_filters(self):
        return self.is_set(QC_FAILED)

    def is_duplicate(self):
        return self.is_set(DUPLICATE)

    def is_supplementary_alignment(self):
        return self.is_set(SUPPLEMENTARY)

    def cigar_iterator(self):
        """Iterate over the elements from the cigar."""
        sl_start = 0
        sl_end = 0
        for idx in range(len(self.cigar)):
            if self.cigar[idx] in CIGAR_OPERATIONS:
                nop = int(self.cigar[sl_start:(sl_end+1)])
                sl_start = idx + 1
                sl_end = idx + 1
                yield nop, self.cigar[idx]
            else:
                sl_end = idx
        return

    def cigar_regions(self):
        """Interpret the cigar regions as separate regions."""
        delta = 0
        cdelta = 0
        for nop, operation in self.cigar_iterator():
            if operation in CIGAR_OPERATIONS_ON_REFERENCE:
                cdelta = delta
                delta += nop
            yield pyngs.interval.DataInterval(
                self.chromosome,
                self.position + cdelta,
                self.position + delta,
                operation)
        return

    def get_tag(self, code):
        """Get a tag with a specific code."""
        for tag in self.tags:
            if tag[0] == code:
                return tag
        return None

    @classmethod
    def format_tag(cls, tag):
        """Format a tag."""
        return tag[0] + ":" + tag[1] + ":" + str(tag[2])

    def __repr__(self):
        tags = "\t".join([self.format_tag(tag) for tag in self.tags])
        return "{name}\t{flag}\t{chromosome}\t{position}\t{mapq}\t{cigar}\t{mate_chr}\t{mate_pos}\t{tlen}\t{sequence}\t{quality}\t{tags}".format(
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


class SAMParser(object):
    """Create alignments from strings."""

    def __init__(self,
                 handle=sys.stdin,
                 chromosome_list=None,
                 allow_append=True):
        """Initialize a new alignment factory."""
        if chromosome_list is None:
            self.chromosome_list = []
        else:
            self.chromosome_list = chromosome_list
        self.header = []
        self.allow_append = allow_append
        self.handle = handle
        self.count = 0

    def sample_map(self):
        """Get the sample map."""
        matcher = re.compile("^.*ID:([^\t]+).*SM:([^\t]+).*$")
        samples = {}
        for line in self.header:
            if line.startswith("@RG"):
                regm = matcher.search(line)
                if regm:
                    samples[regm.group(1)] = regm.group(2)
        return samples

    def from_string(self, sam_string="", sep="\t"):
        """Interpret a SAM string."""
        fields = sam_string.rstrip().split(sep)
        readname = fields[0]
        flag = int(fields[1])
        chromosome_name = fields[2]
        position = int(fields[3])
        mqual = int(fields[4])
        cigar = fields[5]
        mate_chromosome = fields[6]
        mate_position = int(fields[7])
        tlen = int(fields[8])
        sequence = fields[9]
        quality = fields[10]
        tags = self.process_tags(fields[11:])

        # return the interpreted Alignment fields
        return SAMAlignment(
            readname, flag, chromosome_name, position, mqual, cigar,
            mate_chromosome, mate_position, tlen,
            sequence, quality, tags)

    def process_tags(self, tags):
        """Process the tags."""
        rval = []
        for tag in tags:
            tmp = tag.split(":")
            if len(tmp) == 3:
                rval.append((tmp[0], tmp[1], tmp[2]))
        return rval

    def parse_seq_header(self, line):
        """Add a chromosome to the chromosome list."""
        if line.startswith("@SQ"):
            fields = line.split("\t")
            namefld = [fld for fld in fields if fld.startswith("SN:")]
            if not namefld:
                raise ValueError(
                    "Improperly formatted header: {0}".format(line))
            name = namefld[0].replace("SN:", "")
            if name not in self.chromosome_list:
                self.chromosome_list.append(name)

    def __iter__(self):
        """Mark this object as an iterator."""
        return self

    def __next__(self):
        """Get the next entry in the file."""
        for line in self.handle:
            self.count += 1
            line = line.rstrip()
            if not line:
                continue
            if line[0] == "@":
                self.header.append(line)
                self.parse_seq_header(line)
                continue
            return self.from_string(line)
        raise StopIteration

    def next(self):
        """Be compatible with python 2."""
        return self.__next__()

    def chromosome_index(self, name):
        """Get the chromosome index."""
        try:
            return self.chromosome_list.index(name)
        except ValueError:
            return -1


class SAMWriter(object):
    """A class to write SAM files."""

    def __init__(self, outstream=sys.stdout, header=None):
        """Initialize a new SAMWriter."""
        self.outstream = outstream
        self.header = header
        self.header_written = False

    def set_header(self, header=None):
        """Add the header to the object."""
        self.header = header

    def write_header(self):
        """Write the header."""
        if self.header is not None:
            for line in self.header:
                self.outstream.write("{line}\n".format(line=line))
            self.header_written = True

    def write(self, aln):
        """Write a new alignment to the SAM file."""
        # write the header if available
        if not self.header_written and self.header is not None:
            self.write_header()

        self.outstream.write("{aln}\n".format(aln=repr(aln)))
