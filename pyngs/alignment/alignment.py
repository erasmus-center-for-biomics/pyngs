"""This module contains functions to parse and compare alignments."""

import sys
import BioInterval


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


class Alignment(object):
    """A class to represent an alignment."""

    __slots__ = [
        "name", "flag", "chridx", "position",
        "mapping_quality", "cigar", "mate_chromosome",
        "mate_position", "tlen", "sequence", "quality",
        "tags"]

    cigar_operations = ("M", "I", "D", "N", "S", "H", "P", "=", "X")
    cigar_operations_on_reference = ("M", "D", "N", "P", "=", "X")

    def __init__(self, *args):
        """Create a new alignment."""
        self.name = args[0]
        self.flag = args[1]
        self.chridx = args[2]
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

    def is_first_in_template(self):
        return self.is_set(FIRST_IN_PAIR)

    def is_last_in_template(self):
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
            if self.cigar[idx] in Alignment.cigar_operations:
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
            if operation in self.cigar_operations_on_reference:
                cdelta = delta
                delta += nop
            yield BioInterval.interval.DataInterval(
                self.chridx,
                self.position + cdelta,
                self.position + delta,
                operation)
        return


class AlignmentFactory(object):
    """Create alignments from strings."""

    def __init__(self, chromosome_list=None, allow_append=True):
        """Initialize a new alignment factory."""
        assert chromosome_list is not None
        self.chromosome_list = chromosome_list
        self.header = []
        self.allow_append = allow_append

    def from_string(self, sam_string=""):
        """Interpret a SAM string."""
        fields = sam_string.rstrip().split("\t")
        readname = fields[0]
        flag = int(fields[1])
        chromosome_name = fields[2]
        position = int(fields[3]) - 1
        mqual = int(fields[4])
        cigar = fields[5]
        mate_chromosome = fields[6]
        mate_position = int(fields[7]) - 1
        tlen = int(fields[8])
        sequence = fields[9]
        quality = fields[10]
        tags = self.process_tags(fields[11:])

        # add the chromosome to the chromosome list
        chridx = len(self.chromosome_list)
        try:
            chridx = self.chromosome_list.index(chromosome_name)
        except ValueError:
            if self.allow_append:
                self.chromosome_list.append(chromosome_name)
            else:
                raise ValueError("Unknown chromosome %s" % chromosome_name)

        # return the interpreted Alignment fields
        return Alignment(readname, flag, chridx, position, mqual, cigar,
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
            name = fields[1].replace("SN:", "")
            if name not in self.chromosome_list:
                self.chromosome_list.append(name)

    def iterator(self, handle=sys.stdin):
        """Iterate over the lines of a stream."""
        while True:
            line = next(handle)
            if len(line) == 0:
                continue
            if line[0] == '@':
                self.header.append(line)
                self.parse_seq_header(line)
                continue
            yield self.from_string(line)
