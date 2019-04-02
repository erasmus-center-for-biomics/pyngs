"""
A module to objectify BED parsing
"""


class Region:

    __slots__ = [
        "chromosome",
        "start",
        "end",
        "content"]

    def __init__(self, chromosome, start, end, **content):
        """Initialize a region element."""
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.content = content

    def has_overlap(self, other):
        """
        Do 2 regions have overlap.

        :param self: the first region
        :param other: the second region
        :return: True or False
        """
        if self.chromosome != other.chromosome:
            return False
        if self.start > other.end:
            return False
        if self.end < other.start:
            return False
        return True

    def overlap(self, other):
        """Get the number of overlapping bases.

        :param self: the first region
        :param other: the second region
        :return: the number of overlapping bases
        """
        if self.chromosome != other.chromosome:
            return 0
        rstart = max(self.start, other.start)
        rend = min(self.end, other.end)
        return 0 if rend < rstart else rend - rstart

    def contained_in(self, other):
        """Is self contained in another region.

        :param self: the first region
        :param other: the second region
        :return: True or False
        """
        if self.chromosome != other.chromosome:
            return False
        if self.start >= other.start and self.end <= other.end:
            return True
        return False


def region_parser(ftype="BED3"):
    """Get a region parser for regions ."""

    def bed3(tokens):
        return Region(
            tokens[0],
            tokens[1],
            tokens[2]
            ), tokens[3:]

    def bed4(tokens):
        return Region(
            tokens[0],
            tokens[1],
            tokens[2],
            comment=tokens[3]
            ), tokens[4:]

    def bedgraph(tokens):
        value = float(tokens[3]) if tokens[3] is not "." else None
        return Region(
            tokens[0],
            tokens[1],
            tokens[2],
            value=value
            ), tokens[4:]

    def bed6(tokens):
        return Region(
            tokens[0],
            tokens[1],
            tokens[2],
            comment=tokens[3],
            score=tokens[4],
            strand=tokens[5],
            ), tokens[6:]

    def gff(tokens):
        score = None if tokens[5] == "." else float(tokens[5])
        frame = None if tokens[7] == "." else int(tokens[7])
        return Region(
            tokens[0],
            tokens[3],
            tokens[4],
            score=score,
            frame=frame,
            strand=tokens[6],
            comment=tokens[8]), tokens[9:]


    if ftype.upper() == "BED3":
        return bed3
    elif ftype.upper() == "BED4":
        return bed4
    elif ftype.upper() == "BEDGRAPH":
        return bedgraph
    elif ftype.upper() == "BED6":
        return bed6
    elif ftype.upper() == "GFF":
        return gff
    return None