"""A package to work with BED files."""

import typing
from .parsers import parse_gff
from .parsers import parse_bed3, parse_bed4, parse_bed5
from .parsers import parse_bedgraph, parse_bedpe


NFIELDS = {
    "bed3": 3,
    "bed4": 4,
    "bed5": 6,
    "bedpe": 8,
    "bedgraph": 4,
    "gff": 9
}


def reader(instream: typing.TextIO, ftype: str="bed3", sep: str="\t"):
    """Parse the lines from a bed file.

    :param instream: the input stream
    :param format: the format of the input file.
                   The choices here are bed3, bed5, bedpe,
                   bedgraph, gtf and gff
    :param sep: the token separator
    :exception ValueError: raised when the format is
                           not recognized
    :returns: a generator that will yield the
              chromosome, start, end and other fields
              based on the format
    """
    # get the line parser
    parser = None
    if ftype == "bed3":
        parser = parse_bed3
    if ftype == "bed4":
        parser = parse_bed4
    if ftype == "bed5":
        parser = parse_bed5
    if ftype == "bedgraph":
        parser = parse_bedgraph
    if ftype == "bedpe":
        parser = parse_bedpe
    if ftype == "gtf":
        parser = parse_gff
    if ftype == "gff":
        parser = parse_gff

    # raise a ValueError
    if not parser:
        raise ValueError(
            "format {format} could not be identified".format(
                format=format))

    # process the file line by line
    for line in instream:
        tokens = line.rstrip().split(sep)
        result = parser(tokens)
        if not result:
            continue
        yield result


def intersect_reader(instream: typing.TextIO, sep: str="\t",
                     parser_a: typing.Callable=None, fields_a: int=0,
                     parser_b: typing.Callable=None, fields_b: int=0):
    """Parse the output from bedtools intersect.

    :param instream: the input stream
    :param sep: the token separator
    :param parser_a: the token parser for the first part
    :param fields_a: the number of columns of A
    :param parser_b: the token parser for the second part
    :param fields_b: the number of columns of format B or 0 if variable
    :yield: part 1, part 2 and the rest of the line
    """
    # check the input data
    assert parser_a is not None
    assert parser_b is not None
    assert fields_a > 0
    if fields_b > 0:
        assert fields_a < fields_b

    # process the file line by line
    for line in instream:
        # get the tokens from the line
        tokens = line.rstrip().split(sep)

        # get the part_a
        part_a = parser_a(tokens[:fields_a])

        # get part_b
        part_b = None
        toparse = tokens[fields_a:fields_b] \
            if fields_b > 0 else \
            tokens[fields_a:]
        part_b = parser_b(toparse)

        # get any additional fields
        rest = None
        if fields_b > len(tokens):
            rest = tokens[fields_b:]

        # yield the parts
        yield part_a, part_b, rest


def width(region):
    """
    Get the width of a region.

    :param region: a tuple representing the region
                   with chromosome, start, end
                   as the first 3 columns
    :return: an int
    """
    return region[2] - region[1]


def has_overlap(frst, scnd):
    """
    Do 2 regions have overlap.

    :param frst: a tuple representing the first region
                 with chromosome, start, end as the first
                 3 columns
    :param scnd: a tuple representing the second region
                 with chromosome, start, end as the first
                 3 columns
    :return: True or False
    """
    if frst[0] != scnd[0]:
        return False
    if frst[1] > scnd[2]:
        return False
    if frst[2] < scnd[1]:
        return False
    return True


def bases_overlap(frst, scnd):
    """
    Get the number of overlapping bases.

    :param frst: a tuple representing the first region
                 with chromosome, start, end as the first
                 3 columns
    :param scnd: a tuple representing the second region
                 with chromosome, start, end as the first
                 3 columns
    :return: the number of overlapping bases
    """
    if frst[0] != scnd[0]:
        return 0
    rstart = max(frst[1], scnd[1])
    rend = min(frst[2], scnd[2])
    if rend < rstart:
        return 0
    return rend - rstart


def is_contained_in(frst, scnd):
    """
    Is the first region contained in the second.

    :param frst: a tuple representing the first region
                 with chromosome, start, end as the first
                 3 columns
    :param scnd: a tuple representing the second region
                 with chromosome, start, end as the first
                 3 columns
    :return: True or False
    """
    if frst[0] != scnd[0]:
        return False
    if frst[1] >= scnd[1] and frst[2] <= scnd[2]:
        return True
    return False
