"""Parsers for various BED formats."""


def parse_bed3(tokens: list):
    """
    Parse a BED3 line.

    :param tokens: the tokens to parse
    :return: a string, int, int representing the
            chromosome, start and end of the
            region
    """
    # skip track lines
    if tokens[0].startswith("track"):
        return None
    # check the field length
    if len(tokens) != 3:
        return None
    # return the region
    return tokens[0], int(tokens[1]), int(tokens[2])


def parse_bed4(tokens: list):
    """
    Parse a BED4 line.

    :param tokens: the tokens to parse
    :return: a string, int, int, str representing the
            chromosome, start, end of the
            region and the comment
    """
    # skip track lines
    if tokens[0].startswith("track"):
        return None
    # check the field length
    if len(tokens) != 4:
        return None
    # return the region
    return tokens[0], int(tokens[1]), int(tokens[2]), tokens[3]


def parse_bedgraph(tokens: list):
    """
    Parse a BEDGRAPH line.

    :param tokens: the tokens to parse
    :return: a string, int, int, float representing the
            chromosome, start, end and score of the
            region
    """
    if tokens[0].startswith("track"):
        return None
    # check the field length
    if len(tokens) != 4:
        return None
    # yield the region
    return tokens[0], int(tokens[1]), int(tokens[2]), float(tokens[3])


def parse_bed5(tokens: list):
    """
    Parse a BED6 line.

    :param tokens: the tokens to parse
    :yield: a string, int, int, str, int, str
            representing the chromosome, start,
            end, name, score and strand
    """
    # skip track lines
    if tokens[0].startswith("track"):
        return None
    # check the field length
    if len(tokens) != 6:
        return None
    # yield the region
    return tokens[0], int(tokens[1]), int(tokens[2]),\
        tokens[3], int(tokens[4]), tokens[5]


def parse_bedpe(tokens: list):
    """
    Parse a BEDPE line.

    :param tokens: the tokens to parse
    :yield: a string, int, int, string, int, int, str, float
            representing the locations of the first and
            second end (chromosome, start, end), a comment,
            and a score.
    """
    # skip track lines
    if tokens[0].startswith("track"):
        return None
    # check the field length
    if len(tokens) != 8:
        return None
    # yield the region
    return tokens[0], int(tokens[1]), int(tokens[2]), \
        tokens[3], int(tokens[4]), int(tokens[5]), \
        tokens[6], float(tokens[7])


def parse_gff(tokens: list):
    """
    Parse a GFF line.

    :param instream: the input stream
    :return: a str, int, int, str, str, str,
            float, str, int, and str representing
            the 0) chromosome, 1) start, 2) end, 3) source,
            4) type, 5) score, 6) strand, 7) frame, and
            8) attributes of the entry
    """
    # skip track lines
    if tokens[0].startswith("#"):
        return None
    if len(tokens) == 8:
        tokens.append("")
    elif len(tokens) != 9:
        return None
    score = None if tokens[5] == "." else float(tokens[5])
    frame = None if tokens[7] == "." else int(tokens[7])

    # yield the region
    return tokens[0], int(tokens[3]), int(tokens[4]),\
        tokens[1], tokens[2], score, tokens[6], frame,\
        tokens[8]


def parse_gtf_attributes(attributes: str):
    """
    Parse the attributes from a GTF file.

    :param attributes: the attribute string to parse.
    :return: a dict with the parsed values
    """
    retval = {}
    for token in attributes.split(";"):
        pair = token.strip().split(" ")
        if len(pair) == 2:
            retval[pair[0]] = pair[1].strip('"')
    return retval


def parse_gff_attributes(attributes: str):
    """
    Parse the attributes from a GFF file.

    :param attributes: the attribute string to parse.
    :return: a dict with the parsed values
    """
    retval = {}
    for token in attributes.split(";"):
        pair = token.strip().split("=")
        if len(pair) == 2:
            retval[pair[0]] = pair[1].strip('"')
    return retval
