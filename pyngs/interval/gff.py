"""GFF reader classes."""

import sys
import BioInterval


def gff_file(gffstream=sys.stdin):
    """
    Read the contents of a GFF file on by one.

    Arguments:
        gffstream - a stream of the content of a GFF file

    Yields:
        a zipped list with gff content
 
    """
    header = ["chromosome", "source", "type",
              "start", "end", "score", "strand",
              "frame", "comment"]
    while True:
        line = next(gffstream).rstrip()

        # skip empty lines and comments
        if len(line) == 0:
            continue
        if line[0] == "#":
            continue

        # return the data fields annotated
        # with the header
        fields = line.split("\t")
        fields[3] = int(fields[3])
        fields[4] = int(fields[4])
        yield zip(header, fields)
    return


def read_gff_entries(gffstream=sys.stdin, typeids=("exon")):
    """
    Read the entries of type typeid from a GFF file.

    Arguments:
        gffstream - a stream of the content of a GFF file
        typeids - a tuple with typeids to obtain

    Returns:
        get a list of entries of type typeids

    """
    entries = []
    for entry in gff_file(gffstream):
        if entry[2][1] in typeids:
            entries.append(entry)
    return entries


def gff_parse_comment(commentfield):
    """
    Parse the comment field from a GFF file.

    Arguments:
        commentfield - the commentfield to parse

    Returns:
        a list with tuples of key value pairs

    """
    retval = []
    fieldpairs = commentfield.split(";")
    for entry in fieldpairs:
        pairs = entry.lstrip().split(" ")
        if len(pairs) == 2:
            pairs[1] = pairs[1].strip("\"")
            retval.append((pairs[0], pairs[1]))

    # return the list of identifiers
    return retval


def gff_entry_to_interval(entry, chromosomes=None, allow_add=True):
    """
    Convert a GFF entry to an interval.

    Arguments:
        entry - the result from the gfffile iterator
        chromosomes - a list with chromosomes

    Returns:
        An BioInterval.Interval from the entry

    """
    try:
        chridx = chromosomes.index(entry[0][1])
    except ValueError as err:
        if allow_add:
            chromosomes.append(entry[0][1])
            chridx = len(chromosomes) - 1
        else:
            raise err
    return BioInterval.Interval(
        chridx,
        entry[3][1],
        entry[4][1])
