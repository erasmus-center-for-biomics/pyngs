
import sys
import gzip
from pyngs import sam


class TagFilter(object):
    """An object to filter on tag value."""

    def __init__(self):
        """Initialize a tag filter."""
        self.tag = None
        self.value = None
        self.compare = None
        self.discard = False

    def process(self, reader, writer):
        """Process the SAM file."""
        for aln in reader:
            # get the tag
            tag = aln.get_tag(self.tag)

            # if the tag is absent
            if tag is None:
                if not self.discard:
                    writer.write(aln)
                continue

            # otherwise compare the tag to the value
            # and write the alignment if the comparison
            # evaluates to true
            if self.compare(tag, self.value):
                writer.write(aln)

    @classmethod
    def equals(cls, tag, value):
        """Is the tag value equal to value."""
        try:
            if tag[1] == ("Z", "A", "H", "B"):
                return tag[2] == value
            elif tag[1] == "i":
                return int(tag[2]) == int(value)
            elif tag[1] == "f":
                return float(tag[2]) == float(value)
        except ValueError:
            return False

    @classmethod
    def greater(cls, tag, value):
        """Is the tag value greater than value."""
        try:
            if tag[1] == ("Z", "A", "H", "B"):
                return tag[2] > value
            elif tag[1] == "i":
                return int(tag[2]) > int(value)
            elif tag[1] == "f":
                return float(tag[2]) > float(value)
        except ValueError:
            return False

    @classmethod
    def less(cls, tag, value):
        """Is the tag value less than value."""
        try:
            if tag[1] == ("Z", "A", "H", "B"):
                return tag[2] > value
            elif tag[1] == "i":
                return int(tag[2]) < int(value)
            elif tag[1] == "f":
                return float(tag[2]) < float(value)
        except ValueError:
            return False


def filter_by_samtag(args):
    """Filter by samtag."""
    # open the files
    instream = sys.stdin
    if args.sam != "stdin":
        if args.sam.endswith(".gz"):
            instream = gzip.open(args.sam, "rt")
        else:
            instream = open(args.sam, "rt")

    outstream = sys.stdout
    if args.out != "stdout":
        if args.out.endswith(".gz"):
            outstream = gzip.open(args.out, "wt")
        else:
            outstream = open(args.out, "wt")

    # prepare the filter
    tagfilter = TagFilter()
    tagfilter.tag = args.tag
    tagfilter.value = args.value
    if args.discard:
        tagfilter.discard = True

    # set the compare function
    if args.type == "equals":
        tagfilter.compare = tagfilter.equals
    elif args.type == "greater":
        tagfilter.compare = tagfilter.greater
    elif args.type == "less":
        tagfilter.compare = tagfilter.less

    # process the file
    reader = sam.Reader(instream)
    writer = sam.Writer(outstream, reader.header)
    tagfilter.process(reader, writer)

    # close file handles when we are done
    if not outstream.closed and outstream != sys.stdout:
        outstream.close()
    if not instream.closed and instream != sys.stdin:
        instream.close()
