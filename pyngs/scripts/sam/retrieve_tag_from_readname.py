import sys
import gzip
import argparse
import pyngs.sam as sam


def process_tags(instream, outstream, tags, samtags, sep=":"):
    """Add tags to the alignment."""
    reader = sam.Reader(instream)
    writer = sam.Writer(outstream, reader.header)

    for alignment in reader:

        # get the values from the readname
        toadd = []
        toremove = []
        parts = alignment.name.split(sep)
        for pidx, part in enumerate(parts):
            if "=" in part:
                pair = part.split("=", 1)
                if pair[0] in tags:
                    tidx = tags.index(pair[0])
                    toadd.append((samtags[tidx], pair[1]))
                    toremove.append(pidx)
        # add the new tags to the alignment 
        previous = [t[0] for t in alignment.tags] 
        for tag in toadd:
            if tag[0] not in previous:
                alignment.tags.append((tag[0], "Z", tag[1]))

        # construct a new readname with the information removed
        newname = sep.join([p for i, p in enumerate(parts) if i not in toremove])
        alignment.name = newname

        # write the alignment to the output 
        writer.write(alignment)


def retrieve_tag_from_readname(args: argparse.Namespace) -> None:
    """Retriece tags from the readnames and add them to the alignments."""
    instream = sys.stdin
    if args.input != "stdin":
        if args.input.endswith(".gz"):
            instream = gzip.open(args.input, "rt")
        else:
            instream = open(args.input, "rt")

    outstream = sys.stdout
    if args.output != "stdout":
        if args.output.endswith(".gz"):
            outstream = gzip.open(args.output, "wt")
        else:
            outstream = open(args.output, "wt")
    
    if len(args.tags) != len(args.samtags):
        raise ValueError(
            "The length of the tags list ({0}) differs from the sam-tags ({1})".format(
                len(args.tags), len(args.samtags)))
    
    # process the tags from the 
    process_tags(instream, outstream, args.tags, args.samtags)

    # close the streams
    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":

    sparser = argparse.ArgumentParser(
        prog="",
        description="""Extract information from the read names 
        and add them as tags in a SAM file.""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, default="stdin",
        help="""The input SAM file.""")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, default="stdout",
        help="The output SAM file.")
    sparser.add_argument(
        "-t", "--tags", dest="tags",
        type=str, nargs="+",
        help="""The tags to extract from the read name""")
    sparser.add_argument(
        "-s", "--sam-tag", dest="samtags",
        type=str, nargs="+",
        help="""The tag to deposit in the SAM file""")

    sparser.set_defaults(func=retrieve_tag_from_readname)
    args = sparser.parse_args()
    args.func(args)

