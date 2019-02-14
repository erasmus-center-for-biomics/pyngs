import gzip
import sys

from pyngs import sam


def to_bed(instream, outstream, operations, tags):
    """Convert the CIGAR strings to a BED file."""
    bedline = "{chromosome}\t{start}\t{end}\t{comment}\t{score}\t{strand}\n"
    reader = sam.Reader(instream)
    readgroups = reader.readgroups()

    for alignment in reader:

        # get the strand
        strand = "-" if alignment.reverse else "+"

        taglst = []
        if tags:
            tagfnd = alignment.get_tags(tags)
            for idx, tagname in enumerate(tags):
                if tagname == "__sample__":
                    rgrp = alignment.get_tag("RG")
                    if rgrp:
                        taglst.append("__sample__={0}".format(readgroups[rgrp[2]]))
                    else:
                        taglst.append("__sample__=None")
                else:
                    tagres = tagfnd[idx] if tagfnd[idx] else "None"
                    taglst.append("{tag}={result}".format(
                        tag=tagname,
                        result=tagres[2]))
        tagstr = ";".join(taglst)

        # iterate over the cigar operations
        for cigar in alignment.cigar_regions():

            # skip cigar operations that we are not interested in
            if cigar[3] not in operations:
                continue
            comment = "READNAME={name};COP={op};{tags}".format(
                name=alignment.name,
                op=cigar[3],
                tags=tagstr)
            outstream.write(
                bedline.format(
                    chromosome=cigar[0],
                    start=cigar[1],
                    end=cigar[2],
                    comment=comment,
                    score=alignment.mapping_quality,
                    strand=strand))


def cigar_to_bed(args):
    """Convert CIGAR entries to a BED file."""
    instream = sys.stdin
    if args.sam != "stdin":
        if args.sam.endswith(".gz"):
            instream = gzip.open(args.sam, "rt")
        else:
            instream = open(args.sam, "rt")

    outstream = sys.stdout
    if args.bed != "stdout":
        if args.bed.endswith(".gz"):
            outstream = gzip.open(args.bed, "wt")
        else:
            outstream = open(args.bed, "wt")

    # write the BED entries
    to_bed(instream, outstream, args.operations, args.tags)

    # close the streams
    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()
