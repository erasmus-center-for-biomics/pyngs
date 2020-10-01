import gzip
import sys
from operator import itemgetter
from pyngs import sam


def format_all_cigars(alignment, tags, operations):
    """Format all cigars for the alignment."""
    strand = "-" if alignment.reverse else "+"
    for cigar in alignment.cigar_regions():
        if cigar[3] not in operations:
            continue
        
        # prepare the comment
        comment = "READNAME={name};COP={op};{tags}".format(
            name=alignment.name,
            op=cigar[3],
            tags=";".join(tags))
        yield (cigar[0], cigar[1], cigar[2], comment, alignment.mapping_quality, strand)
    

def format_single_cigars(alignment, tags, operations):
    """."""
    strand = "-" if alignment.reverse else "+"
    comment = "READNAME={name};{tags}".format(
        name=alignment.name,
        tags=";".join(tags))
    cigars = [cig for cig in alignment.cigar_regions() if cig[3] in operations]
    cigars.sort(key=itemgetter(0, 1))
    yield (cigars[0][0], cigars[0][1], cigars[len(cigars) - 1][2], comment, alignment.mapping_quality, strand)


def to_bed(instream, outstream, operations, tags, merge_entries=False):
    """Convert the CIGAR strings to a BED file."""
    bedline = "{chromosome}\t{start}\t{end}\t{comment}\t{score}\t{strand}\n"
    reader = sam.Reader(instream)
    readgroups = reader.readgroups()

    # for each alignment in the reader
    for alignment in reader:

        # get the tags
        taglst = []
        if tags:
            tagfnd = alignment.get_tags(tags)
            for idx, tagname in enumerate(tags):
                if tagname == "__sample__":
                    rgrp = alignment.get_tag("RG")
                    if rgrp:
                        taglst.append(
                            "__sample__={0}".format(
                                readgroups[rgrp[2]]))
                    else:
                        taglst.append("__sample__=None")
                else:
                    tagres = tagfnd[idx] if tagfnd[idx] else (tagname, "", "None")
                    taglst.append("{tag}={result}".format(
                        tag=tagname,
                        result=tagres[2]))
        
        # print the cigars
        if not merge_entries:
            for cigar in format_all_cigars(alignment, taglst, operations):
                outstream.write(
                    bedline.format(
                        chromosome=cigar[0],
                        start=cigar[1],
                        end=cigar[2],
                        comment=cigar[3],
                        score=cigar[4],
                        strand=cigar[5]))
        else:
            for cigar in format_single_cigars(alignment, taglst, operations):
                outstream.write(
                    bedline.format(
                        chromosome=cigar[0],
                        start=cigar[1],
                        end=cigar[2],
                        comment=cigar[3],
                        score=cigar[4],
                        strand=cigar[5]))
            

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
    to_bed(
        instream,
        outstream, 
        args.operations, 
        args.tags, 
        merge_entries=args.merge_entries)

    # close the streams
    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()
