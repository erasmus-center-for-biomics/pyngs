import argparse
import sys
import gzip
import pyngs.vcf as vcf


def extract(instream, outstream):
    """."""
    reader = vcf.Reader(instream)
    outline = "{chrom}\t{position}\t{sample}\t{A}\t{B}\n"
    outstream.write(outline.format(
        chrom="chrom",
        position="position",
        sample="sample",
        A="A",
        B="B"))

    # get the field parser for the genotype and alleles
    gtparser = reader.field("FORMAT", "GT")
    patparser = reader.field("FORMAT", "PAT")
    matparser = reader.field("FORMAT", "MAT")
    childparser = reader.field("FORMAT", "OFF")

    for variant in reader:

        try:
            gidx = variant.format.index("GT")
            pidx = variant.format.index("PAT")
            midx = variant.format.index("MAT")
            cidx = variant.format.index("OFF")
        except ValueError:
            continue

        # get the alleles
        alleles = [variant.reference]
        alleles.extend(variant.alternate.split(","))

        #
        for sidx, sample in enumerate(variant.samples):

            # get the relevant fields
            gt = list(gtparser.interpret(sample[gidx]))[0]
            patallele = list(patparser.interpret(sample[pidx]))[0]
            matallele = list(matparser.interpret(sample[midx]))[0]
            childallele = list(childparser.interpret(sample[cidx]))[0]

            # A normal haplotyped child
            if patallele is not "." and matallele is not ".":
                outstream.write(outline.format(
                    chrom=variant.chrom,
                    position=variant.position,
                    sample=reader.samples[sidx],
                    A=patallele,
                    B=matallele))

            # extract haplotypes for a parent that
            # was used to genotype a child
            elif childallele is not "." and gt is not ".":
                allele_a = childallele.split("/")[1]
                # skip erroneous alleles
                if allele_a == "?" or allele_a == "E":
                    outstream.write(outline.format(
                        chrom=variant.chrom,
                        position=variant.position,
                        sample=reader.samples[sidx],
                        A=allele_a,
                        B=allele_a))
                    continue
                #
                genotype = [alleles[int(v)] for v in gt.split("/")]
                other = [g for g in genotype if g not in [allele_a]]
                allele_b = other[0] if other else allele_a

                outstream.write(outline.format(
                    chrom=variant.chrom,
                    position=variant.position,
                    sample=reader.samples[sidx],
                    A=allele_a,
                    B=allele_b))

def extract_haplotypes(args):
    """."""
    instream = sys.stdin
    if args.vcf != "stdin":
        if args.vcf.endswith(".gz"):
            instream = gzip.open(args.vcf, "rt")
        else:
            instream = open(args.vcf, "rt")

    outstream = sys.stdout
    if args.output != "stdout":
        if args.output.endswith(".gz"):
            outstream = gzip.open(args.output, "wt")
        else:
            outstream = open(args.output, "wt")

    #
    extract(instream, outstream)

    # close the input and output streams
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()


if __name__ == "__main__":

    # parse the commandline arguments
    parser = argparse.ArgumentParser(
        prog="extract_haplotypes",
        description=""".""")
    parser.add_argument(
        "-v", "--vcf", dest="vcf",
        type=str, nargs="?", default="stdin",
        help="The input VCF file.")
    parser.add_argument(
        "-o", "--output", dest="output",
        type=str, nargs="?", default="stdout",
        help="The output VCF file.")
    parser.set_defaults(func=extract_haplotypes)

    #
    args = parser.parse_args()
    args.func(args)
