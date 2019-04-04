import sys
import gzip
import argparse
import pyngs.vcf as vcf
from operator import itemgetter


def count(lst):
    """Count the entries in a list."""
    lst.sort()
    score = 0
    pentry = None
    for entry in lst:
        if entry != pentry and score:
            yield pentry, score
            score = 0
        pentry = entry
        score += 1
    if score:
        yield pentry, score


def haplotype(calleles, palleles, malleles):
    """Haplotype a child with the information from its parents.

    return the paternal and maternal alleles in that order
    """

    # homozygous child
    if calleles[0] == calleles[1]:
        return calleles[0], calleles[1]
    # heterozygous child
    else:
        # filter all alleles that are not present in the child
        parental = [a for a in palleles if a in calleles]
        parental.extend([a for a in malleles if a in calleles])

        allele_count = [t for t in count(parental)]
        allele_count.sort(key=itemgetter(1))

        # error case where the child has different alleles than the parents
        if len(allele_count) <= 1:
            return "E", "E"
        # unsolvable case where child is Q/P, father is Q/P and mother is Q/P
        elif allele_count[0][1] == 2:
            return "?", "?"
        elif allele_count[0][1] == 1:
            # case child is Q/P, father is Q/P,
            # mother is P/P
            if allele_count[0][0] in palleles:
                return allele_count[0][0], allele_count[1][0]
            elif allele_count[0][0] in malleles:
                return allele_count[1][0], allele_count[0][0]
    # should never happen
    return None, None


def assign_haplotypes(instream, outstream, assignments, tag="GT"):
    """Assign the haplotypes to a VCF file."""
    reader = vcf.Reader(instream)

    # get the field parser for the genotype
    parser = reader.field("FORMAT", tag)

    # get the sample indexes for the child, father, mother group
    links = []
    for assignment in assignments:
        try:
            entry = (
                reader.samples.index(assignment[0]),
                reader.samples.index(assignment[1]),
                reader.samples.index(assignment[2]))
            links.append(entry)
        except ValueError:
            sys.stdout.write("Error: Could not find the samples for group {0}, skipping\n".format(assignment))

    # add the haplotype fields to the header
    header = reader.header
    header.append(vcf.Header("""FORMAT=<ID=MAT,Number=1,Type=String,Description="The maternal allele">"""))
    header.append(vcf.Header("""FORMAT=<ID=PAT,Number=1,Type=String,Description="The paternal allele">"""))
    header.append(vcf.Header("""FORMAT=<ID=OFF,Number=.,Type=String,Description="The allele transmited to the offspring">"""))

    # prepare the header
    writer = vcf.Writer(outstream, header, reader.samples)

    # foreach variant in the reader
    for variant in reader:
        # skip all
        # if variant.chrom != "chr11" or variant.position != "5265578":
        #     continue
        try:
            fidx = variant.format.index(tag)
        except ValueError:
            # could not find the index of tag
            # so just add the data
            writer.write(variant)
            continue

        # get the alleles
        alleles = [variant.reference]
        alleles.extend(variant.alternate.split(","))

        # prepare the new format fields
        paternal = [None] * len(variant.samples)
        maternal = [None] * len(variant.samples)
        offspring = [None] * len(variant.samples)
        modified = False
        # for each child parent group
        for (cidx, pidx, midx) in links:

            # get the genotypes for child, father and mother
            cgt = list(parser.interpret(variant.samples[cidx][fidx]))[0]
            pgt = list(parser.interpret(variant.samples[pidx][fidx]))[0]
            mgt = list(parser.interpret(variant.samples[midx][fidx]))[0]

            # skip unsolvable cases
            if cgt == "./." or pgt == "./." or mgt == "./.":
                continue

            # get the alleles for a child
            try:
                calleles = [alleles[int(v)] for v in cgt.split("/")]
                palleles = [alleles[int(v)] for v in pgt.split("/")]
                malleles = [alleles[int(v)] for v in mgt.split("/")]
            # skip variants which we cannot identify
            except ValueError:
                continue

            # check our sanity
            assert len(calleles) == 2
            assert len(palleles) == 2
            assert len(malleles) == 2

            # signal that we added format data
            modified = True

            # assign the haplotypes
            paternal[cidx], maternal[cidx] = haplotype(
                calleles, palleles, malleles)

            # add the offspring data
            tmp = "{0}/{1}".format(
                reader.samples[cidx],
                paternal[cidx]
            )
            if offspring[pidx] is None:
                offspring[pidx] = tmp
            else:
                offspring[pidx] += "," + tmp

            tmp = "{0}/{1}".format(
                reader.samples[cidx],
                maternal[cidx])
            if offspring[midx] is None:
                offspring[midx] = tmp
            else:
                offspring[midx] += "," + tmp

        # add the haplotype information
        if modified:
            variant.add_data("PAT", paternal)
            variant.add_data("MAT", maternal)
            variant.add_data("OFF", offspring)
        writer.write(variant)


def haplotype_vcf(args):
    """Haplotype a co-called VCF file."""
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

    # add the assignments
    assignments = []
    if args.relations:
        stream = open(args.relations, "rt")
        for line in stream:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            tokens = line.split("\t")
            assignments.append(
                (tokens[0], tokens[1], tokens[2]))
        stream.close()

    # add the commandline assignments
    if args.offspring and args.fathers and args.mothers:
        offspring = [v for v in args.offspring if v != ""]
        fathers = [v for v in args.fathers if v != ""]
        mothers = [v for v in args.mothers if v != ""]

        assert len(offspring) == len(fathers)
        assert len(offspring) == len(mothers)

        # group the offspring, parent pairs
        for idx, child in enumerate(offspring):
            assignments.append(
                (child, fathers[idx], mothers[idx]))

    # assign the haplotypes to the variants
    assign_haplotypes(instream, outstream, assignments)

    # close the input and output streams
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()


if __name__ == "__main__":

    # parse the commandline arguments
    parser = argparse.ArgumentParser(
        prog="haplotype_vcf",
        description=""".""")
    parser.add_argument(
        "-v", "--vcf", dest="vcf",
        type=str, nargs="?", default="stdin",
        help="The input VCF file.")
    parser.add_argument(
        "-o", "--output", dest="output",
        type=str, nargs="?", default="stdout",
        help="The output VCF file.")
    parser.add_argument(
        "-r", "--relations", dest="relations",
        type=str, nargs="?",
        help="""A file with the trio relations in the
        order child, father, mother.""")
    parser.add_argument(
        "--offspring", dest="offspring",
        type=str, nargs="*",
        help="The sample names of the offspring.")
    parser.add_argument(
        "--fathers", dest="fathers",
        type=str, nargs="*",
        help="The sample names of the fathers.")
    parser.add_argument(
        "--mothers", dest="mothers",
        type=str, nargs="*",
        help="The sample names of the mothers.")
    parser.set_defaults(func=haplotype_vcf)

    #
    args = parser.parse_args()
    args.func(args)
