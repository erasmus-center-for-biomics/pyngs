import sys
import gzip
import math
import pyngs.vcf


class VCFToTable:

    def __init__(self, tags, sep="\t"):
        self.fields = []
        fields = []
        for idx, tag in enumerate(tags):
            if "/" in tag:
                [section, name] = tag.split('/')
                if section in ('INFO', 'FORMAT'):
                    fields.append((tag, section, name))
                    continue
            else:
                fields.append((tag, tag, '.'))
        
        self.sep = sep

    def __prepare_header__(self):
        """Prepare the header."""
        extra = [f[0] for f in self.fields]
        return self.sep.join([
            "chrom",
            "position",
            "reference",
            "alternate",
            "sample",
            self.sep.join(extra)])

    def __set_fields__(self, variant, values):
        """Set the field dependent values based on a variant."""        
        infod = dict(variant.info)
        for fidx, field in enumerate(self.fields):
            if field[1] == "INFO":
                try:
                    values[fidx] = str(infod[field[2]])
                except KeyError:
                    values[fidx] = "."
            elif field[1] == 'ID':
                values[fidx] = variant.id
            elif field[1] == 'QUAL':
                values[fidx] = str(variant.quality)
            elif field[1] == 'Filter':
                values[fidx] = str(variant.filter)


    def __call__(self, instream, outstream):
        """Convert a VCF file to a table."""
        reader = pyngs.vcf.Reader(instream)
        outstream.write("{0}\n".format(self.__prepare_header__()))
        
        for variant in reader:

            # prepare the default output
            values = ["."] * len(self.fields)

            # parse the other fields
            self.__set_fields__(variant, values)

            # parse the format field per sample
            for sidx, sample in enumerate(variant.samples):
                for fidx, field in enumerate(self.fields):
                    if field[1] == 'FORMAT':
                        try:
                            values[fidx] = sample[variant.format.index(field[2])]
                        except ValueError:
                            values[fidx] = "."

                # write the output per sample here
                line = self.sep.join([
                    variant.chrom,
                    str(variant.position),
                    variant.reference,
                    variant.alternate,
                    sample,
                    self.sep.join(values)])
                outstream.write("{0}\n".format(line))


def vcf_to_table(args):
    """Convert a VCF file to a tab delimited file."""
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
    
    # parse the VCF file
    vcftotable = VCFToTable(tags=args.tags, sep=args.sep)
    vcftotable(instream, outstream)

    # close the in and output files
    if instream is not sys.stdin:
        instream.close()
    if outstream is not sys.stdout:
        outstream.close()
        