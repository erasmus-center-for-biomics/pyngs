import sys
import gzip
import argparse
import logging
from typing import Dict, List, Union, TextIO, Generator
import pyngs.vcf as vcf
from pyngs import open_stream
from pyngs.vcf.tools import FilterAlt


class FilterAlleles:

    def __init__(self, samples: List[str], alt:int = 0, freq: float=0.0, digits: int=5, ploidy: int=2) -> None:
        self.samples = [] if samples is None else samples
        self.alternate = alt
        self.frequency = freq
        self.afdigits = digits
        self.alternate_filter = FilterAlt(ploidy=ploidy)

    def filter(self, instream: TextIO, outstream: TextIO) -> None:
        """Filter the alleles."""
        reader = vcf.Reader(instream)
        metas = reader.meta

        # check if AD is present
        if "AD" not in reader.format.keys():
            raise ValueError("FORMAT/AD not found in the VCF file")

        # add the DP meta if not already present
        if "DP" not in reader.format.keys():
            dpstr = """##FORMAT=<ID=DP,Number=1,Type=Integer,Description="The total depth at this position">"""
            dpmeta = vcf.Meta.from_str(dpstr)
            dpformat = vcf.Format.from_meta(dpmeta)
            metas.append(dpmeta)

        # add the AF tag if not present
        if "AF" not in reader.format.keys():
            afstr = """##FORMAT=<ID=AF,Type=Float,Number=R,Description="The allele frequencies of the alleles">"""
            afmeta = vcf.Meta.from_str(afstr)
            afformat = vcf.Format.from_meta(afmeta)
            metas.append(afmeta)
        writer = vcf.Writer(outstream, metas, reader.samples)

        # determine which samples to check
        if self.samples:
            sindices = [i for i,v in enumerate(reader.samples) if v in self.samples]
        else:
            sindices = list(range(len(reader.samples)))

        # for each variant
        for variant in reader:
            to_keep = []

            # check whether the AD tag is present
            if not variant.fstore.has_tag("AD"):
                logging.warning(
                    "Variant %s has no FORMAT/AD tag, skipping",
                    variant.to_simple_repr())
                continue

            # check which formats to add
            formats_to_add = []
            if not variant.fstore.has_tag("DP"):
                formats_to_add.append(dpformat)
            if not variant.fstore.has_tag("AF"):
                formats_to_add.append(afformat)

            # add the formats
            for fmt in formats_to_add:
                variant.fstore.add_format(fmt)

            # for each of the samples we need to check
            for sidx, store in enumerate(variant.fstore.stores):

                # get the allelic depths
                adval = store["AD"]
                if adval is None:
                    store["DP"] = None
                    store["AF"] = None
                    continue

                # add the read depth and frequency
                dp = sum(adval)
                store["DP"] = [dp]
                if dp > 0:
                    store["AF"] = [round(v/dp, self.afdigits) for v in adval]
                else:
                    store["AF"] = None

                # skip non target samples
                if sidx not in sindices:
                    continue

                # check the allelic depths
                if self.alternate > 0:
                    altkeep = [i for i,v in enumerate(adval[1:]) if v >= self.alternate]
                    to_keep.extend(altkeep)

                # check the frequencies
                if store["AF"] is not None and self.frequency > 0.0:
                    frqkeep = [i for i,f in enumerate(store["AF"][1:]) if f >= self.frequency]
                    to_keep.extend(frqkeep)

            # summarise the alleles to keep
            to_keep = list(set(to_keep))

            # if we keep all alleles we do not need to filter. Otherwise,
            # if we do not keep any alleles, we do not have to write the variant
            if len(to_keep) == len(variant.alternates):
                writer.write(variant)
            elif not to_keep:
                continue

            # determine the alts to remove
            remove = [v for i,v in enumerate(variant.alternates) if i not in to_keep]
            try:
                nvariant = self.alternate_filter(variant, remove=remove)
            except TypeError:
                logging.warning(
                    "An error occured while processing variant {0}",
                    variant.to_simple_repr())
                continue
            writer.write(nvariant)


def allele_filter(args):
    """Run the script."""
    fltalleles = FilterAlleles(
        samples=args.samples,
        alt=args.alternate,
        freq=args.frequency,
        digits=args.digits,
        ploidy=args.ploidy)

    with open_stream(args.input, "rt") as instream:
        with open_stream(args.output, "wt") as outstream:
            fltalleles.filter(instream, outstream)


if __name__ == "__main__":
    sparser = argparse.ArgumentParser(
        prog="allele_filter",
        description="""Keep only alleles which fit the contraints.""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, default="stdin",
        help="""The input VCF file.""")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, default="stdout",
        help="The output VCF file.")
    sparser.add_argument(
        "-p", "--ploidy", dest="ploidy",
        type=int, default=2,
        help="The ploidy of the organism")
    sparser.add_argument(
        "-s", "--samples", dest="samples",
        type=str,
        help="The samples on which to base the filtering.")
    sparser.add_argument(
        "-a", "--alternate", dest="alternate",
        type=int, nargs="?", default=0,
        help="The minimum number of alternate alleles.")
    sparser.add_argument(
        "-f", "--frequency", dest="frequency",
        type=float, nargs="?", default=0.0,
        help="The minimum frequency of alternate alleles.")
    sparser.add_argument(
        "-d", "--digits", dest="digits",
        type=int, nargs="?", default=5,
        help="The number of digits in the AF field.")
    sparser.set_defaults(func=allele_filter)
    args = sparser.parse_args()
    args.func(args)
