import sys
import gzip
import argparse
import logging

from itertools import groupby
from typing import Dict, List, Union, TextIO, Generator
from operator import attrgetter, itemgetter

import pyngs.vcf as vcf
from pyngs import open_stream


"""A map between the ANN fields and the VCF header."""
FIELD_MAP = {
    "allele": ("ann_allele", ""),
    "annotation": ("ann_annotation", ""),
    "impact": ("ann_impact",""),
    "gene_name": ("ann_gene_name", ""),
    "gene_id": ("ann_gene_id", ""),
    "feature_type": ("ann_feature_type", ""),
    "feature_id": ("ann_feature_id", ""),
    "transcript_biotype": ("ann_transcript_biotype", ""),
    "rank": ("ann_rank", ""),
    "total": ("ann_total", ""),
    "hgvs_c": ("ann_hgvs_c", ""),
    "hgvs_p": ("ann_hgvs_p", ""),
    "cdna_pos": ("ann_cdna_pos", ""),
    "cdna_len": ("ann_cdna_len", ""),
    "cds_pos": ("ann_cds_pos", ""),
    "cds_len": ("ann_cds_len", ""),
    "prot_pos": ("ann_prot_pos", ""),
    "prot_len": ("ann_prot_len", ""),
    "dist_feature": ("ann_dist_feature", ""),
    "messages": ("ann_messages","")
}

AnnotationTypes = Union[List[str], str, int]

class AnnotationParser:

    def __init__(self, fields: List[str], impacts: List[str], unordered: bool= False, anntag: str="ANN"):
        """."""
        self.anntag = anntag
        self.number = "." if unordered else "A"
        self.impacts = [i for i in impacts]

        self.attrgetters = []
        self.fields = [f for f in fields]
        for fld in self.fields:
            self.attrgetters.append(attrgetter(fld))
            if fld not in FIELD_MAP.keys():
                raise ValueError("field {0} not found".format(fld))

    def convert(self, instream: TextIO, outstream: TextIO) -> None:
        """Extract the required fields from the annotation field."""
        reader = vcf.Reader(instream)

        # add the fields to extract from the ANN field
        metas = reader.meta
        infos = []
        for fld in self.fields:

            # check whether the tag is already present
            tag = FIELD_MAP[fld][0]
            if tag in reader.info.keys():
                raise ValueError("INFO/{0} already present in the header".format(tag))

            # create a new header for the meta data
            hdrstr = """##INFO=<ID={0},Number={1},Type=String,Description="{2}">""".format(
                tag, self.number, FIELD_MAP[fld][1])
            meta = vcf.Meta.from_str(hdrstr)
            metas.append(meta)
            infos.append(vcf.Info.from_meta(meta))

        # check if the annotation tag is present
        if self.anntag not in reader.info.keys():
            raise ValueError("INFO/{0} not found in the VCF file".format(self.anntag))

        # prepare the writer
        writer = vcf.Writer(outstream, metas, reader.samples)

        # foreach variant
        for variant in reader:

            # if the variant has no annotations
            if self.anntag not in variant.istore.info.keys():
                writer.write(variant)
                continue

            annlst = variant.istore.data[self.anntag]
            if annlst is None:
                writer.write(variant)
                continue

            # get the annotations per allele
            annotations: Dict[str, AnnotationTypes] = {}
            for ann in annlst:
                anno = vcf.Annotation.from_str(ann)

                # skip annotations with a non-specified impact
                if self.impacts and anno.impact not in self.impacts:
                    continue

                # add alleles if required
                if anno.allele not in annotations.keys():
                    annotations[anno.allele] = []
                annotations[anno.allele].append(anno)

            if not annotations:
                writer.write(variant)
                continue

            # deduplicate the values
            deduplicated = {}
            for allele, annos in annotations.items():
                temp = []
                for anno in annos:
                    entry = []
                    for getter in self.attrgetters:
                        val = getter(anno)
                        if isinstance(val, list):
                            val = "&".join(val)
                        entry.append(val)
                    temp.append(entry)
                temp.sort()
                deduplicated[allele] = [temp[0]]
                for idx, values in enumerate(temp):
                    if idx == 0:
                        continue
                    if values != deduplicated[allele][-1]:
                        deduplicated[allele].append(values)

            # add the alleles
            for idx, info in enumerate(infos):

                entries = []
                if self.number == "A":

                    for alt in variant.alternates:
                        if alt not in deduplicated.keys():
                            entries.append(".")
                        else:
                            val = [str(l[idx]) for l in deduplicated[alt]]
                            entries.append("|".join(val))
                else:
                    # ignore allele order and print in the order obtained
                    for dd in deduplicated.values():
                        val = [str(l[idx]) for l in dd]
                        entries.append("|".join(val))

                # add the info to the variant
                if entries:
                    variant.istore.add_info(info)
                    variant.istore.data[info.code] = entries

            # write the values
            writer.write(variant)


def parse_annotations(args):
    """."""
    ap = AnnotationParser(
        fields=args.fields,
        impacts=args.impacts)
    with open_stream(args.input, "rt") as instream:
        with open_stream(args.output, "wt") as outstream:
            ap.convert(instream, outstream)


if __name__ == "__main__":
    sparser = argparse.ArgumentParser(
        prog="parse_annotations",
        description="""
# Annotation parser: a parser for the ANN field from VCF files

This parser is based on the explanation at http://snpeff.sourceforge.net/SnpEff_manual.html#ann.

The parser can be used as a python library to parse individual ANN fields or as a stand-alone
tool.

The stand-alone tool supports VCF files. It putsthe selected fields in the VCF INFO
field. For the VCF parsing the tool requires access to the pyngs library. In the default mode,
the output fields are reported per alternative allele (A). If the snpEff tool was run in cancer mode,
set the --unordered, so the tool will set the Number parameter in header as "." (instead of R).

To select specific ANN fields, a syntax is used like bcftools annotations with the following fields.

| Field | CLI identifier     | field in VCF           |
|:------|:------------------:|:-----------------------|
|  1    | allele             | ann_allele             |
|  2    | annotation         | ann_annotation         |
|  3    | impact             | ann_impact             |
|  4    | gene_name          | ann_gene_name          |
|  5    | gene_id            | ann_gene_id            |
|  6    | feature_type       | ann_feature_type       |
|  7    | feature_id         | ann_feature_id         |
|  8    | transcript_biotype | ann_transcript_biotype |
|  9    | rank               | ann_rank               |
|  9    | total              | ann_total              |
| 10    | hgvs_c             | ann_hgvs_c             |
| 11    | hgvs_p             | ann_hgvs_p             |
| 12    | cdna_pos           | ann_cdna_pos           |
| 12    | cdna_len           | ann_cdna_len           |
| 13    | cds_pos            | ann_cds_pos            |
| 13    | cds_len            | ann_cds_len            |
| 14    | prot_pos           | ann_prot_pos           |
| 14    | prot_len           | ann_prot_len           |
| 15    | dist_feaure        | ann_dist_feature       |
| 16    | messages           | ann_messages           |

""")
    sparser.add_argument(
        "-i", "--input", dest="input",
        type=str, default="stdin",
        help="""The input VCF file.""")
    sparser.add_argument(
        "-o", "--output", dest="output",
        type=str, default="stdout",
        help="The output VCF file.")
    sparser.add_argument(
        "-f", "--fields", dest="fields",
        type=str, nargs="+", default=[],
        help="The fields from the ANN field to select")
    sparser.add_argument(
        "--unordered", dest="unordered",
        type=bool, default=False,
        help="Do not match the output to the alternate alleles")
    sparser.add_argument(
        "--impacts", dest="impacts",
        type=str, nargs="*", default=[],
        help="Only report annotations with the specified impacts")
    sparser.set_defaults(func=parse_annotations)
    args = sparser.parse_args()
    args.func(args)
