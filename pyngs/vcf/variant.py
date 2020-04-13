from collections import namedtuple, OrderedDict
from typing import List, Union, Tuple, Optional, Dict, TypeVar

from .utils import genotypes
from .values import ParseVcfTypes, ReprVcfTypes, VcfValue,  VCFTYPES, VCFNUMBERS
from .format import Format
from .info import Info
from .stores import FormatStore, InfoStore


V = TypeVar("V")

class Variant:

    def __init__(self,
            chrom: str, pos: int, ids: str, ref: str, alt: List[Optional[str]],
            qual: Optional[float], filters: Optional[List[str]],
            infostore: InfoStore, formatstore: FormatStore) -> None:
        self.chrom = chrom
        self.position = pos
        self.id = ids
        self.reference = ref
        self.alternates = alt
        self.quality = qual
        self.filter = filters
        self.istore = infostore
        self.fstore = formatstore

    @classmethod
    def from_str(cls: V, strval: str, info_dct: Dict[str, Info], format_dct: Dict[str,Format]) -> V:
        """Parse a Variant from a string."""
        fields = strval.rstrip().split("\t")

        # prepare the alternates
        alternates = fields[4].split(",")
        quality = float(fields[5]) if fields[5] != "." else None
        flts = None if fields[6] == "." else fields[6].split(",")

        # prepare the info store
        infostore = InfoStore()
        if fields[7] != ".":
            infostore.add_data(fields[7], info_dct)

        # prepare the format store
        formatstore = FormatStore()
        if len(fields) > 8:
            formatstore.add_formats_from_dict(fields[8], format_dct)

            formatstore.init_stores(len(fields) - 9)
            for sidx, value in enumerate(fields[9:]):
                try:
                    formatstore.add_sample_data(sidx, value)
                except ValueError:
                    formatstore.add_sample_data(
                        sidx,
                        ":".join(["."] * len(formatstore.format)))

        # return the variant
        return cls(
            chrom=fields[0],
            pos=int(fields[1]),
            ids=fields[2],
            ref=fields[3],
            alt=alternates,
            qual=quality,
            filters=flts,
            infostore=infostore,
            formatstore=formatstore)

    def __str__(self) -> str:
        """Create a str from the variant."""
        fields = [
            self.chrom,
            str(self.position),
            self.id,
            self.reference,
            ",".join(self.alternates),
            str(self.quality) if self.quality is not None else ".",
            ",".join(self.filter) if self.filter is not None else ".",
            str(self.istore)]
        if self.fstore.format:
            fields.append(str(self.fstore))
        return "\t".join(fields) + "\n"

    def genotypes(self, ploidy: int=2) -> List[List[str]]:
        """Get the genotypes for the current variant."""
        genoint = [g for g in genotypes(ploidy, len(self.alternate))]
        alleles = [self.reference] + self.alternate

        retval = []
        for geno in genoint:
            retval.apped([alleles[i] for i in geno])
        return retval

    # @classmethod
    # def from_row(cls, row: Row, index: HeaderIndex):

    #     # parse the alternate alleles
    #     alternates = [a for a in quote_tokenizer(row.alternate, ",")]

    #     # parse the info fields
    #     info: List[Tuple[str, str, VcfValue]] = []
    #     for key, val in row.info:
    #         value = None
    #         pa = index.get("INFO", key)
    #         if pa.type == "Flag":
    #             info.append((key, pa.type, pa.number, []))
    #             continue

    #         value = []
    #         for v in quote_tokenizer(val, ","):
    #             value.append(ParseVcfTypes[pa.type](v))
    #         info.append((key, pa.type, pa.number, value))

    #     # parse the sample and format fields
    #     formats = None
    #     samples = None
    #     if row.format is not None:
    #         formats = []
    #         samples = [[None] * len(row.format)] * len(row.samples)
    #         for fidx, key in enumerate(row.format):
    #             pa  = index.get("FORMAT", key)
    #             func = ParseVcfTypes[pa.type]
    #             formats.append((key, pa.type, pa.number))

    #             # parse the values per sample
    #             for sidx, sample in enumerate(row.samples):
    #                 value = []
    #                 for v in quote_tokenizer(sample[fidx], ","):
    #                     value.append(func(v))
    #                 samples[sidx][fidx] = value
    #     # return a variant
    #     return Variant(
    #             row.chrom,
    #             row.position,
    #             row.id,
    #             row.reference,
    #             alternates,
    #             row.quality,
    #             row.filter,
    #             info,
    #             formats,
    #             samples)

    # def filter_alternate(self, altidx: int, ploidy: int=2):
    #     """Return a variant with alternate index altidx removed."""

    #     gt_idx_to_rm = []
    #     for idx, geno in enumerate(genotypes(ploidy, len(self.alternate))):
    #         if altidx in geno:
    #             gt_idx_to_rm.append(idx)

    #     alternates = [a for i,a in enumerate(self.alternate) if i != altidx]

    #     # filter the info
    #     info = []
    #     for entry in self.info:
    #         if entry[1] == "Flag":
    #             info.append(entry)
    #             continue

    #         values = []

    #         if entry[2] == "G":
    #             for idx, val in enumerate(entry[3]):
    #                 if idx not in gt_idx_to_rm:
    #                     values.append(val)

    #         elif entry[2] == "R":
    #             for idx, val in enumerate(entry[3]):
    #                 if idx != altidx + 1:
    #                     values.append(val)

    #         elif entry[2] == "A":
    #             for idx, val in enumerate(entry[3]):
    #                 if idx != altidx:
    #                     values.append(val)

    #         else:
    #             values = entry[3]
    #         info.append((entry[0], entry[1], entry[2], values))

    #     # filter the samples
    #     samples = [[None] * len(self.format)] * len(self.samples)
    #     for fidx, ftype in enumerate(self.format):
    #         for sidx, svalue in enumerate(self.samples):
    #             values = []
    #             if ftype[2] == "G":

    #                 for idx, val in enumerate(svalue[fidx]):
    #                     if idx not in gt_idx_to_rm:
    #                         values.append(val)

    #             elif ftype[2] == "R":
    #                 for idx, val in enumerate(svalue[fidx]):
    #                     if idx  != altidx + 1:
    #                         values.append(val)
    #                 samples[sidx][fidx] = values

    #             elif ftype[2] == "A":
    #                 for idx, val in enumerate(svalue[fidx]):
    #                     if idx  != altidx:
    #                         values.append(val)
    #                 samples[sidx][fidx] = values
    #             else:
    #                 values = svalue[fidx]

    #             samples[sidx][fidx] = values

    #     # return the modified variant
    #     return Variant(
    #         self.chrom,
    #         self.position,
    #         self.id,
    #         self.reference,
    #         alternates,
    #         self.quality,
    #         self.filter,
    #         info,
    #         self.format,
    #         samples)

    # def __repr__(self):
    #     """Return a string representation of this variant."""

    #     info = []
    #     for entry in self.info:
    #         if entry[1] == "Flag":
    #             info.append(entry[0])
    #             continue

    #         func = ReprVcfTypes[entry[1]]
    #         entry = "{0}={1}".format(
    #             entry[0],
    #             ",".join([func(v) for v in entry[3]]))
    #         info.append(entry)

    #     alt = "."
    #     if self.alternate:
    #         alt = ",".join(self.alternate)
    #     #
    #     bstr = "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}"
    #     vstr = bstr.format(
    #         chrom=self.chrom,
    #         pos=self.position,
    #         id=self.id if self.id is not None else ".",
    #         ref=self.reference,
    #         alt=alt,
    #         qual=self.quality if self.quality is not None else ".",
    #         filter=";".join(self.filter),
    #         info=";".join(info))

    #     if self.format is not None:
    #         funcs = []
    #         for fentry in self.format:
    #             funcs.append(ReprVcfTypes[fentry[1]])

    #         formatstr = ":".join([f[0] for f in self.format])
    #         vstr += "\t{0}".format(formatstr)

    #         for sample in self.samples:
    #             data = []
    #             for idx, fentry in enumerate(self.format):
    #                 val = ",".join([funcs[idx](v) for v in sample[idx]])
    #                 data.append(val)
    #             vstr += "\t{0}".format(":".join(data))
    #     #
    #     return vstr

    # def guess_ploidy(self):
    #     """Guess the ploidy based on a the length of a genotype info field."""
    #     total
