from collections import namedtuple, OrderedDict
from typing import List, Union, Tuple, Optional, Dict, TypeVar
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

    def to_simple_repr(self) -> str:
        """Convert the variant to a simple string."""
        return "{0}:{1}:{2}-{3}".format(
            self.chrom,
            str(self.position),
            self.reference,
            ",".join(self.alternates))

    def samples(self) -> int:
        """Return the number of samples."""
        return len(self.fstore.stores)

    def genotypes(self, ploidy: int=2) -> List[List[str]]:
        """Get the genotypes for the current variant."""
        genoint = [g for g in genotypes(ploidy, len(self.alternate))]
        alleles = [self.reference] + self.alternate

        retval = []
        for geno in genoint:
            retval.apped([alleles[i] for i in geno])
        return retval

