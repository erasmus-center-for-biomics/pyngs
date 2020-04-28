import logging
from typing import List, Tuple, Optional
import pyngs.vcf as vcf


class Haplotyper:

    def __init__(self, idx_offspring: int, idx_paternal: int, idx_maternal: int, tag: str="GT") -> None:
        """."""
        self.offspring = idx_offspring
        self.paternal = idx_paternal
        self.maternal = idx_maternal

    def haplotype(self, variant: vcf.Variant) -> Optional[Tuple[str,str]]:
        """Haplotype the variant."""
        if not variant.fstore.has_tag(self.tag):
            return None
        if variant.fstore.stores[self.offspring]["GT"] == None:
            return None
        if variant.fstore.stores[self.paternal]["GT"] == None:
            return None
        if variant.fstore.stores[self.maternal]["GT"] == None:
            return None

        #
        offspring = variant.fstore.stores[self.offspring]["GT"][0].split("/")
        paternal = variant.fstore.stores[self.paternal]["GT"][0].split("/")
        maternal = variant.fstore.stores[self.maternal]["GT"][0].split("/")

        if offspring == [".", "."]:
            return None
        if paternal == [".", "."]:
            return None
        if maternal == [".", "."]:
            return None

        if offspring[0] == offspring[1]:
            if offspring[0] in paternal and offspring[0] in maternal:
                return ["P", "M"]
            else:
                return None
        else:
            if offspring[0] in paternal and offspring[1] in maternal:
                return ["P", "M"]
            elif offspring[1] in paternal and offspring[0] in maternal:
                return ["M", "P"]
            else:
                return None

class FilterAlt:

    def __init__(self, ploidy: int=2):
        """Initialize the object."""
        self.ploidy = ploidy

    def __call__(self, variant: vcf.Variant, remove: List[str]) -> vcf.Variant:
        """Filter the other alt alleles."""

        # get the alternates to keep
        atokeep = [i for i,v in enumerate(variant.alternates) if v not in remove]
        if len(atokeep) == len(variant.alternates):
            return variant

        # get the indexes of the R types to keep
        rtokeep = [0] + [v + 1 for v in atokeep]

        # get the genotypes
        genotypes = vcf.genotypes(self.ploidy, alleles=len(variant.alternates))
        genotokeep = []
        for idx, geno in enumerate(genotypes):
            keep = True
            for allele in geno:
                if allele not in rtokeep:
                    keep = False
            if keep:
                genotokeep.append(idx)

        # copy and filter the info
        nistore = vcf.InfoStore()
        for info_p in variant.istore.info.values():
            nistore.add_info(info_p)
        for code, value in variant.istore.data.items():
            try:
                if nistore.info[code].number == "G":
                    nistore.data[code] = [value[i] for i in genotokeep]
                elif nistore.info[code].number == "R":
                    nistore.data[code] = [value[i] for i in rtokeep]
                elif nistore.info[code].number == "A":
                    nistore.data[code] = [value[i] for i in atokeep]
                elif value is None:
                    nistore.data[code] = value
                else:
                    nistore.data[code] = value
            except TypeError:
                logging.warning(
                    "Malformed field %s for variant %s with value %s",
                    code, variant.to_simple_repr(), str(value))

        # copy and filter the format
        nfstore = vcf.FormatStore()
        nfstore.init_stores(variant.samples())
        for fmt in variant.fstore.format:
            nfstore.add_format(fmt)
            for sidx, store in enumerate(variant.fstore.stores):
                value = store[fmt.code]
                if fmt.number == "G":
                    nfstore.stores[sidx][fmt.code] = [value[i] for i in genotokeep]
                elif fmt.number == "R":
                    nfstore.stores[sidx][fmt.code] = [value[i] for i in rtokeep]
                elif fmt.number == "A":
                    nfstore.stores[sidx][fmt.code] = [value[i] for i in atokeep]
                else:
                    nfstore.stores[sidx][fmt.code] = value

        alternates = [variant.alternates[i] for i in atokeep]
        if not alternates:
            alternates = ["."]

        return vcf.Variant(
            chrom=variant.chrom,
            pos=variant.position,
            ids=variant.id,
            ref=variant.reference,
            alt=alternates,
            qual=variant.quality,
            filters=variant.filter,
            infostore=nistore,
            formatstore=nfstore)
