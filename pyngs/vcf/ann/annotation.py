

from typing import List
from ..utils import quote_tokenizer


class Annotation:
    """
    Represent an annotation field as specified at
    http://snpeff.sourceforge.net/SnpEff_manual.html#ann.
    """
    __slots__ = [
        "allele",
        "annotation",
        "impact",
        "gene_name",
        "gene_id",
        "feature_type",
        "feature_id",
        "transcript_biotype",
        "rank",
        "total",
        "hgvs_c",
        "hgvs_p",
        "cdna_pos",
        "cdna_len",
        "cds_pos",
        "cds_len",
        "prot_pos",
        "prot_len",
        "dist_feature",
        "messages"
    ]

    def __init__(self,
            allele: str, annotation: List[str], impact: str,
            gene_name: str, gene_id: str, feature_type: str,
            feature_id: str, transcript_biotype: str,
            rank: int, total: int, hgvs_c: str,
            hgvs_p: str,
            cdna_pos: int, cdna_len: int,
            cds_pos: int, cds_len: int,
            prot_pos: int, prot_len: int,
            dist_feature: int, messages: str) -> None:
        self.allele = allele
        self.annotation = annotation
        self.impact = impact
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.feature_type = feature_type
        self.feature_id = feature_id
        self.transcript_biotype = transcript_biotype
        self.rank = rank
        self.total = total
        self.hgvs_c = hgvs_c
        self.hgvs_p = hgvs_p
        self.cdna_pos = cdna_pos
        self.cdna_len = cdna_len
        self.cds_pos = cds_pos
        self.cds_len = cds_len
        self.prot_pos = prot_pos
        self.prot_len = prot_len
        self.dist_feature = dist_feature
        self.messages = messages

    @classmethod
    def from_str(cls, text: str) -> "Annotation":
        """Get an annotation from a string."""
        flds = [t for t in quote_tokenizer(text, sep="|")]
        if len(flds) != 16:
            raise ValueError("{0} instead of 16 fields found in annotation".format(len(flds)))

        allele = flds[0]
        anno = flds[1].split("&")
        impact = flds[2]
        gene_name = flds[3]
        gene_id = flds[4]
        feature_type = flds[5]
        feature_id = flds[6]
        transcript_biotype = flds[7]
        rank_total = flds[8].split("/")
        rank = int(rank_total[0]) if flds[8] != "" else -1
        total = int(rank_total[1]) if flds[8] != "" else -1
        hgvs_c = flds[9]
        hgvs_p = flds[10]
        cdna = flds[11].split("/")
        cdna_pos = int(cdna[0]) if flds[11] != "" else -1
        cdna_len = int(cdna[1]) if flds[11] != "" else -1
        cds = flds[12].split("/")
        cds_pos = int(cds[0]) if flds[12] != "" else -1
        cds_len = int(cds[1]) if flds[12] != "" else -1
        prot = flds[13].split("/")
        prot_pos = int(prot[0]) if flds[13] != "" else -1
        prot_len = int(prot[1]) if flds[13] != "" else -1
        dist_feature = int(flds[14]) if flds[14] != "" else -1
        messages = flds[15]
        return cls(
            allele, anno, impact, gene_name, gene_id,
            feature_type, feature_id, transcript_biotype,
            rank, total, hgvs_c, hgvs_p,
            cdna_pos, cdna_len,
            cds_pos, cds_len,
            prot_pos, prot_len,
            dist_feature, messages)

    def __repr__(self) -> str:
        """Create a string representation of the Annotation."""
        rank_total = "" if self.rank == -1 else "{0}/{1}".format(self.rank, self.total)
        cdna = "" if self.cdna_pos == -1 else "{0}/{1}".format(self.cdna_pos, self.cdna_len)
        cds = "" if self.cds_pos == -1 else "{0}/{1}".format(self.cds_pos, self.cds_len)
        prot = "" if self.prot_pos == -1 else "{0}/{1}".format(self.prot_pos, self.prot_len)
        return "{allele}|{annotation}|{impact}|{gene_name}|{gene_id}|{feature_type}|{feature_id}|{transcript_biotype}|{rank_total}|{hgvs_c}|{hgvs_p}|{cdna}|{cds}|{prot}|{dist_feature}|{messages}".format(
            allele = self.allele,
            annotation  = "&".join(self.annotation),
            impact = self.impact,
            gene_name = self.gene_name,
            gene_id = self.gene_id,
            feature_type = self.feature_type,
            feature_id = self.feature_id,
            transcript_biotype = self.transcript_biotype,
            rank_total = rank_total,
            hgvs_c = self.hgvs_c,
            hgvs_p = self.hgvs_p,
            cdna = cdna,
            cds = cds,
            prot = prot,
            dist_feature = self.dist_feature if self.dist_feature != -1 else "",
            messages = self.messages)

