
# FastQ methods
from .fastq.trim import trim_fastq
from .fastq.extract_subreads import extract_subreads
from .fastq.add_umis import add_umis_to_readname
from .fastq.iupac_subreads import iupac_subreads

# SAM to BED methods
from .cigar_to_bed import cigar_to_bed
from .merge_bed_entries import merge_bed_entries

# SAM methods
from .consensus import run_consensus
from .filter_by_samtag import filter_by_samtag
from .sam.retrieve_tag_from_readname import retrieve_tag_from_readname

# VCF methods
# from .vcf.haplotype_vcf import haplotype_vcf
# from .vcf.b_allel_frequency import add_b_allel_frequency
# from .vcf.log_r_ratio import add_log_r_ratio
# from .vcf.mpileup_frequency_caller import mpileup_frequency_caller
# from .vcf.clean_caller import clean_caller
from .vcf.clean_alt import clean_alt
from .vcf.add_frequency import add_frequency
from .vcf.allele_filter import allele_filter
from .vcf.divide_ann import parse_annotations

# Tabular methods
from .extract_from_columns import extract_from_columns
from .count_last_column import count_last_column

