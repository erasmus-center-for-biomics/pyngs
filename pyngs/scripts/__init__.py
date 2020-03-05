
# FastQ methods
from .fastq.trim import trim_fastq
from .fastq.extract_subreads import extract_subreads
from .fastq.add_umis import add_umis_to_readname

# SAM to BED methods
from .cigar_to_bed import cigar_to_bed
from .merge_bed_entries import merge_bed_entries

# SAM methods
from .consensus import run_consensus
from .filter_by_samtag import filter_by_samtag
from .sam.retrieve_tag_from_readname import retrieve_tag_from_readname

# VCF methods
from .vcf.haplotype_vcf import haplotype_vcf
from .vcf.b_allel_frequency import add_b_allel_frequency
from .vcf.log_r_ratio import add_log_r_ratio
from .vcf.tabulate_vcf import tabulate_vcf
from .vcf.vcf_to_table import vcf_to_table
from .vcf.mpileup_frequency_caller import mpileup_frequency_caller

# Tabular methods
from .extract_from_columns import extract_from_columns
from .count_last_column import count_last_column

