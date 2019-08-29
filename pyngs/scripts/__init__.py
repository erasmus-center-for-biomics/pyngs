
# FastQ methods
from .trim_fastq import trim_fastq

# SAM to BED methods
from .cigar_to_bed import cigar_to_bed
from .merge_bed_entries import merge_bed_entries

# SAM methods
from .consensus import run_consensus
from .filter_by_samtag import filter_by_samtag

# VCF methods
from .haplotype_vcf import haplotype_vcf
from .b_allel_frequency import add_b_allel_frequency
from .log_r_ratio import add_log_r_ratio
from .tabulate_vcf import tabulate_vcf

# Tabular methods
from .extract_from_columns import extract_from_columns
from .count_last_column import count_last_column

