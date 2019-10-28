
# FastQ methods
from .fastq.trim import trim_fastq
from .fastq.extract_subreads import extract_subreads

# SAM to BED methods
from .cigar_to_bed import cigar_to_bed
from .merge_bed_entries import merge_bed_entries

# SAM methods
from .consensus import run_consensus
from .filter_by_samtag import filter_by_samtag
from .sam.retrieve_tag_from_readname import retrieve_tag_from_readname

# VCF methods
from .haplotype_vcf import haplotype_vcf
from .b_allel_frequency import add_b_allel_frequency
from .log_r_ratio import add_log_r_ratio
from .tabulate_vcf import tabulate_vcf

# Tabular methods
from .extract_from_columns import extract_from_columns
from .count_last_column import count_last_column

