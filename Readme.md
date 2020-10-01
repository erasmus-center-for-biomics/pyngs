# pyngs

## Description

*pyngs* is a Python package to work with NGS data. The package contains parsers for commonly used
formats, such as FastA, FastQ, SAM and VCF. Furthermore several algorithms and scripts are
implemented in this package. The package is meant to remove some of the boilerplate code in
writing custom scripts that work with NGS data and to be widely useable. A single executable
is added to this packages (*pyngs_tools*) which contains several handy
tools for thematically diverse applications (see below).

## pyngs_tools

*pyngs_tools has a collection of tools which can be more widely applied in NGS workflows. Furthermore,
these tools provide an example or proof-of-concept for one or more of *pyngs* functions.

- *pyngs_tools cigar-to-bed* converts one or more CIGAR operations from [SAM alignments](https://github.com/samtools/hts-specs)
to a BED file. In this conversion tags from the SAM file can be added to the comment section of the BED file. The comment
section of the resulting BED is encoded in a GFF like style. The BED file can be used by [BEDtools](https://bedtools.readthedocs.io/en/latest/) to
generate a genome coverage for a sample or specific subset of the data. Furthermore these files can be used in
bedtools intersect to count the reads per annotation.

- *pyngs_tools merge-bed-entries* takes BED entries, such as those from *pyngs_tools cigar-to-bed*, and merges them
into the spanning region of the entire entry. This is especially handy in ChIP and ATAC-seq applications where one
would like to proceed with smaller DNA fragments.

- *pyngs_tools extract-from-columns* extracts fields from the columns of a tab-delimited file. This tool is able to
parse comment fields in the GTF and GFF format and extract particular tags from these. This tool is used in conjunction
with *bedtools intersect* and *pyngs_tools count-last-column* to quantify reads over genomic areas.

- *pyngs_tools count-last-column* counts the last column in a sorted tab-delimited file. This tool is useful for
quantifying reads over regions.

- *pngs_tools umi-consensus* creates consensus alignments from alignments annotated with the same (unique) molecular
identifier. This tool is especially useful in variant analyses where alignments are deduplicated prior to variant
calling. Unlike other tools, pngs_tools umi-consensus allows alignments to be merged that are not aligned to the
same locations. This feature allows for alignments to be merged that align on (slightly) different locations based on
mismatches in their first bases.

- *pngs_tools filter-by-samtag* filters alignments based on the value (or presence/absence) of a tag. With this tool,
alignments can be discarded prior to a processing step.

- *pyngs_tools haplotype-vcf* assigns bases to an allele in VCF files in which mother, father and child are present.
Multiple trios can be assigned in a single file and the alleles of the parents are separated per child. The VCF files
should be genotyped and thus need to have a GT tag. The maternal allele in the child is indicated in the *MAT* column
in the format. The paternal allele is in the *PAT* column. Both parents will be annotated with the *OFF* (from offspring) column in
which the transmitted base is indicated together with the sample-id in the form of *sample_id*/*base*.

### Examples

A quick quantification workflow with *pyngs_tools cigar-to-bed*, *pyngs_tools extract-from-columns*, and
*pyngs_tools count-last-column*. This workflow makes use of *samtools*, *bgzip*, *bedtools* and common
Linux command-line tools, such as *grep* and *sort*. This is an example and options like strandedness and
alignment filtering can easily be added at the appropriate steps (bedtools intersect and samtools view).

```bash
samtools view -h -F 0x0100 -F 0x0800 -F 0x0004 $input_bam_file | \
    pyngs_tools cigar-to-bed --cigar-operations M | \
    sort -k 1,1 -k 2,2n | \
    bgzip -c > $intermediate_bed_file

bedtools intersect -wao -sorted -a $intermediate_bed_file -b $gtf_reference_file | \
    pyngs_tools extract-from-columns --columns 8 14 3 --encodings none gtf gff --fields X gene_biotype,gene_name,gene_id,transcript_id READNAME --output $intermediate_tsv_file

cat $intermediate_tsv_file | \
    grep '^exon' | \
    pyngs_tools extract-from-columns --columns 0 1 2 4 --encodings none none none none --fields X X X X \
    sort | \
    pyngs_tools count-last-column --output $output_tsv_file
```

## Installation instructions

To install this python module from source do

```bash
# install with all scripts
python setup.py install

# install only the libraries
python setup.py install_lib
```
