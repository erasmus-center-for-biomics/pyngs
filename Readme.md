# PyNGSVCF

## Description

PyNGS is a package to work with NGS data formats. It contains multiple
submodules for the various formats, such as sam, vcf, bed and bases.

## Sripts

This library provides several scripts to work with SAM alignments,
CIGAR operations and VCF files.

### SAM scripts

#### cigar_to_bed.py

`cigar_to_bed.py` takes the alignments in a SAM file and interprets their
CIGAR operations. The desired CIGAR operations are written to a BED file.
The BED entries are decorated with at least the READNAME and CIGAR operation, but
other tags from the alignments can be added as well. `__sample__` is used to
add the sample identifier to the BED file.

### BED scripts

#### merge_bed_entries.py

`merge_bed_entries.py` merges BED entries from the same DNA fragment to a
single larger BED entry. Mostly this will be performed based on the read
name. Entries over multiple chromosomes will not be merged.

The BED entries should be sorted on the factor on which they should be merged.

## Installation instructions

To install this python module from source do

```bash
# install with all scripts
python setup.py install

# install only the libraries
python setup.py install_lib
```
