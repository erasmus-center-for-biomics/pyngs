pyngs
=====

Tools
-----

* `scripts/add_umi_to_readname.py`
* `scripts/add_umi_as_samtag.py`
* `scripts/umi_consensus.py`

### UMI consensus

The UMI consensus tool deduplicates paired-end alignments that are coupled
to Unique Molecular Identifiers (UMIs). These UMIs indicate individual DNA
molecules that were tagged before PCR amplifications. Disagreements in the
sequence of DNA fragments with the same UMI are therefore likely due to
errors introduced in the sample preparation or the sequencing.

This tool generates consensus alignments from reads that align to the same
locus and have the same UMI. The consensus alignments are fault tolerant and
can merge alignments with different start positions.

The alignments to be merged are required to have the same flag as these properties
cannot reliably be merged without performing a new alignment. Supplementary and
secondary alignments are also merged, but from these the mate chromosome, position
and tlen fields are not added.

Alignments for molecules that are represented by only a single read or read-pair
are printed as is.

The `umi_consensus` script can only be applied to SAM files or streams
that are sorted on the UMI tag first and then on the alignment position. To sort
alignments on the UMI tag (default: `um`), use the following command:

```bash
samtools sort -t um -o output.bam input.bam
```

To subsequently calculate the consensus from the UMIs using an inline BAM to
SAM conversion, use the following command:

```bash
samtools view -h output.bam | python3 umi_consensus.py -t um --output consensus.sam
```
