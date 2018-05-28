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
can merge alignments with different start positions. Alignments for molecules
that are represented by only a single read are printed as is.

