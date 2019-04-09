
def consensus_header(header, rgmapping):
    """Create a header for the consensus writer."""
    # get the read group mapping
    rgids = list(rgmapping.keys())
    sids = list(set(list(rgmapping.values())))

    # make sure there is only one sample in the BAM file
    if len(sids) != 1:
        raise ValueError("Wrong number of samples in the BAM file")

    # check where to add the header
    rgidx = len(header)
    for idx, line in enumerate(header):
        if line.startswith("@RG"):
            rgidx = idx
            break

    # get a unique consensus code
    crgid = sids[0] + "_CNS"
    while crgid in rgids:
        crgid += "_1"

    # add an RG line to the header
    smline = "@RG\tID:{rgid}\tCN:ECB\tLB:{sample}\tSM:{sample}\tPL:ILLUMINA\tDS:NULL".format(
        rgid=crgid,
        sample=sids[0])
    header.insert(rgidx, smline)

    # return the new header and the consensus id
    return header, crgid
