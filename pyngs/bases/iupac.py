"""Methods to deal with IUPAC encode bases."""

IUPAC_BASES = (
    "A", "C",
    "T", "G",
    "U",
    "R", "Y",
    "S", "W",
    "K", "M",
    "B", "D",
    "H", "V",
    "N")


def complement(base: str) -> (str):
    """
    Complement the base with respect to the IUPAC code.
    :param base: the base to complement
    :return: the complemented base
    """
    # assert base in IUPAC_BASES
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'G':
        return 'C'
    elif base == 'C':
        return 'G'
    elif base == 'U':
        return 'A'
    elif base == 'N':
        return 'N'
    elif base == 'R':
        return 'Y'
    elif base == 'Y':
        return 'R'
    elif base == 'W':
        return 'S'
    elif base == 'S':
        return 'W'
    elif base == 'K':
        return 'M'
    elif base == 'M':
        return 'K'
    elif base == 'B':  # CGT
        return 'V'  # GCT
    elif base == 'D':  # AGT
        return 'H'  # TCA
    elif base == 'H':  # ACT
        return 'D'  # TGA
    elif base == 'V':  # ACG
        return 'B'  # TGC
    else:
        raise ValueError(
            "base {0} is not present in the IUPAC alphabet".format(
                base))


def reverse_complement(seq: str) -> (str):
    """
    Reverse and complement a sequence.
    :param seq: the sequence to reverse complement
    :return: the reverse complemented sequence
    """
    retval = [complement(b) for b in seq]
    retval.reverse()
    return ''.join(retval)


def to_regexp(seq: str) -> (str):
    """
    Convert IUPAC to regular expresions.

    Decodes a sequence which is IUPAC and convert
    this to a regular expression friendly sequence.

    :param seq: the sequence to encode
    :return: the regular expression
    """
    # convert IUPAC bases
    seq = seq.replace('R', '[A,G]')
    seq = seq.replace('Y', '[C,T]')
    seq = seq.replace('S', '[G,C]')
    seq = seq.replace('W', '[A,T]')
    seq = seq.replace('K', '[G,T]')
    seq = seq.replace('M', '[A,C]')
    seq = seq.replace('B', '[C,G,T]')
    seq = seq.replace('D', '[A,G,T]')
    seq = seq.replace('H', '[A,C,T]')
    seq = seq.replace('V', '[A,C,G]')
    seq = seq.replace('N', '.')
    # return the sequence
    return seq


def expand_iupac(base: str, fill_n: bool=False) -> (str):
    """
    Expand the IUPAC base
    :param base: the IUPAC base to expand
    :param fill_n: should we fill N or leave it empty
    :return: a string with all the primary
             bases that are encoded by the
             IUPAC bases
    """
    if base == "N":
        if not fill_n:
            return ""
        else:
            return "ACGT"
    elif base == "V":
        return "ACG"
    elif base == "H":
        return "ACT"
    elif base == "D":
        return "AGT"
    elif base == "B":
        return "CGT"
    elif base == "M":
        return "AC"
    elif base == "K":
        return "GT"
    elif base == "S":
        return "GC"
    elif base == "W":
        return "AT"
    elif base == "R":
        return "AG"
    elif base == "Y":
        return "CT"
    else:
        return base


def to_list(seq: str, fill_n: bool= False) -> (list):
    """Convert a sequence to a list."""
    retval = []
    for base in seq:
        retval.append(expand_iupac(base, fill_n))
    return retval
