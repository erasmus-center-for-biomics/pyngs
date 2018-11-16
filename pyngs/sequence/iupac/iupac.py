"""Methods to deal with IUPAC encode bases."""

IUPAC_BASES = (
    "A", "C",
    "T", "G",
    "R", "Y",
    "S", "W",
    "K", "M",
    "B", "D",
    "H", "V",
    "N")


def complement(base: str) -> (str):
    """Complement the base with respect to the IUPAC code."""
    assert base in IUPAC_BASES
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'G':
        return 'C'
    elif base == 'C':
        return 'G'
    elif base == 'N':
        return 'N'
    elif base == 'R':
        return 'Y'
    elif base == 'Y':
        return 'R'
    elif base == 'W':
        return 'W'
    elif base == 'S':
        return 'S'
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


def to_regexp(seq: str) -> (str):
    """
    Convert IUPAC to regular expresions.

    Decodes a sequence which is IUPAC and convert
    this to a regular expression friendly sequence.
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
    # returns the sequence
    return seq


def to_list(seq: str) -> (list):
    """Convert a sequence to a list."""
    retval = []
    for base in seq:
        if base == "N":
            retval.append("")
        elif base == "V":
            retval.append("ACG")
        elif base == "H":
            retval.append("ACT")
        elif base == "D":
            retval.append("AGT")
        elif base == "B":
            retval.append("CGT")
        elif base == "M":
            retval.append("AC")
        elif base == "K":
            retval.append("GT")
        elif base == "S":
            retval.append("GC")
        elif base == "W":
            retval.append("AT")
        elif base == "R":
            retval.append("AG")
        elif base == "Y":
            retval.append("CT")
        else:
            retval.append(base)
    return retval


def match(sequence: str, recognition: list) -> (bool):
    """Match the sequence with a recognition sequence."""
    if len(recognition) > len(sequence):
        return False
    for idx, bases in enumerate(recognition):
        if bases:
            if not sequence[idx] in bases:
                return False
    return True


