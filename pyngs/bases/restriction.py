"""
A package to deal with restriction data.

Restriction enzyme recognition sites recognize a
specific sequence and cut on the leading (/) and/or lagging (|)
strand of the DNA. In some cases the cut-site is downstream
from the recognized sequence.

A special class of restriction enzymes specifically cuts
methylated DNA. The methylated base is indicated with a ^
character before it.

"""
from .iupac import complement, expand_iupac, IUPAC_BASES

"""
the possible sequence modifiers that we recognize
^ next base is methylated
/ cut on the leading strand
| cut on the lagging strand
"""
SEQUENCE_MODIFIERS = ["^", "/", "|"]


def complement_modifiers(char: str):
    """Translate the additional sequence modifiers."""
    if char == "/":
        return "|"
    elif char == "|":
        return "/"
    else:
        return char


def over_recognition(sequence: str):
    """Move over the sequence and pick up any
    modifiers that are present."""
    offset = 0
    for char in sequence:
        yield offset, char
        if char not in SEQUENCE_MODIFIERS:
            offset += 1


def partition(sequence: str):
    """Partition the sequence in SEQUENCES and MODIFIERS."""
    current = []
    previous = None
    for base in sequence:
        if base not in IUPAC_BASES:
            ctype = "__SEQUENCE__"
        else:
            ctype = "__MODIFIER__"

        if not previous:
            previous = ctype
        if ctype != previous:
            yield previous, current
            current = []
        current.append(base)
        previous = ctype
    yield previous, current


class Enzyme:
    """
    A representation of a restriction enzyme.

    As cuts occur between bases, we record 2 positions
    between which DNA is cut.
    """

    def __init__(self, name: str="", sequence: str=""):
        self.name = name
        self.sequence = sequence
        self.bases = []
        self.recognition = []
        self.sites = []
        self.marks = []
        self.__interpret__()

    def __interpret__(self):
        for idx, base in over_recognition(self.sequence):
            if base in "/|":
                self.sites.append((idx-1, idx, base))
            elif base in IUPAC_BASES:
                self.bases.append(base)
                self.recognition.append(
                    expand_iupac(base))
            else:
                self.marks.append((idx, base))

    def reverse_complement(self):
        """Get the reverse complement of the enzyme."""
        bases = [complement(b) for b in self.bases]
        bases.reverse()
        recognition = [expand_iupac(b) for b in bases]
        sites = []
        for site in self.sites:
            sites.append((
                len(bases) - site[1] - 1,
                len(bases) - site[0] - 1,
                complement_modifiers(site[2])
                ))
        sites.reverse()
        marks = []
        for mark in self.marks:
            marks.append(
                (len(bases) - mark[0] - 1, mark[1]))
        marks.reverse()

        enzyme = Enzyme()
        enzyme.bases = bases
        enzyme.recognition = recognition
        enzyme.sites = sites
        enzyme.marks = marks
        return enzyme

    def to_sequence(self):
        """Generate a sequence from the enzyme."""
        bases = [b for b in self.bases]
        for site in self.sites:
            if site[0] >= 0:
                bases[site[0]] += site[2]
            else:
                bases[site[1]] = site[2] + bases[site[1]]
        for mark in self.marks:
            bases[mark[0]] = mark[1] + bases[mark[0]]
        return "".join(bases)

    def leading_site(self):
        for site in self.sites:
            if site[2] == "/":
                yield site

    def lagging_site(self):
        for site in self.sites:
            if site[2] == "|":
                yield site

    def __repr__(self):

        marks = ["-"] * len(self.bases)
        for mrk in self.marks:
            marks[mrk[0]] = mrk[1]

        sites = ["-"] * len(self.bases)
        for site in self.sites:
            if site[1] >= 0 and site[1] < len(sites):
                sites[site[1]] = "<"
            if site[0] >= 0 and site[0] < len(sites):
                sites[site[0]] = ">"

        return "sequence: {0}\nsites:    {1}\nmarks:    {2}\n".format(
            "".join(self.bases),
            "".join(sites),
            "".join(marks))
