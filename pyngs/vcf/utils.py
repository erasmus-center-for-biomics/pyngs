from typing import Generator, List

def quote_tokenizer(sstring: str="", sep: str=",", quote='"') -> Generator[str, None, None]:
    """Separate the fields while ignoring quotes."""
    tokens = sstring.split(sep)
    keep = False
    current = []
    for token in tokens:
        if keep:
            current.append(token)
            if token[-1] == quote:
                keep = False
                yield sep.join(current)
        elif quote in token:
            # if the quote is ended in the same token,
            # just yield the token
            if token[-1] == quote:
                yield token
            # otherwise accumulate until the ending quote.
            else:
                keep = True
                current = [token]
        # not quotes, no issue
        else:
            yield token


def genotypes(ploidy: int=2, alleles: int=1, thusfar: List[int]=None) -> Generator[List[int], None, None]:
    """
    Get the genotype compbinations according to https://samtools.github.io/hts-specs/VCFv4.3.pdf.

    Note that alleles is equivalent to N and indicates the largest index of the alleles present.
    """

    # initialize the allele
    if thusfar is None:
        thusfar = []

    for allele in range(0, alleles + 1):
        #
        if ploidy == 1:
            yield [allele] + thusfar

        elif ploidy > 1:
            yield from genotypes(ploidy-1, allele, [allele] + thusfar)

