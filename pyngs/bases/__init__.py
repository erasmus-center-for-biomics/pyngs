
import typing


def fasta(instream: typing.TextIO, toupper: bool=True, fullnames: bool=False):
    """
    Get the name, sequence pairs from a FastA formatted file.

    :param instream: the input file stream
    :param toupper: convert the sequence to uppercase
    :param fullnames: keep the full names
    :yield: a tuple with the name and the sequence
    """
    name = None
    sequence = []
    for line in instream:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, "".join(sequence)
            name = line[1:].split(" ")[0] if not fullnames else line[1:]
            sequence = []
        else:
            if toupper:
                sequence.append(line.upper())
            else:
                sequence.append(line)
    if name:
        yield name, "".join(sequence)
