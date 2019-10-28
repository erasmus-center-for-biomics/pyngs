"""
A generator to parse FastQ files.
"""
import sys
from typing import TextIO, Generator, Tuple, Iterable


def fastq(instream: TextIO) -> Generator[Tuple[str, str, str], None, None]:
    """
    Parse a FastQ file.
    :param instream: the input stream with FastQ formatted data.
    :return: a generator with the read name, sequence and quality string
    """
    try:
        while True:
            name = next(instream).rstrip()[1:]
            sequence = next(instream).rstrip()
            next(instream).rstrip()
            quality = next(instream).rstrip()
            yield name, sequence, quality
    except StopIteration:
        pass


def format(name: str, sequence: str, quality: str) -> str:
    """
    Format a FastQ entry.
    :param name: the read name
    :param sequence: the read sequence
    :param quality: the read quality
    :return: a formatted fastq entry
    """
    return "@{name}\n{seq}\n+\n{qual}\n".format(
        name=name,
        seq=sequence,
        qual=quality)


def clean_readname(readname: str) -> str:
    """
    Removes all characters after the first space.
    :param readname: the name of a read
    :return: the string of the readname up to the first space
    """
    return readname.split(" ")[0]


def encode_in_readname(readname: str, data: Iterable[str], sep: str=":") -> str:
    """
    Encode data in the read-name.
    :param readname: the name of a read
    :param data: a list of data to encode in the read name
    :param sep: the field separator in the read name, default is :
    :return: a new readname with the data encode
    """
    return readname + sep + sep.join([str(d) for d in data])
