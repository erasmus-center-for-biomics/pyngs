from typing import TextIO, List, Iterator, List
from .header import Header
from .row import Row, HeaderIndex, Variant
from .utils import quote_tokenizer, genotypes


class Reader:

    def __init__(self, instream: TextIO) -> None:
        """Initialize the VCF reader object."""
        self.header = []
        self.samples = []
        self.stream = instream
        self.info = {}
        self.format = {}
        self.filter = {}
        self.version = None
        self.__parse_header__()

    def __parse_header__(self) -> None:
        """Parse the header of the VCF file."""
        for line in self.stream:
            line = line.rstrip()
            if line.startswith("##"):
                entry = Header(line[2:])
                if entry.section == "fileformat":
                    self.version = entry.value
                self.header.append(entry)
            elif line.startswith("#"):
                tokens = line.split("\t")
                if len(tokens) > 9:
                    self.samples = tokens[9:]
                break

        # index the parser
        for idx, entry in enumerate(self.header):
            if entry.section == "FORMAT" and entry.id:
                self.format[entry.id] = idx
            elif entry.section == "FILTER" and entry.id:
                self.filter[entry.id] = idx
            elif entry.section == "INFO" and entry.id:
                self.info[entry.id] = idx

    def field(self, section: str="INFO", name: str="") -> Header:
        """Get the relevant VCF header for parsing."""
        if section == "FORMAT":
            return self.header[self.format[name]]
        elif section == "INFO":
            return self.header[self.info[name]]
        elif section == "FILTER":
            return self.header[self.filter[name]]
        raise KeyError("section '{0}' is not registered".format(section))

    def has_id(self, name: str, section: str=None) -> bool:
        """Check whether ID is already present in the header."""
        for entry in self.header:
            if entry.id == name:
                if section is None:
                    return True
                elif section == section:
                    return True
        return False

    def __iter__(self) -> Iterator[Row]:
        """Return a new iterator over the VCF file."""
        for line in self.stream:
            line = line.rstrip()
            if not line:
                continue
            tokens = line.split("\t")
            yield Row(tokens)


class Writer:
    """A class to write a VCF file."""
    def __init__(self, stream: TextIO, header: List[Header], samples: List[str]) -> None:
        """Initialize a VCF writer."""
        self.stream = stream
        self.header = header
        for definition in header:
            self.stream.write("##{0}\n".format(repr(definition)))
        preset = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        if samples:
            preset.append("FORMAT")
            preset.extend(samples)
        self.stream.write(
            "{0}\n".format(
                "\t".join(preset)))

    def write(self, row: Row) -> None:
        """Write a VCF row."""
        self.stream.write("{0}\n".format(repr(row)))
