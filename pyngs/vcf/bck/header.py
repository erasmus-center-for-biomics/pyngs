from .utils import quote_tokenizer
from typing import Generator, Union

class Header:

    def __init__(self, hstring: str=None):
        """Initialize a new VCF header."""
        self.section = None
        self.value = None
        self.id = None
        # G: genotypes, A: alternate alleles,
        # R: reference and alternate alleles,
        # .: unbound,
        # [0-9]+: specific number
        self.number = None
        self.real_number = None

        # Integer, Float, Flag, Character, and String
        self.type = None
        self.description = None
        self.other = {}
        if hstring:
            self.__parse__(hstring)

    def __parse__(self, hstring: str="") -> None:
        """Parse the header line."""
        tokens = hstring.split("=", 1)
        self.section = tokens[0]
        self.value = tokens[1]
        if self.value.startswith("<"):
            self.__parse_value__(self.value.lstrip("<").rstrip(">"))

    def __parse_value__(self, value: str="") -> None:
        """Parse the value."""
        for token in quote_tokenizer(value, ","):
            fields = token.split("=", 1)
            if fields[0].lower() == "id":
                self.id = fields[1]
            elif fields[0].lower() == "number":
                self.number = fields[1]
                if self.number not in ("G", "A", "R", "."):
                    self.real_number = int(self.number)
            elif fields[0].lower() == "type":
                self.type = fields[1]
            elif fields[0].lower() == "description":
                self.description = fields[1]
            else:
                self.other[fields[0]] = fields[1]

    def interpret(self, value: str="") -> Generator[Union[None, str, float,int], None, None]:
        """Parse a field with the specified data."""
        if self.type == "Flag":
            yield True
        # get the field separator
        sep = ","
        maxidx = 0
        for token in quote_tokenizer(value, sep):
            maxidx += 1
            if self.type == "String":
                yield token
            elif self.type == "Character":
                yield token[0]
            elif self.type == "Float":
                yield float(token) if token != "." else None
            elif self.type == "Integer":
                yield int(token) if token != "." else None

        # check specific numbers
        if self.real_number is not None and maxidx != self.real_number:
            raise ValueError("expected {0} entries obtained {1}".format(
                self.real_number, maxidx))

    def __repr__(self):
        """Return a string representation of the VCF header."""
        tokens = []
        if self.id:
            tokens.append("ID={0}".format(self.id))
        if self.number:
            tokens.append("Number={0}".format(self.number))
        if self.type:
            tokens.append("Type={0}".format(self.type))
        if self.description:
            tokens.append("Description={0}".format(self.description))
        for key, value in self.other.items():
            tokens.append("{0}={1}".format(key, value))
        value = self.value
        if tokens:
            value = "<{0}>".format(",".join(tokens))
        return "{0}={1}".format(self.section, value)
