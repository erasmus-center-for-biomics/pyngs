import typing
from .utils import quote_tokenizer


class Row:

    # __slots__ = [
    #     "chrom", "position", "id",
    #     "reference", "alternate",
    #     "quality", "filter",
    #     "info", "format", "samples"]

    def __init__(self, tokens: list):
        """Initialize a new VCF row."""
        self.chrom = tokens[0]
        self.position = tokens[1]
        self.id = tokens[2]
        self.reference = tokens[3]
        self.alternate = tokens[4]
        self.quality = tokens[5]

        # tokenize the filter
        self.filter = [x for x in quote_tokenizer(tokens[6], ";")]
        self.info = []
        for token in quote_tokenizer(tokens[7], ";"):
            vals = token.split("=", 1)
            if len(vals) == 2:
                key, value = vals[0], vals[1]
            elif len(vals) == 1:
                key, value = vals[0], True
            self.info.append((key, value))

        # parse the format
        self.format = None
        if len(tokens) > 8:
            self.format = tokens[8].split(":")

        # parse the samples
        self.samples = []
        if len(tokens) > 9:
            for token in tokens[9:]:
                spart = [x for x in quote_tokenizer(token, ":")]
                self.samples.append(spart)
                if len(spart) != len(self.format):
                    raise ValueError("Expected {0} fields, found {1}".format(
                        len(self.format), len(spart)))

    @property
    def position(self):
        return self.__pos__

    @position.setter
    def position(self, value):
        self.__position__ = value
        try:
            self.__pos__ = int(value)
        except ValueError:
            self.__pos__ = None

    @property
    def quality(self):
        return self.__qual__

    @quality.setter
    def quality(self, value):
        self.__quality__ = value
        try:
            self.__qual__ = int(self.__quality__)
        except ValueError:
            self.__qual__ = None

    def add_data(self, fcode: str, data: typing.List):
        """Adds new data to the format and samples of a variant."""
        assert(len(data) == len(self.samples))

        # add the code to the format
        self.format.append(fcode)

        # add the data
        for sidx, value in enumerate(data):
            if value is None:
                value = "."
            self.samples[sidx].append(repr(value))


    def __repr__(self):
        """return a string representation of the data."""
        info = []
        for entry in self.info:
            if entry[1] is True:
                info.append(entry[0])
            else:
                info.append("{0}={1}".format(entry[0], entry[1]))
        #
        bstr = "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}"
        vstr = bstr.format(
            chrom=self.chrom,
            pos=self.__position__,
            id=self.id,
            ref=self.reference,
            alt=self.alternate,
            qual=self.__quality__,
            filter=";".join(self.filter),
            info=";".join(info))
        if self.format is not None:
            vstr += "\t{0}".format(":".join(self.format))
        for sample in self.samples:
            vstr += "\t{0}".format(":".join(sample))
        #
        return vstr
