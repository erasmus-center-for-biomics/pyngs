import logging
from typing import List, Union, Tuple, Optional
from .utils import quote_tokenizer, genotypes
from .header import Header

class Row:

    # __slots__ = [
    #     "chrom", "position", "id",
    #     "reference", "alternate",
    #     "quality", "filter",
    #     "info", "format", "samples"]

    def __init__(self, tokens: List[Union[str, int]]) -> None:
        """Initialize a new VCF row."""
        self.chrom: str = tokens[0]
        self.position: int = tokens[1]
        self.id: str = tokens[2]
        self.reference: str = tokens[3]
        self.alternate: str = tokens[4]
        self.quality: float = tokens[5]

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
        self.format: List[str] = None
        if len(tokens) > 8:
            self.format = tokens[8].split(":")

        # parse the samples
        self.samples: List[str] = []
        if len(tokens) > 9:
            for sidx, token in enumerate(tokens[9:]):
                spart = [x for x in quote_tokenizer(token, ":")]

                # check the format partitions
                if len(spart) != len(self.format):
                    logging.warn(
                        "For variant %s:%d sample %d: Expected %d fields, found %d",
                        self.chrom, self.position, sidx,
                        len(self.format), len(spart))
                    if len(self.format) > len(spart):
                        for _ in range(len(self.format) - len(spart)):
                            spart.append(".")
                    elif len(self.format) < len(spart):
                        spart = spart[0:len(self.format)]

                # add the sample
                self.samples.append(spart)

    @property
    def position(self) -> Union[None, int]:
        return self.__pos__

    @position.setter
    def position(self, value: Union[str, int]) -> None:
        self.__position__ = value
        try:
            self.__pos__ = int(value)
        except ValueError:
            self.__pos__ = None

    @property
    def quality(self) -> Union[None, int]:
        return self.__qual__

    @quality.setter
    def quality(self, value: Union[str, int]) -> None:
        self.__quality__ = value
        try:
            self.__qual__ = int(self.__quality__)
        except ValueError:
            self.__qual__ = None

    def add_data(self, fcode: str, data: List[Union[None, str, float, int]]) -> None:
        """Adds new data to the format and samples of a variant."""
        assert(len(data) == len(self.samples))

        # add the code to the format
        self.format.append(fcode)

        # add the data
        for sidx, value in enumerate(data):
            if value is None:
                value = "."
            self.samples[sidx].append(value)


    def __repr__(self) -> str:
        """return a string representation of the data."""
        info = []
        for entry in self.info:
            if entry[1] is True:
                info.append(entry[0])
            elif entry[1] is None or entry[1] is False:
                pass
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
            sval = ":".join([str(s) for s in sample])
            vstr += "\t{0}".format(sval)
        #
        return vstr

class HeaderIndex:

    def __init__(self, headers: List[Header]) -> None:
        """."""
        self.headers = headers
        self.info_index = {}
        self.format_index = {}

        for idx, header in enumerate(self.headers):
            if header.section == "INFO":
                self.info_index[header.id] = idx
            if header.section == "FORMAT":
                self.format_index[header.id] = idx

    def get(self, section: str, id: str) -> Header:
        if section == "INFO":
            return self.headers[self.info_index[id]]
        elif section == "FORMAT":
            return self.headers[self.format_index[id]]


# possible values in a VCF
VcfValue = Union[Optional[str], Optional[float], Optional[int], Optional[bool], List[Optional[str]], List[Optional[float]], List[Optional[int]], List[Optional[bool]]]

# a dict with VCF type parsers
ParseVcfTypes = {
    "String": lambda x: x,
    "Character": lambda x: x[0],
    "Integer": lambda x: int(x),
    "Float": lambda x: float(x),
    "Flag": lambda x: True
}

def repr_str(x):
    if x is None:
        return "."
    elif "," in x:
        return '"' + x + '"'
    else:
        return x

# A dict with VCF type writers
ReprVcfTypes = {
    "String": repr_str,
    "Character": lambda x: x[0] if x is not None else ".",
    "Integer": lambda x: str(x) if x is not None else ".",
    "Float": lambda x: str(x) if x is not None else ".",
    "Flag": lambda x: ""
}

class Variant:

    def __init__(self,
            chrom: str, pos: int, ids: str, ref: str, alt: List[Optional[str]],
            qual: Optional[float], filters: List[Optional[str]],
            info: List[Tuple[str, str, str, VcfValue]],
            formats: Optional[List[Tuple[str, str, str]]], samples: Optional[List[List[VcfValue]]]) -> None:
        self.chrom = chrom
        self.position = pos
        self.id = ids
        self.reference = ref
        self.alternate = alt
        self.quality = qual
        self.filter = filters
        self.info = info
        self.format = formats
        self.samples = samples

    @classmethod
    def from_row(cls, row: Row, index: HeaderIndex):

        # parse the alternate alleles
        alternates = [a for a in quote_tokenizer(row.alternate, ",")]

        # parse the info fields
        info: List[Tuple[str, str, VcfValue]] = []
        for key, val in row.info:
            value = None
            pa = index.get("INFO", key)
            if pa.type == "Flag":
                info.append((key, pa.type, pa.number, []))
                continue

            value = []
            for v in quote_tokenizer(val, ","):
                value.append(ParseVcfTypes[pa.type](v))
            info.append((key, pa.type, pa.number, value))

        # parse the sample and format fields
        formats = None
        samples = None
        if row.format is not None:
            formats = []
            samples = [[None] * len(row.format)] * len(row.samples)
            for fidx, key in enumerate(row.format):
                pa  = index.get("FORMAT", key)
                func = ParseVcfTypes[pa.type]
                formats.append((key, pa.type, pa.number))

                # parse the values per sample
                for sidx, sample in enumerate(row.samples):
                    value = []
                    for v in quote_tokenizer(sample[fidx], ","):
                        value.append(func(v))
                    samples[sidx][fidx] = value
        # return a variant
        return Variant(
                row.chrom,
                row.position,
                row.id,
                row.reference,
                alternates,
                row.quality,
                row.filter,
                info,
                formats,
                samples)

    def filter_alternate(self, altidx: int, ploidy: int=2):
        """Return a variant with alternate index altidx removed."""

        gt_idx_to_rm = []
        for idx, geno in enumerate(genotypes(ploidy, len(self.alternate))):
            if altidx in geno:
                gt_idx_to_rm.append(idx)

        alternates = [a for i,a in enumerate(self.alternate) if i != altidx]

        # filter the info
        info = []
        for entry in self.info:
            if entry[1] == "Flag":
                info.append(entry)
                continue

            values = []

            if entry[2] == "G":
                for idx, val in enumerate(entry[3]):
                    if idx not in gt_idx_to_rm:
                        values.append(val)

            elif entry[2] == "R":
                for idx, val in enumerate(entry[3]):
                    if idx != altidx + 1:
                        values.append(val)

            elif entry[2] == "A":
                for idx, val in enumerate(entry[3]):
                    if idx != altidx:
                        values.append(val)

            else:
                values = entry[3]
            info.append((entry[0], entry[1], entry[2], values))

        # filter the samples
        samples = [[None] * len(self.format)] * len(self.samples)
        for fidx, ftype in enumerate(self.format):
            for sidx, svalue in enumerate(self.samples):
                values = []
                if ftype[2] == "G":

                    for idx, val in enumerate(svalue[fidx]):
                        if idx not in gt_idx_to_rm:
                            values.append(val)

                elif ftype[2] == "R":
                    for idx, val in enumerate(svalue[fidx]):
                        if idx  != altidx + 1:
                            values.append(val)
                    samples[sidx][fidx] = values

                elif ftype[2] == "A":
                    for idx, val in enumerate(svalue[fidx]):
                        if idx  != altidx:
                            values.append(val)
                    samples[sidx][fidx] = values
                else:
                    values = svalue[fidx]

                samples[sidx][fidx] = values

        # return the modified variant
        return Variant(
            self.chrom,
            self.position,
            self.id,
            self.reference,
            alternates,
            self.quality,
            self.filter,
            info,
            self.format,
            samples)

    def __repr__(self):
        """Return a string representation of this variant."""

        info = []
        for entry in self.info:
            if entry[1] == "Flag":
                info.append(entry[0])
                continue

            func = ReprVcfTypes[entry[1]]
            entry = "{0}={1}".format(
                entry[0],
                ",".join([func(v) for v in entry[3]]))
            info.append(entry)

        alt = "."
        if self.alternate:
            alt = ",".join(self.alternate)
        #
        bstr = "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}"
        vstr = bstr.format(
            chrom=self.chrom,
            pos=self.position,
            id=self.id if self.id is not None else ".",
            ref=self.reference,
            alt=alt,
            qual=self.quality if self.quality is not None else ".",
            filter=";".join(self.filter),
            info=";".join(info))

        if self.format is not None:
            funcs = []
            for fentry in self.format:
                funcs.append(ReprVcfTypes[fentry[1]])

            formatstr = ":".join([f[0] for f in self.format])
            vstr += "\t{0}".format(formatstr)

            for sample in self.samples:
                data = []
                for idx, fentry in enumerate(self.format):
                    val = ",".join([funcs[idx](v) for v in sample[idx]])
                    data.append(val)
                vstr += "\t{0}".format(":".join(data))
        #
        return vstr

    # def guess_ploidy(self):
    #     """Guess the ploidy based on a the length of a genotype info field."""
    #     total
