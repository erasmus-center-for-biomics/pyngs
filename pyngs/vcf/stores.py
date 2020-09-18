from collections import OrderedDict
from typing import List, Dict, Callable

from .values import VcfValue
from .format import Format
from .info import Info
from .utils import quote_tokenizer


class FormatStore:
    """A store to represent and hold the values for samples."""

    def __init__(self) -> None:
        """Initialize the object."""
        self.format: List[Format] = []
        self.stores: List[Dict[str, VcfValue]] = []
        self.__changed__ = False
        self.__tags__ = set()

    def __bool__(self) -> bool:
        """Check if the store is empty."""
        return bool(self.format)

    def has_tag(self, code: str) -> bool:
        """Is the tag defined in the current store."""
        # determine the tag set only once
        if self.__changed__:
            self.__tags__ = set([f.code for f in self.format])
            self.__changed__ = False
        return code in self.__tags__

    def add_format(self, fmt: Format) -> None:
        """Add a format to the store."""
        # append to the Format list
        self.format.append(fmt)
        self.__changed__ = True

    def add_formats_from_dict(self, strval: str, fmt_d: Dict[str,Format]) -> None:
        """Add a format to the store."""
        # append to the Format list
        for code in strval.split(":"):
            self.add_format(fmt_d[code])

    def init_stores(self, n_samples: int) -> None:
        """Initialize stores for n samples."""
        self.stores = []
        for _ in range(0, n_samples):
            self.stores.append({})

    def add_sample_data(self, sidx: int, strdata: str) -> None:
        """Add data for a sample in from strings."""

        # check whether the sample is present
        if sidx >= len(self.stores):
            raise ValueError("only {0} samples defined not {1}", len(self.stores), sidx)

        # check the number of fields
        data = strdata.split(":")
        if len(data) != len(self.format):
            raise ValueError("expeccted {0} fields not {1}", len(self.format), len(data))

        # add the data to the store
        for idx, fmt in enumerate(self.format):
            self.stores[sidx][fmt.code] = fmt.valueparser(data[idx])

    def __sample_str__(self, store) -> str:
        """Convert the stored values to a valid format string."""
        values = []
        for fmt in self.format:
            value = store[fmt.code]
            values.append(fmt.representer(value))
        return ":".join(values)

    def __str__(self) -> str:
        """Represent the data in the store as a string."""
        fmtbase = ":".join([fmt.code for fmt in self.format])
        smpstrs = [self.__sample_str__(st) for st in self.stores]
        return "{0}\t{1}".format(fmtbase, "\t".join(smpstrs))


class InfoStore:

    def __init__(self) -> None:
        super().__init__()
        self.info: Dict[str, Info] = OrderedDict()
        self.data: Dict[str, VcfValue] = {}

    def __bool__(self) -> bool:
        """Check if the store is empty."""
        return bool(self.data)

    def add_info(self, info_p: Info) -> None:
        """Add an info parser to the obejct."""
        self.info[info_p.code] = info_p

    def add_data(self, strval: str, info_p: Dict[str, Info]) -> None:
        """Add data to the store."""
        if strval == ".":
            return
        parts = [p for p in quote_tokenizer(strval, sep=";")]
        for part in parts:
            if "=" in part:
                [code, value] = part.split("=", 1)
                self.add_info(info_p[code])
                self.data[code] = self.info[code].valueparser(value)
            else:
                code = part
                self.add_info(info_p[code])
                self.data[code] = True

    def __str__(self) -> str:
        """Convert the data to a str."""
        parts = []
        for code, info_p in self.info.items():
            if not code in self.data.keys():
                continue
            if info_p.type == "Flag":
                parts.append(code)
            else:
                strval = info_p.representer(self.data[code])
                parts.append("{0}={1}".format(code, strval))
        return ";".join(parts)
