from typing import List, Dict, Callable

from .values import VcfValue
from .format import Format
from .info import Info
from .utils import quote_tokenizer


class ValueStore:

    def __init__(self) -> None:
        super().__init__()
        self.storage: Dict[str, VcfValue] = {}

    def add(self, code: str, value: VcfValue) -> None:
        """Add or set a value to the per store."""
        self.storage[code] = value

    def add_str(self, code: str, value: str, converter: Callable[[str], VcfValue]) -> None:
        """Add or set a value based on a str."""
        self.add(code, converter(value))

    def get(self, code: str) -> VcfValue:
        """Get the value for field code."""
        try:
            return self.storage[code]
        except KeyError:
            return None

    def to_format_str(self, formats: List[Format]) -> str:
        """Convert the stored values to a valid format string."""
        values = []
        for fmt in formats:
            value = self.get(fmt.code)
            values.append(fmt.representer(value))
        return ":".join(values)

    def to_info_str(self, info: List[Info]):
        """Convert the stored values to a valid info string."""
        parts = []
        for fld in info:
            value = self.get(fld.code)
            if fld.type == "Flag":
                if value:
                    parts.append(fld.code)
            else:
                strval = fld.representer(value)
                parts.append("{-}={1}".format(fld.code, strval))
        return ";".join(parts)

    def __bool__(self) -> bool:
        """Check if the store is empty."""
        return bool(self.storage)

class FormatStore:
    """A store to represent and hold the values for samples."""

    def __init__(self) -> None:
        """Initialize the object."""
        self.format: List[Format] = []
        self.stores: List[ValueStore] = []
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
            self.stores.append(ValueStore())

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
            self.stores[sidx].add_str(fmt.code, data[idx], fmt.valueparser)

    def __str__(self) -> str:
        """Represent the data in the store as a string."""
        fmtbase = ":".join([fmt.code for fmt in self.format])
        smpstrs = [st.to_format_str() for st in self.stores]
        return "{0}\t{1}".format(fmtbase, "\t".join(smpstrs))


class InfoStore:

    def __init__(self) -> None:
        super().__init__()
        self.info: Dict[str, Info] = {}
        self.store: ValueStore = ValueStore()

    def __bool__(self) -> bool:
        """Check if the store is empty."""
        return bool(self.store)

    def add_info(self, info: Info) -> None:
        """Add an info parser to the obejct."""
        self.info[info.code] = info

    def add_data(self, strval: str) -> None:
        """Add data to the store."""
        if strval == ".":
            return
        parts = [p for p in quote_tokenizer(strval, sep=";")]
        for part in parts:
            if "=" in part:
                [code, value] = part.split("=", 1)
                self.store.add_str(code, value, self.info[code].valueparser)
            else:
                self.store.add_str(code, ".", self.info[code].valueparser)

    def __str__(self) -> str:
        """Convert the data to a str."""
        if not self.store.storage:
            return "."
        return self.store.to_info_str(self.info)
