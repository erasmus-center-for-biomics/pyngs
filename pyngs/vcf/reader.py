from collections import OrderedDict
from typing import TextIO, List, Dict, Generator

from .variant import Variant
from .meta import Meta
from .info import Info
from .format import Format


class Reader:

    def __init__(self, instream: TextIO) -> None:
        """Initialize the reader object."""
        self.instream = instream

        self.meta: List[Meta] = []
        self.samples: List[str] = []
        for line in self.instream:
            if line.startswith("##"):
                self.meta.append(Meta.from_str(line))
            elif line.startswith("#"):
                tokens = line.rstrip().split("\t")
                if len(tokens) > 9:
                    self.samples = tokens[9:]
                break
        # add the help dictionaries for the variants
        self.info: OrderedDict[str, Info] = OrderedDict()
        self.infokeys: Dict[str, int] = {}
        self.format: Dict[str, Format] = {}
        for meta in self.meta:
            if meta.key.lower() == "info":
                obj = Info.from_meta(meta)
                self.info[obj.code] = obj
            elif meta.key.lower() == "format":
                obj = Format.from_meta(meta)
                self.format[obj.code] = obj


    def __iter__(self) -> Generator[Variant, None, None]:
        """Iterate over the variants in a VCF file."""
        for line in self.instream:
            line = line.rstrip()
            if not line:
                continue
            yield Variant.from_str(
                line,
                self.info,
                self.format)
