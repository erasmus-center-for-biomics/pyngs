from typing import TextIO, List
from .meta import Meta
from .variant import Variant


class Writer:

    def __init__(self, outstream: TextIO, metas: List[Meta], samples: List[str]) -> None:
        super().__init__()
        self.outstream = outstream
        for meta in meta:
            self.outstream.write(meta)

        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        if samples:
            header.append("FORMAT")
            header.extend(samples)

    def write(self, variant: Variant) -> None:
        """Write the variant to the output stream."""
        self.outstream.write(str(variant))
