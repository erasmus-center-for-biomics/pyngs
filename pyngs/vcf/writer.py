from typing import TextIO, List
from .meta import Meta
from .variant import Variant


class Writer:

    def __init__(self, outstream: TextIO, metas: List[Meta], samples: List[str]) -> None:
        super().__init__()
        self.outstream = outstream
        for meta in metas:
            self.outstream.write(str(meta))

        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        if samples:
            header.append("FORMAT")
            header.extend(samples)
        outstream.write("\t".join(header) + "\n")

    def write(self, variant: Variant) -> None:
        """Write the variant to the output stream."""
        self.outstream.write(str(variant))
