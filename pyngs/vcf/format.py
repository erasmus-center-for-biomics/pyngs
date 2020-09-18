from typing import TypeVar
from .values import ParseVcfTypes, ReprVcfTypes, VcfValue, VCFTYPES, VCFNUMBERS
from .meta import Meta


F = TypeVar("F")

class Format:

    def __init__(self, code: str, vtype: str, number: str) -> None:
        """Initialize a new Format object."""
        # check the type
        if vtype not in VCFTYPES:
            raise ValueError("invalid value type {0} for format code {1}".format(vtype, code))

        # check the value number
        self.numval = -1
        if number not in VCFNUMBERS:
            try:
                self.numval = int(number)
            except ValueError:
                raise ValueError("invalid value number {0} for format code {1}".format(number, code))
        self.code = code
        self.type = vtype
        self.number = number
        self.valueparser = ParseVcfTypes[self.type]
        self.representer = ReprVcfTypes[self.type]

    @classmethod
    def from_meta(cls: F, meta: Meta) -> F:
        """Create a Format line from a format meta object."""
        if meta.key.lower() != "format":
            raise ValueError("key is {0} not a format field".format(meta.key))

        mtype = None
        mnumber = None
        for pair in meta.value.content:
            if pair[0] == "Type":
                mtype = pair[1]
            if pair[0] == "Number":
                mnumber = pair[1]

        return cls(
            meta.value.id,
            mtype,
            mnumber)

