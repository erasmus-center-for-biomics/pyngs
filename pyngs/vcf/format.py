from .utils import VCFTYPES, VCFNUMBERS
from .values import ParseVcfTypes, ReprVcfTypes, VcfValue
from .meta import Meta

class Format:

    def __init__(self, code: str, vtype: [str], number: str) -> None:
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
    def from_header(cls, header: Header) -> Format:
        """Create a Format object from a Header."""
        return Format(header.id, header.type, header.number)

    @classmethod
    def from_meta(cls, meta: Meta) -> Format:
        """Create a Format line from a format meta object."""
        if meta.key != "format":
            raise ValueError("key is {0} not format".format(meta.key))

        return Format(
            meta.value.id,
            meta.value.content["Type"],
            meta.value.content["Number"])

