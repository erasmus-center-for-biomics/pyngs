from typing import List, Union, Optional, Callable, Dict
from .utils import quote_tokenizer

# the types in a VCF file
VCFTYPES = ("Flag", "Character", "String", "Float", "Integer")
VCFNUMBERS = (".", "G", "R", "A")

# possible python types from a VCF
VcfValue = Union[None, bool, List[Optional[str]], List[Optional[float]], List[Optional[int]]]
VcfPrintable = Union[None, List[Optional[str]], List[Optional[float]], List[Optional[int]]]
VcfBaseValues = Union[None, str, float, int]

def parse_types(convert: Callable[[str], Union[bool, str, float, int]]) -> Callable[[str], VcfPrintable]:
    """Parse types from a string to a list."""
    def inner(value: str) -> VcfPrintable:
        if value == ".":
            return None
        try:
            retval = []
            for s in quote_tokenizer(value, sep = ","):
                if s is ".":
                    retval.append(None)
                else:
                    retval.append(convert(s))
            return retval
        except ValueError:
            return None
    #
    return inner


# a dict with VCF type parsers
ParseVcfTypes: Dict[str, Callable[[str], VcfValue]] = {
    "String": parse_types(str),
    "Character": parse_types(lambda x: x[0]),
    "Integer": parse_types(int),
    "Float": parse_types(float),
    "Flag": lambda x: True
}


def to_str_types(convert: Callable[[VcfBaseValues], str]) -> Callable[[VcfValue], str]:
    """Represent data types for VCF."""
    def inner(values: VcfValue):
        if values is None:
            return "."

        fields = []
        for s in values:
            if s is None:
                fields.append(".")
            else:
                fields.append(convert(s))
        return ",".join(fields)

    return inner

# A dict with VCF type writers
ReprVcfTypes = {
    "String": to_str_types(str),
    "Character": to_str_types(lambda x: x[0]),
    "Integer": to_str_types(str),
    "Float": to_str_types(str),
    "Flag": lambda x: ""
}

