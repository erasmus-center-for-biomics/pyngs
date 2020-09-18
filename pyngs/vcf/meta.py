from typing import Dict, Tuple, List, TypeVar, Union
from .utils import quote_tokenizer

S = TypeVar("S")

class Structured:

    def __init__(self, idval: str, content: List[Tuple[str,str]]) -> None:
        """Initialize a new Structured object."""
        self.id = idval
        self.content = content

    @classmethod
    def from_str(cls: S, strval: str) -> S:
        """Get structured information from a string."""
        parts = [p for p in quote_tokenizer(strval, ",")]

        idval = None
        content: List[Tuple[str,str]] = []

        for idx, part in enumerate(parts):
            tokens = part.split("=", 1)
            if tokens[0] == "ID":
                idval = tokens[1]
            else:
                content.append((tokens[0], tokens[1]))
        return cls(idval, content)

    def __str__(self) -> str:
        """Represent the data as a string."""
        parts = ["{0}={1}".format(v[0], v[1]) for v in self.content]
        return "<ID={0},{1}>".format(self.id, ",".join(parts))



M = TypeVar("M")

class Meta:

    def __init__(self, key: str, value: Union[str, Structured]) -> None:
        """Initialize a new Meta object."""
        super().__init__()
        self.key = key
        self.value = value

    @classmethod
    def from_str(cls: M, strval: str) -> M:
        """Get a meta line from a string."""
        [key, value] = strval.rstrip()[2:].split("=", 1)

        if value.startswith("<"):
            if value.endswith(">"):
                value = value.rstrip(">")
            value = Structured.from_str(value.lstrip("<"))

        return cls(key, value)

    def __str__(self) -> str:
        return "##{key}={value}\n".format(
            key=str(self.key),
            value=str(self.value))
