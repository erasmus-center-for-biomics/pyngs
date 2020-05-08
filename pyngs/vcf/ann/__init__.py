
from typing import Generator, List
from ..utils import quote_tokenizer
from .annotation import Annotation


def parse_ann(text:str) -> Generator[Annotation, None, None]:
    """Parse annotations from an ANN field."""
    for txt in quote_tokenizer(text, sep=","):
        yield Annotation.from_str(txt)

