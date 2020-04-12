from typing import List
from .header import Header


class HeaderIndex:

    def __init__(self, headers: List[Header]) -> None:
        """."""
        self.headers = headers
        self.info_index = {}
        self.format_index = {}

        for idx, header in enumerate(self.headers):
            if header.section == "INFO":
                self.info_index[header.id] = idx
            if header.section == "FORMAT":
                self.format_index[header.id] = idx

    def get(self, section: str, id: str) -> Header:
        if section == "INFO":
            return self.headers[self.info_index[id]]
        elif section == "FORMAT":
            return self.headers[self.format_index[id]]
