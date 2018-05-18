"""Test the current package."""

import sys
import unittest

from .interval import IntervalTest
from .gff import GFFTest
from .alignment import AlignmentTest, AlignmentParserTest

if __name__ == "__main__":
    print(sys.version)
    print(sys.executable)
    unittest.main(verbosity=2)
