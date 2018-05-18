"""Test the current package."""

import sys
import unittest

from tests.interval import IntervalTest
from tests.gff import GFFTest
from tests.alignment import AlignmentTest, AlignmentParserTest

if __name__ == "__main__":
    print(sys.version)
    print(sys.executable)
    unittest.main(verbosity=2)
