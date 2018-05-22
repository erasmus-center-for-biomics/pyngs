"""Test the current package."""

import unittest
import sys
import os
# sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../"))
from tests.interval import IntervalTest
from tests.gff import GFFTest
from tests.alignment import AlignmentTest, AlignmentParserTest
from tests.consensus import ConsensusTest

if __name__ == "__main__":
    print(sys.version)
    print(sys.executable)
    unittest.main(verbosity=2)
