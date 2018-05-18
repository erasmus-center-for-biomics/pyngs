"""Test the current package."""

import unittest
import sys
import os

for key in sorted(os.environ.keys()):
    value = os.environ[key]
    sys.stderr.write("{key}\t{value}\n".format(key=key, value=value))

from tests.interval import IntervalTest
from tests.gff import GFFTest
from tests.alignment import AlignmentTest, AlignmentParserTest

if __name__ == "__main__":
    print(sys.version)
    print(sys.executable)
    unittest.main(verbosity=2)
