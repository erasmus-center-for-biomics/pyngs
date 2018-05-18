"""Test the current package."""

import sys
import os.path
import unittest
from tests.test_interval import IntervalTest
from tests.test_gff import GFFTest

if __name__ == "__main__":
    print(sys.version)
    print(sys.executable)
    unittest.main(verbosity=2)
