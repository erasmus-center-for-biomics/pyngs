"""Test the interval class."""

import unittest
import pyngs.interval


class IntervalTest(unittest.TestCase):
    """Test the interval class."""

    def setUp(self):
        """Prepare the test cases."""
        self.interval_a = pyngs.interval.Interval(0, 10, 1000)
        self.interval_b = pyngs.interval.Interval(0, 900, 1100)
        self.interval_c = pyngs.interval.Interval(0, 5, 15)
        self.interval_d = pyngs.interval.Interval(0, 11, 999)
        self.interval_e = pyngs.interval.Interval(0, 1001, 2000)
        
    def test_overlap(self):
        """Test the overlap between intervals."""
        self.assertEqual(self.interval_a.overlaps_with(self.interval_a), True)
        self.assertEqual(self.interval_a.overlaps_with(self.interval_b), True)
        self.assertEqual(self.interval_a.overlaps_with(self.interval_c), True)
        self.assertEqual(self.interval_a.overlaps_with(self.interval_d), True)
        self.assertEqual(self.interval_a.overlaps_with(self.interval_e), False)

    def test_contains(self):
        """Test wether one interval is contained in another interval."""
        self.assertEqual(self.interval_a.contained_in(self.interval_b), False)
        self.assertEqual(self.interval_a.contained_in(self.interval_c), False)
        self.assertEqual(self.interval_a.contained_in(self.interval_a), True)
        self.assertEqual(self.interval_a.contained_in(self.interval_d), False)
        self.assertEqual(self.interval_a.contained_in(self.interval_e), False)

if __name__ == "__main__":
    unittest.main()
