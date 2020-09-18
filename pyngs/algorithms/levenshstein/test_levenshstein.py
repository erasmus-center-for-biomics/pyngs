import unittest

import pyngs.algorithms.levenshstein as levenshstein


class TestLevenshstein(unittest.TestCase):

    def test_matching(self):
        for case in self.cases:
            dst = levenshstein.Distance(
                len(case["source"]),
                len(case["target"]))
            dst(case["source"], case["target"])
            self.assertEqual(case["expected"], dst.score)

    def setUp(self):
        self.cases = [
            {
                # all matches
                "source": "ACGTA",
                "target": "ACGTA",
                "expected": 0
            },
            {
                # all matches, single N
                "source": "ACGTA",
                "target": "ACNTA",
                "expected": 0
            },
            {
                # all N
                "source": "NNNNNNNNN",
                "target": "NNNNNNNNN",
                "expected": 0
            },
            {
                # single mismatch
                "source": "ACGTA",
                "target": "AAGTA",
                "expected": 1
            },
            {
                # deletion
                "source": "ACGTA",
                "target": "ACTA",
                "expected": 1
            },
            {
                # insertion
                "source": "ACTA",
                "target": "ACGTA",
                "expected": 1
            }
        ]