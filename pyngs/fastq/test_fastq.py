import unittest

import pyngs.fastq as fastq


def per_line(text):
    """Return the data per line."""
    for line in text.split("\n"):
        yield line


class TestFastqGenerator(unittest.TestCase):
    """Test whether the methods in fastq.py behave as expected."""

    def setUp(self):
        """Add a default FastQ to test with."""
        self.data = """@read 1
ACTGACTG
+
!!!!!!!!
@read 2
TGACAAAAAAA
+
IIIIIIIIIII
"""

    def test_fastq(self):
        """Test the fastq generator."""
        gen = per_line(self.data)

        # explicitly create the fastq generator
        # so that we can iterate through it slowly
        fq = fastq.fastq(gen)
        name, seq, qual = next(fq)
        self.assertEquals("read 1", name)
        self.assertEquals("ACTGACTG", seq)
        self.assertEquals("!!!!!!!!", qual)

        name, seq, qual = next(fq)
        self.assertEquals("read 2", name)
        self.assertEquals("TGACAAAAAAA", seq)
        self.assertEquals("IIIIIIIIIII", qual)

    def test_format(self):
        """Test the format method."""
        result = fastq.format("name", "AAAAAA", "######")
        self.assertEquals("@name\nAAAAAA\n+\n######\n", result)

    def test_clean_readname(self):
        """Test the clean_readname method."""
        totest = "part1 part2"
        result = fastq.clean_readname(totest)
        self.assertEquals("part1", result)

    def test_encode_data(self):
        """Test the encode_in_readname method."""
        totest = ["apple", "pear", 1]
        result = fastq.encode_in_readname("pineapple:berry", totest)
        self.assertEquals("pineapple:berry:apple:pear:1", result)


if __name__ == "__main__":
    unittest.main()
