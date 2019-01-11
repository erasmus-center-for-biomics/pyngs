import unittest

from pyngs.bases import iupac


class TestIUPAC(unittest.TestCase):
    """Test the IUPAC methods."""

    def test_complement(self):
        """Test the complement method."""
        sequence = "ACTGURYSWKMBDHVN"
        expected = "NBDHVKMSWRYACAGT"
        observed = iupac.reverse_complement(sequence)
        self.assertEquals(expected, observed)

    def test_to_regexp(self):
        """Test the regular expression conversion."""
        sequence = "ACTGURYSWKMBDHVN"
        expected = "ACTGU[A,G][C,T][G,C][A,T][G,T][A,C]" \
                   "[C,G,T][A,G,T][A,C,T][A,C,G]."
        observed = iupac.to_regexp(sequence)
        self.assertEquals(expected, observed)

    def test_to_list(self):
        """Test the to list conversion."""
        sequence = "ACTGURYSWKMBDHVN"
        expected = [
            "A",
            "C",
            "T",
            "G",
            "U",
            "AG",
            "CT",
            "GC",
            "AT",
            "GT",
            "AC",
            "CGT",
            "AGT",
            "ACT",
            "ACG",
            ""]
        observed = iupac.to_list(sequence)
        self.assertListEqual(expected, observed)

        expected[-1] = "ACGT"
        observed = iupac.to_list(sequence, fill_n=True)
        self.assertListEqual(expected, observed)


if __name__ == "__main__":
    unittest.main()
