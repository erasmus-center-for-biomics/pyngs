import unittest

import pyngs.sam as sam
import pyngs.sam.consensus as consensus


class ConsensusTest(unittest.TestCase):
    """A test for the Consensus method."""

    def setUp(self):
        """Prepare alignments for the test."""
        self.alignments = [
            sam.from_string(
                "A 0 1 1003 60 3S3M1I6M1D4M * 0 0 NNNAACTAAAAAATAAA " +
                "IIIIIIIIIIIIIIIII", sep=" "),
            sam.from_string(
                "B 0 1 1005 60 1M1I10M2S * 0 0 CTAAAAAAGTCANN " +
                "IIIIIIIIIIII??", sep=" "),
            sam.from_string(
                "B1 0 1 1005 60 11M * 0 0 CAAAAAAGTAA " +
                "IIIIIIIIIII", sep=" "),
            sam.from_string(
                "C 0 1 1012 60 10M4S * 0 0 GTAAAAAAAANNNN " +
                "IIIIIIIIII????", sep=" ")]

    def test_consensus(self):
        """Test the consensus alignments."""
        consobj = consensus.Consensus()
        consop = consobj.operations(self.alignments)
        for idx, cns in enumerate(consobj):
            print(idx, cns)

if __name__ == "__main__":
    unittest.main()
