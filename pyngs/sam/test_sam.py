import unittest

import pyngs.sam as sam


def per_line(text):
    """Return the data per line."""
    for line in text.split("\n"):
        yield line


class AlignmentTest(unittest.TestCase):
    """Test the alignment class."""
    def setUp(self):
        tokens = [
            "HISEQ:351:HCG3TBCXX:2:1216:13689:15138:sA1_r38_c58:ACGGGGGGCC",
            "0",
            "1",
            "12069",
            "1",
            "100M",
            "*",
            "0",
            "0",
            "AGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCAACTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAACCAGGCAT",
            "<C1?HHG@1EEHH@HHHCGGHEHEGC1E<DHF1EHECCHFCG?G?EHC?@<<HGFEHHCDGFIHIHGGGHGHG?HHFEHHFFD1FF@HHG?HHHGD@?@B",
            "AS:i:-4", "ZS:i:-4", "XN:i:0"]
        self.alignment = sam.from_tokens(tokens)

    def test_flags_paired(self):
        # test the paired
        self.assertEqual(self.alignment.paired, False)
        self.alignment.paired = True
        self.assertEqual(self.alignment.paired, True)
        self.alignment.paired = False
        self.assertEqual(self.alignment.paired, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_proper_pair(self):
        # test the proper_pair
        self.assertEqual(self.alignment.proper_pair, False)
        self.alignment.proper_pair = True
        self.assertEqual(self.alignment.proper_pair, True)
        self.alignment.proper_pair = False
        self.assertEqual(self.alignment.proper_pair, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_unmapped(self):
        # test the unmapped
        self.assertEqual(self.alignment.unmapped, False)
        self.alignment.unmapped = True
        self.assertEqual(self.alignment.unmapped, True)
        self.alignment.unmapped = False
        self.assertEqual(self.alignment.unmapped, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_mate_unmapped(self):
        # test the mate_unmapped
        self.assertEqual(self.alignment.mate_unmapped, False)
        self.alignment.mate_unmapped = True
        self.assertEqual(self.alignment.mate_unmapped, True)
        self.alignment.mate_unmapped = False
        self.assertEqual(self.alignment.mate_unmapped, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_reverse(self):
        # test the reverse
        self.assertEqual(self.alignment.reverse, False)
        self.alignment.reverse = True
        self.assertEqual(self.alignment.reverse, True)
        self.alignment.reverse = False
        self.assertEqual(self.alignment.reverse, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_mate_reverse(self):
        # test the mate_reverse
        self.assertEqual(self.alignment.mate_reverse, False)
        self.alignment.mate_reverse = True
        self.assertEqual(self.alignment.mate_reverse, True)
        self.alignment.mate_reverse = False
        self.assertEqual(self.alignment.mate_reverse, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_first_in_pair(self):
        # test the first_in_pair
        self.assertEqual(self.alignment.first_in_pair, False)
        self.alignment.first_in_pair = True
        self.assertEqual(self.alignment.first_in_pair, True)
        self.alignment.first_in_pair = False
        self.assertEqual(self.alignment.first_in_pair, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_last_in_pair(self):
        # test the last_in_pair
        self.assertEqual(self.alignment.last_in_pair, False)
        self.alignment.last_in_pair = True
        self.assertEqual(self.alignment.last_in_pair, True)
        self.alignment.last_in_pair = False
        self.assertEqual(self.alignment.last_in_pair, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_secondary_alignment(self):
        # test the secondary_alignment
        self.assertEqual(self.alignment.secondary_alignment, False)
        self.alignment.secondary_alignment = True
        self.assertEqual(self.alignment.secondary_alignment, True)
        self.alignment.secondary_alignment = False
        self.assertEqual(self.alignment.secondary_alignment, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_does_not_pass_filters(self):
        # test the does_not_pass_filters
        self.assertEqual(self.alignment.does_not_pass_filters, False)
        self.alignment.does_not_pass_filters = True
        self.assertEqual(self.alignment.does_not_pass_filters, True)
        self.alignment.does_not_pass_filters = False
        self.assertEqual(self.alignment.does_not_pass_filters, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_duplicate(self):
        # test the duplicate
        self.assertEqual(self.alignment.duplicate, False)
        self.alignment.duplicate = True
        self.assertEqual(self.alignment.duplicate, True)
        self.alignment.duplicate = False
        self.assertEqual(self.alignment.duplicate, False)
        self.assertEqual(self.alignment.flag, 0)

    def test_flags_supplementary_alignment(self):
        # test the supplementary_alignment
        self.assertEqual(self.alignment.supplementary_alignment, False)
        self.alignment.supplementary_alignment = True
        self.assertEqual(self.alignment.supplementary_alignment, True)
        self.alignment.supplementary_alignment = False
        self.assertEqual(self.alignment.supplementary_alignment, False)
        self.assertEqual(self.alignment.flag, 0)


class ReaderTest(unittest.TestCase):
    """Test alignment reading."""

    def setUp(self):
        """Prepare the test."""
        data = """@HD\tVN:1.0\tSO:coordinate
@SQ\tSN:1\tLN:248956422
@SQ\tSN:10\tLN:133797422
@SQ\tSN:11\tLN:135086622
@RG\tID:chip\tCN:ECB\tLB:Wafergen\tSM:chip\tPL:ILLUMINA\tDS:/research_t/Analysis/Wafergen/tests/20161216_test_hiseq/per_sample/chip_73952_hiseq_sA1.srt.bam
@PG\tID:hisat2\tPN:hisat2\tVN:2.0.4\tCL:"/usr/local/hisat2-2.0.4/hisat2-align-s --wrapper basic-0 -p 5 -x /research_t/Reference/indexes/GRCh38/hisat2/grch38_tran/genome_tran -U chip_73952_hiseq_sA1.tr.fastq"
HISEQ:351:HCG3TBCXX:2:1216:13689:15138:sA1_r38_c58:ACGGGGGGCC\t16\t1\t12069\t1\t100M\t*\t0\t0\tAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCAACTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAACCAGGCAT\t<C1?HHG@1EEHH@HHHCGGHEHEGC1E<DHF1EHECCHFCG?G?EHC?@<<HGFEHHCDGFIHIHGGGHGHG?HHFEHHFFD1FF@HHG?HHHGD@?@B\tAS:i:-4\tZS:i:-4\tXN:i:0\tXM:i:1\tXO:i:0\tXG:i:0\tNM:i:1\tMD:Z:51T48\tYT:Z:UU\tNH:i:4\tum:Z:ACGGGGGGCC\twb:Z:sA1_r38_c58\tRG:Z:chip
HISEQ:352:HCFVCBCXX:1:1203:18109:34114:sA1_r38_c58:ACGGGGGGCC\t272\t1\t12069\t1\t100M\t*\t0\t0\tAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAACCAGGCAT\t1IIIIIHHGEFHHFIHGIIHHHHHHEECIIHGHIHEF@EHGHHHE@CC@GHD??HHGIIHHEEHHHFHHHGIHHHHEGIHHHEHHHHIIIIIIIHDCDCD\tAS:i:0\tZS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:100\tYT:Z:UU\tNH:i:4\tum:Z:ACGGGGGGCC\twb:Z:sA1_r38_c58\tRG:Z:chip
HISEQ:352:HCFVCBCXX:2:2115:5066:18397:sA1_r38_c58:ACGGGGGGCC\t272\t1\t12069\t1\t100M\t*\t0\t0\tAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAACCAGGCAT\t1IHGHIIIHGIIIHHEHHCHIIHHGHHECGHHIHHIHHHEIHHG?EHIHIIGIIIIHHIIHHHHHIHHHHIIIGIHIHHIHF?EHHHHIIIIHIIABBBD\tAS:i:0\tZS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:100\tYT:Z:UU\tNH:i:4\tum:Z:ACGGGGGGCC\twb:Z:sA1_r38_c58\tRG:Z:chip
"""
        self.stream = per_line(data)

    def test_parser(self):
        """Parse the SAM file."""
        reader = sam.Reader(self.stream)

        entries = 0
        for _ in reader:
            entries += 1
        self.assertEqual(len(reader.header), 6)
        self.assertEqual(entries, 3)


if __name__ == "__main__":
    unittest.main()
