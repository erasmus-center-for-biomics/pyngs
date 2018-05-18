"""Test the Alignment class."""

import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import pyngs.interval


class AlignmentTest(unittest.TestCase):
    """A class to test the Alignment class."""

    def setUp(self):
        """Prepare the test cases."""
        self.alignment_a = pyngs.alignment.Alignment(
            "HISEQ:351:HCG3TBCXX:2:2113:6295:13667:sC2_r18_c39:AATTAGCGTC",
            16, 1, 24217, 60,
            "86M",
            "*", 0, 0,
            "GGTTGTGTACTGTTTATATGAATATTAGTGATTACAAAATATTATCAATAGATCTTTTCACATCTGTTATTGAAGAACTAATCATC",
            "<0HD1<111111EEF1F1D<11G1GD11<11C<111111G1HD<<111<111D1G<1<1<1111F1HD1HD<111111GC@DB0D<",
            [
                ("AS", "i", "-17"), ("XN", "i", "0"), ("XM", "i", "5"),
                ("XO", "i", "0"), ("XG", "i", "0"), ("NM", "i", "5"),
                ("MD", "Z", "7C2G41C3G5A23"), ("YT", "Z", "UU"),
                ("NH", "i", "1"), ("um", "Z", "AATTAGCGTC"),
                ("wb", "Z", "sC2_r18_c39"), ("RG", "Z", "chip")]
        )
        self.alignment_b = pyngs.alignment.Alignment(
            "HISEQ:351:HCG3TBCXX:2:2113:6295:13667:sC2_r18_c39:AATTAGCGTC",
            0, 1, 24217, 60,
            "86M",
            "*", 0, 0,
            "GGTTGTGTACTGTTTATATGAATATTAGTGATTACAAAATATTATCAATAGATCTTTTCACATCTGTTATTGAAGAACTAATCATC",
            "<0HD1<111111EEF1F1D<11G1GD11<11C<111111G1HD<<111<111D1G<1<1<1111F1HD1HD<111111GC@DB0D<",
            [
                ("AS", "i", "-17"), ("XN", "i", "0"), ("XM", "i", "5"),
                ("XO", "i", "0"), ("XG", "i", "0"), ("NM", "i", "5"),
                ("MD", "Z", "7C2G41C3G5A23"), ("YT", "Z", "UU"),
                ("NH", "i", "1"), ("um", "Z", "AATTAGCGTC"),
                ("wb", "Z", "sC2_r18_c39"), ("RG", "Z", "chip")]
        )

    def test_alignment(self):
        """Test the overlap between intervals."""
        self.assertEqual(self.alignment_a.is_duplicate(), False)
        self.assertEqual(self.alignment_a.is_unmapped(), False)
        self.assertEqual(self.alignment_a.is_reverse(), True)

        self.assertEqual(self.alignment_b.is_duplicate(), False)
        self.assertEqual(self.alignment_b.is_unmapped(), False)
        self.assertEqual(self.alignment_b.is_reverse(), False)


class AlignmentParserTest(unittest.TestCase):
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
        self.handle = StringIO(data)

    def test_parser(self):
        """Parse the SAM file."""
        factory = pyngs.interval.AlignmentFactory()

        entries = 0
        for _ in factory.iterator(self.handle):
            entries += 1

        self.assertEqual(len(factory.header), 6)
        self.assertEqual(len(factory.chromosome_list), 3)
        self.assertEqual(entries, 3)
