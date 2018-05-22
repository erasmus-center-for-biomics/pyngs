"""Test the Alignment class."""

import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import pyngs.interval


class GFFTest(unittest.TestCase):
    """A class to test the GFF methods."""

    def setUp(self):
        """Prepare the test cases."""
        data = """#!genome-build GRCh38.p5
#!genome-version GRCh38
#!genome-date 2013-12
#!genome-build-accession NCBI:GCA_000001405.20
#!genebuild-last-updated 2015-10
1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_so
1\thavana\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"
1\thavana\texon\t12613\t12721\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"
1\thavana\texon\t13221\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"
1\thavana\ttranscript\t12010\t13670\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000450305"; transcript_version "2"; gene_name "DDX11L1"; gene_so
1\thavana\texon\t12010\t12057\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000450305"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"
1\thavana\texon\t12179\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000450305"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"
1\thavana\texon\t12613\t12697\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000450305"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"
1\thavana\texon\t12975\t13052\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000450305"; transcript_version "2"; exon_number "4"; gene_name "DDX11L1"
1\thavana\texon\t13221\t13374\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000450305"; transcript_version "2"; exon_number "5"; gene_name "DDX11L1"
1\thavana\texon\t13453\t13670\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000450305"; transcript_version "2"; exon_number "6"; gene_name "DDX11L1"
1\thavana\tgene\t14404\t29570\t.\t-\t.\tgene_id "ENSG00000227232"; gene_version "5"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"
"""
        self.handle = StringIO(data)

    def test_gff_methods(self):
        """Test the GFF methods."""

        # Get the first entry
        parser = pyngs.interval.gff_file(self.handle)
        entry = list(next(parser))
        self.assertEqual(entry[0][0], "chromosome")
        self.assertEqual(entry[0][1], "1")
        self.assertEqual(entry[1][0], "source")
        self.assertEqual(entry[1][1], "havana")
        self.assertEqual(entry[2][0], "type")
        self.assertEqual(entry[2][1], "gene")
        self.assertEqual(entry[3][0], "start")
        self.assertEqual(entry[3][1], 11869)
        self.assertEqual(entry[4][0], "end")
        self.assertEqual(entry[4][1], 14409)
        self.assertEqual(entry[5][0], "score")
        self.assertEqual(entry[5][1], ".")
        self.assertEqual(entry[6][0], "strand")
        self.assertEqual(entry[6][1], "+")
        self.assertEqual(entry[7][0], "frame")
        self.assertEqual(entry[7][1], ".")
        self.assertEqual(entry[8][0], "comment")

        # Parse the comment
        commentfields = pyngs.interval.gff_parse_comment(entry[8][1])
        self.assertEqual(len(commentfields), 5)
        self.assertEqual(commentfields[0][0], "gene_id")
        self.assertEqual(commentfields[0][1], "ENSG00000223972")
        self.assertEqual(commentfields[1][0], "gene_version")
        self.assertEqual(commentfields[1][1], "5")
        self.assertEqual(commentfields[2][0], "gene_name")
        self.assertEqual(commentfields[2][1], "DDX11L1")
        self.assertEqual(commentfields[3][0], "gene_source")
        self.assertEqual(commentfields[3][1], "havana")
        self.assertEqual(commentfields[4][0], "gene_biotype")
        self.assertEqual(commentfields[4][1], "transcribed_unprocessed_pseudogene")

        # convert to Interval test
        chromosomes = []
        with self.assertRaises(ValueError):
            pyngs.interval.gff_entry_to_interval(
                entry,
                chromosomes=chromosomes,
                allow_add=False)

        # Test the interval
        interval = pyngs.interval.gff_entry_to_interval(
            entry,
            chromosomes=chromosomes,
            allow_add=True)
        self.assertEqual(interval.chromosome, 0)
        self.assertEqual(interval.start, 11869)
        self.assertEqual(interval.end, 14409)

        # Count the rest of the entries
        count = 0
        for _ in parser:
            count += 1
        self.assertEqual(count, 12)



