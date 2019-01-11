import unittest

import pyngs.bed as bed


def per_line(text):
    """Return the data per line."""
    for line in text.split("\n"):
        yield line


class TestReader(unittest.TestCase):
    """Test the bed parsers."""

    def setUp(self):
        self.data = {
            "bed3": """track id="bed3"
chr1\t213941196\t213942363
chr1\t213942363\t213943530
chr1\t213943530\t213944697
chr2\t158364697\t158365864
chr2\t158365864\t158367031
chr3\t127477031\t127478198
chr3\t127478198\t127479365
chr3\t127479365\t127480532
chr3\t127480532\t127481699""",
            "bedgraph": """chr1\t213941196\t213942363\t1.0
chr1\t213942363\t213943530\t1.0
chr1\t213943530\t213944697\t1.0
chr2\t158364697\t158365864\t0.0
chr2\t158365864\t158367031\t1
chr3\t127477031\t127478198\t3
chr3\t127478198\t127479365\t1.0
chr3\t127479365\t127480532\t1.0
chr3\t127480532\t127481699\t1.0""",
            "bed5": """chr7\t127471196\t127472363\tPos1\t0\t+
chr7\t127472363\t127473530\tPos2\t0\t+
chr7\t127473530\t127474697\tPos3\t0\t+
chr7\t127474697\t127475864\tPos4\t0\t+
chr7\t127475864\t127477031\tNeg1\t0\t-
chr7\t127477031\t127478198\tNeg2\t0\t-
chr7\t127478198\t127479365\tNeg3\t0\t-
chr7\t127479365\t127480532\tPos5\t0\t+
chr7\t127480532\t127481699\tNeg4\t0\t-""",
            "gtf": """1\ttranscribed_unprocessed_pseudogene\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
1\tprocessed_transcript\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_sourc e "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";""",
            "gff": """X\tEnsembl\tRepeat\t2419108\t2419128\t42\t.\t.\thid=trf; hstart=1; hend=21
X\tEnsembl\tRepeat\t2419108\t2419410\t2502\t-\t.\thid=AluSx; hstart=1; hend=303
X\tEnsembl\tRepeat\t2419108\t2419128\t0\t.\t.\thid=dust; hstart=2419108; hend=2419128
X\tEnsembl\tPred.trans.\t2416676\t2418760\t450.19\t-\t2\tgenscan=GENSCAN00000019335
X\tEnsembl\tVariation\t2413425\t2413425\t.\t+\t.\t
X\tEnsembl\tVariation\t2413805\t2413805\t.\t+\t.\t"""
        }

    def test_bed3(self):
        entries = 0
        for entry in bed.reader(per_line(self.data["bed3"]), ftype="bed3"):
            entries += 1
            self.assertEqual(3, len(entry))
        self.assertEqual(entries, 9)

    def test_bedgraph(self):
        entries = 0
        for entry in bed.reader(per_line(self.data["bedgraph"]), ftype="bedgraph"):
            entries += 1
            self.assertEqual(4, len(entry))
        self.assertEqual(entries, 9)

    def test_bed5(self):
        entries = 0
        for entry in bed.reader(per_line(self.data["bed5"]), ftype="bed5"):
            entries += 1
            self.assertEqual(6, len(entry))
        self.assertEqual(entries, 9)

    def test_gtf(self):
        entries = 0
        for entry in bed.reader(per_line(self.data["gtf"]), ftype="gtf"):
            entries += 1
            self.assertEqual(9, len(entry))
        self.assertEqual(entries, 2)

    def test_gff(self):
        entries = 0
        for entry in bed.reader(per_line(self.data["gff"]), ftype="gff"):
            entries += 1
            self.assertEqual(9, len(entry))
        self.assertEqual(entries, 6)


if __name__ == "__main__":
    unittest.main()
