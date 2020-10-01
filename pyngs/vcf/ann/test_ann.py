

import unittest
import pyngs.vcf.ann as ann


class AnnTest(unittest.TestCase):

    def test_ann(self):
        annflds = [a for a in ann.parse_ann(self.ann_field)]

        # make sure we have 2 fields
        self.assertEquals(len(annflds), 2)
        self.assertEquals(annflds[0].allele, "T")
        self.assertEquals(annflds[0].annotation, ["missense_variant"])
        self.assertEquals(annflds[0].impact, "MODERATE")
        self.assertEquals(
            repr(annflds[0]),
            """T|missense_variant|MODERATE|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.1406G>A|p.Gly469Glu|1666/2034|1406/1674|469/557||""")
        self.assertEquals(
            repr(annflds[1]),
            """T|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|transcript|ENST00000430910|processed_pseudogene||n.*397G>A|||||3944|""")

    def setUp(self):
        self.ann_field = """T|missense_variant|MODERATE|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.1406G>A|p.Gly469Glu|1666/2034|1406/1674|469/557||,T|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|transcript|ENST00000430910|processed_pseudogene||n.*397G>A|||||3944|"""

if __name__ == "__main__":
    unittest.main()

