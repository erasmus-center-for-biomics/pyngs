import unittest
import pyngs.vcf as vcf
from pyngs.vcf.utils import quote_tokenizer


def per_line(text):
    """Return the data per line."""
    for line in text.split("\n"):
        yield line


class TestReader(unittest.TestCase):

    def test_tokenizer(self):
        text = 'A;B;"V D";x="XYZ"'
        tokens = [q for q in quote_tokenizer(text, ";")]
        self.assertEqual(len(tokens), 4)
        self.assertEqual(tokens[0], "A")
        self.assertEqual(tokens[1], "B")
        self.assertEqual(tokens[2], '"V D"')
        self.assertEqual(tokens[3], 'x="XYZ"')

    def test_reader(self):
        stream = per_line(self.vcf_data)
        reader = vcf.Reader(stream)
        # test the header length
        self.assertEqual(len(reader.header), 42)

        # test the header format index
        self.assertEqual(len(reader.format.keys()), 10)

        # test the header filter index
        self.assertEqual(len(reader.filter.keys()), 2)

        # test the header info index
        self.assertEqual(len(reader.info.keys()), 20)

        # test the samples
        self.assertEqual(len(reader.samples), 1)

        # test the line parsing
        entries = []
        for entry in reader:
            entries.append(entry)
        self.assertEqual(len(entries), 3)

        vardata = []
        for entry in entries:
            sdata = []
            for sample in entry.samples:
                sentry = []
                for idx, code in enumerate(entry.format):
                    parser = reader.field("FORMAT", code)
                    values = list(parser.interpret(sample[idx]))
                    sentry.append((code, values))
                sdata.append(sentry)
            vardata.append(sdata)

        self.assertEqual(len(vardata), 3)
        self.assertEqual(len(vardata[0]), 1)
        # GT:AD:DP:GQ:PGT:PID:PL
        #  1  2  3  4   5   6  7
        # 0/0:14,0:14:17:.:.:0,17,357
        # print(vardata[0][0])
        self.assertEqual(len(vardata[0][0]), 7)
        self.assertEqual(vardata[0][0][0], ("GT", ["0/0"]))
        self.assertEqual(vardata[0][0][1], ("AD", [14, 0]))
        self.assertEqual(vardata[0][0][2], ("DP", [14]))
        self.assertEqual(vardata[0][0][3], ("GQ", [17]))
        self.assertEqual(vardata[0][0][4], ("PGT", ["."]))
        self.assertEqual(vardata[0][0][5], ("PID", ["."]))
        self.assertEqual(vardata[0][0][6], ("PL", [0, 17, 357]))

        # GT:AD:DP:GQ:PL
        #  1  2  3  4  5
        # 0/0:36,0:36:39:0,39,1052
        self.assertEqual(len(vardata[1][0]), 5)
        self.assertEqual(vardata[1][0][0], ("GT", ["0/0"]))
        self.assertEqual(vardata[1][0][1], ("AD", [36, 0]))
        self.assertEqual(vardata[1][0][2], ("DP", [36]))
        self.assertEqual(vardata[1][0][3], ("GQ", [39]))
        self.assertEqual(vardata[1][0][4], ("PL", [0, 39, 1052]))

        # GT:AD:DP:GQ:PL
        #  1  2  3  4  5
        # 1/1:0,2:2:6:51,6,0
        self.assertEqual(len(vardata[2][0]), 5)
        self.assertEqual(vardata[2][0][0], ("GT", ["1/1"]))
        self.assertEqual(vardata[2][0][1], ("AD", [0, 2]))
        self.assertEqual(vardata[2][0][2], ("DP", [2]))
        self.assertEqual(vardata[2][0][3], ("GQ", [6]))
        self.assertEqual(vardata[2][0][4], ("PL", [51, 6, 0]))

    def setUp(self):
        self.vcf_data = """##fileformat=VCFv4.2
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group"
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
##GATKCommandLine=<ID=GenomicsDBImport,CommandLine="GenomicsDBImport  --genomicsdb-workspace-path gatk/genomicsdb/chrM:1-16571 --variant gatk/haplotypecaller/I18-1012-13.mrg.gvcf.gz --variant
##GATKCommandLine=<ID=GenotypeGVCFs,CommandLine="GenotypeGVCFs  --output gatk/genotypegvcfs/chrM:1-16571.vcf.gz --use-new-qual-calculator true --annotation-group StandardAnnotation --dbsnp /da
##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller  --emit-ref-confidence GVCF --max-alternate-alleles 3 --contamination-fraction-to-filter 0.0 --output gatk/haplotypecaller/I1
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as l
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=RAW_MQ,Number=1,Type=Float,Description="Raw data for RMS Mapping Quality">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
##contig=<ID=chrM,length=16571>
##contig=<ID=chr1,length=249250621>
##source=GenomicsDBImport
##source=GenotypeGVCFs
##source=HaplotypeCaller
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmy-sample
chr1\t10474\t.\tT\tA\t35.77\t.\tAC=2;AF=0.071;AN=28;DP=49;ExcessHet=0.0812;FS=0.000;InbreedingCoeff=0.1803;MLEAC=2;MLEAF=0.071;MQ=22.52;QD=30.62;SOR=1.609\tGT:AD:DP:GQ:PGT:PID:PL\t0/0:14,0:14:17:.:.:0,17,357
chr1\t13116\trs201725126\tT\tG\t2080.88\t.\tAC=16;AF=0.333;AN=48;BaseQRankSum=1.35;ClippingRankSum=0.00;DB;DP=557;ExcessHet=6.8926;FS=33.636;InbreedingCoeff=-0.2152;MLEAC=2;MLEAF=0.019;MQ=22.83;MQRankSum=-1.036e+00;QD=4.72;ReadPosRankSum=1.64;SOR=1.609\tGT:AD:DP:GQ:PL\t0/0:36,0:36:39:0,39,1052
chr1\t16298\trs200451305\tC\tT\t99.86\t.\tAC=4;AF=0.667;AN=6;DB;DP=8;ExcessHet=0.4576;FS=0.000;MLEAC=10;MLEAF=1.00;MQ=13.75;QD=30.02;SOR=2.833\tGT:AD:DP:GQ:PL\t1/1:0,2:2:6:51,6,0
"""
