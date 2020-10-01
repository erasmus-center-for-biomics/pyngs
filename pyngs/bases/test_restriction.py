import unittest

from pyngs.bases import restriction


class TestRestriction(unittest.TestCase):

    def setUp(self):
        self.enzymes = {
            "AluI": ("AG/|CT", "AG/|CT"),
            "BccI": ("CCATCNNNN/N|", "/N|NNNNGATGG"),
            "SfcI": ("C/TRYA|G", "C/TRYA|G"),
            "LpnPI": ("C^CDGNNNNNNNNNN/NNNN|", "/NNNN|NNNNNNNNNNCH^GG")
        }

    def test_restriction_sites(self):
        """
        Test inverting a restriction enzyme.

        LpnPI
        0  1 2 3 4 5        10       14      18
        C ^C D G N N N N N N N N N N /N N N N |
        C  C D G N N N N N N N N N N  N N N N ""
           ^                          /       |

        """
        obs = self.enzymes["LpnPI"]
        for part in restriction.partition(obs[0]):
            print(part)

        #
        enzyme = restriction.Enzyme("LpnPI", obs[0])
        print()
        print("Ori")
        print(enzyme)
        print(enzyme.to_sequence())

        rcenzyme = enzyme.reverse_complement()
        print()
        print("RC")
        print(rcenzyme)
        print(rcenzyme.to_sequence())

        # print("recognition", rcenzyme.recognition)
        # print("marks", rcenzyme.marks)
        # print("lead", rcenzyme.lead)
        # print("lag", rcenzyme.lag)
        # print(rcenzyme.to_string())

        # benzyme = restriction.Enzyme(obs[1])
        # print()
        # print("RC(b)")
        # print("recognition", benzyme.recognition)
        # print("marks", benzyme.marks)
        # print("lead", benzyme.lead)
        # print("lag", benzyme.lag)
        # print(benzyme.to_string())

        # bases, modifiers = restriction.to_recognition(obs[0])

        # self.assertEqual(bases[0], "C")
        # self.assertEqual(bases[1], "C")
        # self.assertEqual(bases[2], "D")
        # # self.assertEqual(bases[2], "AGT")
        # self.assertEqual(bases[3], "G")
        # self.assertEqual(bases[4], "N")
        # self.assertEqual(bases[5], "N")
        # self.assertEqual(len(bases), 19)

        # self.assertEqual(modifiers[0], [])
        # self.assertEqual(modifiers[1], ["^"])
        # self.assertEqual(modifiers[2], [])
        # self.assertEqual(modifiers[3], [])
        # self.assertEqual(modifiers[4], [])
        # self.assertEqual(modifiers[5], [])
        # self.assertEqual(modifiers[14], ["/"])
        # self.assertEqual(modifiers[18], ["|"])
        # self.assertEqual(len(modifiers), 19)

        # # get up to the cut sites
        # ldcuts = [c for c in restriction.index(modifiers, "/")]
        # lgcuts = [c for c in restriction.index(modifiers, "|")]
        # msite = [c for c in restriction.index(modifiers, "^")]

        # self.assertEqual(msite[0], 1)
        # self.assertEqual(ldcuts[0], 14)
        # self.assertEqual(lgcuts[0], 18)

        # up_to_lead = "".join(bases[:ldcuts[0]])
        # self.assertEqual(up_to_lead, "CCDGNNNNNNNNNN")

        # up_to_lag = "".join(bases[:lgcuts[0]])
        # self.assertEqual(up_to_lag, "CCDGNNNNNNNNNNNNNN")

        # # check the reverse complement
        # rcbases, rcmods = restriction.reverse_complement(bases, modifiers)
        # rcsequence = restriction.to_sequence(rcbases, rcmods)
        # self.assertEquals(rcsequence, obs[1])

        # # print the reverse complement
        # print("Reverse complement")
        # for idx, base in enumerate(rcbases):
        #     print(idx, base, rcmods[idx])

        # print("Forward")
        # for idx, base in enumerate(bases):
        #     print(idx, base, modifiers[idx])

        # self.assertEqual(restriction.to_sequence(bases, modifiers), obs[0])


if __name__ == "__main__":
    unittest.main()
