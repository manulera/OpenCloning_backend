from opencloning.cre_lox import cre_loxP_overlap
from pydna.dseqrecord import Dseqrecord
import unittest


loxP_sequence = 'ATAACTTCGTATAGCATACATTATACGAAGTTAT'
loxP_sequence2 = 'ATAACTTCGTATAGTTTACATTATACGAAGTTAT'


class TestCreLox(unittest.TestCase):
    def test_cre_loxP_overlap(self):
        # Works with consensus sequence and other sequence
        for loxP in [loxP_sequence, loxP_sequence2]:
            seqA = Dseqrecord('aa' + loxP + 'acgt')
            seqB = Dseqrecord('ccc' + loxP + 'acgt')
            self.assertEqual(cre_loxP_overlap(seqA, seqB), [(2, 3, 34)])
            self.assertEqual(cre_loxP_overlap(seqA.reverse_complement(), seqB.reverse_complement()), [(4, 4, 34)])
            self.assertEqual(cre_loxP_overlap(seqA, seqB.reverse_complement()), [])
            self.assertEqual(cre_loxP_overlap(seqA.reverse_complement(), seqB), [])

        # It does not mix different loxP sites, even if they match the consensus
        seqA = Dseqrecord('aa' + loxP_sequence + 'acgt')
        seqB = Dseqrecord('ccc' + loxP_sequence2 + 'acgt')
        self.assertEqual(cre_loxP_overlap(seqA, seqB), [])
        self.assertEqual(cre_loxP_overlap(seqA.reverse_complement(), seqB.reverse_complement()), [])
        self.assertEqual(cre_loxP_overlap(seqA, seqB.reverse_complement()), [])
        self.assertEqual(cre_loxP_overlap(seqA.reverse_complement(), seqB), [])
