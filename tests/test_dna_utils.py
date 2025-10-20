import os
import unittest
from unittest.mock import patch
from opencloning.dna_utils import sum_is_sticky, align_sanger_traces, permutate_trace
from pydna.dseq import Dseq
from pydna.parsers import parse
from pydna.dseqrecord import Dseqrecord
from Bio.Seq import reverse_complement

test_files = os.path.join(os.path.dirname(__file__), 'test_files')


class PartialStickyTest(unittest.TestCase):
    # General test functions
    def expectTrue(self, seq_left, seq_right, partial):
        with self.subTest():
            self.assertTrue(sum_is_sticky(seq_left.three_prime_end(), seq_right.five_prime_end(), partial))

    def expectFalse(self, seq_left, seq_right, partial):
        with self.subTest():
            self.assertFalse(sum_is_sticky(seq_left.three_prime_end(), seq_right.five_prime_end(), partial))

    # Specific cases
    def test_blunt_ends(self):
        for partial in [False, True]:
            self.expectFalse(Dseq('ACGT'), Dseq('ACGT'), partial)

    def test_sticky_ends_full_overlap_3(self):
        seq1 = Dseq('ACGTAAA', 'ACGT', ovhg=0)
        seq2 = Dseq('ACGT', 'ACGTTTT', ovhg=3)

        for partial in [False, True]:
            self.expectTrue(seq1, seq2, partial)

    def test_sticky_ends_full_overlap_5(self):
        seq1 = Dseq('ACGT', 'TTTACGT', ovhg=0)
        seq2 = Dseq('AAAACGT', 'ACGT', ovhg=-3)

        for partial in [False, True]:
            self.expectTrue(seq1, seq2, partial)

    def test_sticky_ends_partial_overlap_3(self):
        seq1 = Dseq('ACGTAA', 'ACGT', ovhg=0)
        seq2 = Dseq('ACGT', 'ACGTTTT', ovhg=3)

        self.expectTrue(seq1, seq2, True)
        self.expectFalse(seq1, seq2, False)

        seq3 = Dseq('ACGTAAA', 'ACGT', ovhg=0)
        seq4 = Dseq('ACGT', 'ACGTTT', ovhg=2)

        self.expectTrue(seq3, seq4, True)
        self.expectFalse(seq3, seq4, False)

    def test_sticky_ends_partial_overlap_5(self):
        seq1 = Dseq('ACGT', 'TTACGT', ovhg=0)
        seq2 = Dseq('AAAACGT', 'ACGT', ovhg=-3)

        self.expectTrue(seq1, seq2, True)
        self.expectFalse(seq1, seq2, False)

        seq3 = Dseq('ACGT', 'TTTACGT', ovhg=0)
        seq4 = Dseq('AAACGT', 'ACGT', ovhg=-2)

        self.expectTrue(seq3, seq4, True)
        self.expectFalse(seq3, seq4, False)

    def test_sticky_ends_max_len(self):
        # Ensures that all possible overlapping lengths are covered
        seq1 = Dseq('ACGT', 'GTACGT', ovhg=0)
        seq2 = Dseq('ACAACGT', 'ACGT', ovhg=-3)

        self.expectTrue(seq1, seq2, True)
        self.expectFalse(seq1, seq2, False)


class AlignSangerTrackTest(unittest.TestCase):

    def test_align_sanger_traces(self):
        seq = parse(os.path.join(test_files, 'GIN11M86.gb'))[0]
        trace = 'ttgcagcattttgtctttctataaaaatgtgtcgttcctttttttcattttttggcgcgtcgcctcggggtcgtatagaatat'
        alignment = align_sanger_traces(seq, [trace])
        self.assertTrue(alignment[1].startswith('-' * 152 + trace.upper()))

        # Fails with non-nucleotide sequences
        trace_wrong = 'helloworld'
        with self.assertRaises(ValueError):
            align_sanger_traces(seq, [trace_wrong])

        # Works with degenerate sequences
        trace_degenerate = trace.replace('g', 'n')
        alignment_degenerate = align_sanger_traces(seq, [trace_degenerate])
        self.assertTrue(alignment_degenerate[1].startswith('-' * 152 + trace_degenerate.upper()))

        # If the trace aligns to the reverse complement, it returns the trace in the opposite orientation
        trace_rc = reverse_complement(trace)
        alignment_rc = align_sanger_traces(seq, [trace_rc])
        self.assertEqual(alignment_rc[1], alignment[1])

        # Works with circular sequences
        seq = seq.looped()
        end_of_seq = 'ggtggtgggattggtataaagtggtagggtaagtatgtgtgtattatttacgatc'
        start_of_seq = 'gatcaataacagtgtttgtggagca'
        trace_across_origin = end_of_seq + start_of_seq
        alignment = align_sanger_traces(seq, [trace_across_origin])
        # The alignment returns the sequence without shifting it
        self.assertEqual(alignment[0].upper(), str(seq.seq))
        self.assertTrue(
            alignment[1],
            end_of_seq + (len(seq) - len(trace_across_origin)) * '-' + start_of_seq,
        )

        # Works with circular sequences reversed
        trace_across_origin_rc = reverse_complement(trace_across_origin)
        alignment_rc = align_sanger_traces(seq, [trace_across_origin_rc])
        self.assertEqual(alignment_rc[0].upper(), str(seq.seq))
        self.assertEqual(alignment_rc[1], alignment[1])

    def test_align_sanger_traces_multiple(self):
        seq = parse(os.path.join(test_files, 'GIN11M86.gb'))[0].looped()
        trace = 'ttgcagcattttgtctttctataaaaatgtgtcgttcctttttttcattttttggcgcgtcgcctcggggtcgtatagaatatg'
        trace_rc = reverse_complement(trace)

        # Works with multiple traces
        alignments = align_sanger_traces(seq, [trace, trace, trace_rc])
        self.assertEqual(alignments[0].upper().replace('-', ''), str(seq.seq))
        self.assertEqual(len(alignments), 4)
        self.assertEqual(alignments[1], alignments[2])
        # TODO: this has to do with https://github.com/manulera/OpenCloning_frontend/issues/336
        self.assertEqual(alignments[1], alignments[3])

    def test_permutate_trace_error(self):
        with self.assertRaises(RuntimeError):
            permutate_trace('hello', 'world')

        with self.assertRaises(RuntimeError):
            permutate_trace('a', 'aa')

    def test_binaries_missing(self):
        seq = Dseqrecord('ACGT')

        # Test mars missing
        with patch('shutil.which', lambda x: None if x == 'mars' else True):
            with self.assertRaises(RuntimeError) as cm:
                align_sanger_traces(seq, ['ACGT'])
            self.assertEqual(str(cm.exception), "'mars' executable not found in PATH")

        # Test mafft missing
        with patch('shutil.which', lambda x: None if x == 'mafft' else True):
            with self.assertRaises(RuntimeError) as cm:
                align_sanger_traces(seq, ['ACGT'])
            self.assertEqual(str(cm.exception), "'mafft' executable not found in PATH")
