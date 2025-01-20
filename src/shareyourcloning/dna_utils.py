"""
Utility functions moved here to avoid circular imports.
"""

from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
import tempfile
import subprocess
import os
import shutil
from pydna.parsers import parse
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
aligner.open_gap_score = -1.0
aligner.extend_gap_score = 0.0
aligner.mode = 'local'


def sum_is_sticky(three_prime_end: tuple[str, str], five_prime_end: tuple[str, str], partial: bool = False) -> int:
    """Return the overlap length if the 3' end of seq1 and 5' end of seq2 ends are sticky and compatible for ligation.
    Return 0 if they are not compatible."""
    type_seq1, sticky_seq1 = three_prime_end
    type_seq2, sticky_seq2 = five_prime_end

    if 'blunt' != type_seq2 and type_seq2 == type_seq1 and str(sticky_seq2) == str(reverse_complement(sticky_seq1)):
        return len(sticky_seq1)

    if not partial:
        return 0

    if type_seq1 != type_seq2 or type_seq2 == 'blunt':
        return 0
    elif type_seq2 == "5'":
        sticky_seq1 = str(reverse_complement(sticky_seq1))
    elif type_seq2 == "3'":
        sticky_seq2 = str(reverse_complement(sticky_seq2))

    ovhg_len = min(len(sticky_seq1), len(sticky_seq2))
    # [::-1] to try the longest overhangs first
    for i in range(1, ovhg_len + 1)[::-1]:
        if sticky_seq1[-i:] == sticky_seq2[:i]:
            return i
    else:
        return 0


def get_alignment_shift(alignment: Dseq, shift: int) -> int:
    """Shift the alignment by the given number of positions, ignoring gap characters (-).

    Parameters
    ----------
    alignment : Dseq
        The alignment sequence that may contain gap characters (-)
    shift : int
        Number of positions to shift the sequence by

    """

    nucleotides_shifted = 0
    positions_shifted = 0
    corrected_shift = shift if shift >= 0 else len(alignment) + shift
    alignment_str = str(alignment)

    while nucleotides_shifted != corrected_shift:
        if alignment_str[positions_shifted] != '-':
            nucleotides_shifted += 1
        positions_shifted += 1

    return positions_shifted


def align_with_mafft(inputs: list[str], orientation_known: bool) -> list[str]:
    """Align a sanger track to a dseqr sequence"""

    with tempfile.TemporaryDirectory() as tmpdir:
        input_file = os.path.join(tmpdir, 'input.fa')
        with open(input_file, 'w') as f:
            for i, input_seq in enumerate(inputs):
                f.write(f">trace-{i+1}\n{input_seq}\n")

        result = subprocess.run(
            ['mafft', '--nuc'] + (['--adjustdirection'] if orientation_known else []) + [input_file],
            capture_output=True,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(f'MAFFT alignment failed:\n{result.stderr}')

    return [str(s.seq) for s in parse(result.stdout)]


def permutate_trace(reference: str, sanger_trace: str) -> str:
    """Permutate a trace with respect to the reference using MARS"""
    # As an input for MARS, we need the reference + all traces
    # We include traces in both directions, since MARS does not handle
    # reverse complements - see https://github.com/lorrainea/MARS/issues/17#issuecomment-2598314356
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = os.path.join(tmpdir, 'input.fa')
        with open(input_path, 'w') as f:
            f.write(f">ref\n{reference}\n")
            f.write(f">trace\n{sanger_trace}\n")

        output_path = os.path.join(tmpdir, 'output.fa')
        result = subprocess.run(['mars', '-a', 'DNA', '-m', '0', '-i', input_path, '-o', output_path, '-q', '5', '-l', '20', '-P', '1'], capture_output=True, text=True)  # fmt: skip

        if result.returncode != 0:
            raise RuntimeError(f'MARS failed:\n{result.stderr}')

        # read permutated trace
        return str(parse(output_path, 'fasta')[1].seq)


def align_sanger_traces(dseqr: Dseqrecord, sanger_traces: list[str]) -> list[str]:
    """Align a sanger track to a dseqr sequence"""
    query_str = str(dseqr.seq)
    # Check that required executables exist in PATH
    if not shutil.which('mars'):
        raise RuntimeError("'mars' executable not found in PATH")
    if not shutil.which('mafft'):
        raise RuntimeError("'mafft' executable not found in PATH")

    # If the sequence is circular, use MARS to permutate the traces
    if dseqr.circular:
        permutated_traces = []
        for trace in sanger_traces:
            permutated_traces.append(permutate_trace(query_str, trace))
            permutated_traces.append(permutate_trace(query_str, reverse_complement(trace)))

        traces_oriented = []
        # Pairwise-align and keep the best alignment, to decide which orientation to keep
        for fwd, rvs in zip(permutated_traces[::2], permutated_traces[1::2]):
            fwd_score = aligner.score(*align_with_mafft([query_str, fwd], True))
            rvs_score = aligner.score(*align_with_mafft([query_str, rvs], True))

            traces_oriented.append(fwd if fwd_score > rvs_score else rvs)
        sanger_traces = traces_oriented

    return align_with_mafft([query_str, *sanger_traces], True)
