from opencloning.dna_utils import align_sanger_traces
from pydna.dseqrecord import Dseqrecord

# Linear Sequence Example
linear_record = Dseqrecord('ATGCGTACGTAGCTAGC', circular=False)

traces = [
    'GTACGTAG',  # exact match with middle of reference
    'GCTCGC',  # partial match with end of reference
]

aligned_linear = align_sanger_traces(linear_record, traces)

for i, trace in enumerate(aligned_linear):
    print(f"Sequence {i}: {trace}")

# Circular Sequence Example
long_seq = (
    'ATGATTAGTATATACCTACACAATTACAATATGTAGATATATAGTATGGTAGTAATGTTACATTTTCGAT'
    'TTAGACTAGCCAACCTGCGTATTTACTTGGAATGCTTACATGATTATATTGAATATTACAAATTTTATTA'
    'TGTCCATATTAATCCAGAGTCACGGTGCAATGGAAAAAGGAAAGGGAGAATGATAGAATAAAAGAAGGAA'
    'AGAAAGAGACAATGTAGCGAAGGCTAGAAAGTGATGTGAAAAAAAAAACCACAAAAAAAATAAATAAGAG'
    'ATCAGAGGGTTAAATGAATGCGCTTTTAAGAAATTCAAATATCTGTAAAGGAGAATCCATTCAATTTCGA'
)

circ_record = Dseqrecord(long_seq, circular=True)

trace1 = long_seq[0:50]
trace2 = long_seq[100:150]
# One trace aligns across the origin
trace3 = long_seq[-30:] + long_seq[:30]

sanger_traces = [trace1, trace2, trace3]

aligned_circ = align_sanger_traces(circ_record, sanger_traces)

for i, trace in enumerate(aligned_circ):
    print(f"Sequence {i}: {trace}")
