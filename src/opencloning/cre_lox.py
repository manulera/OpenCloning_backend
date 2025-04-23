from itertools import product
from pydna.dseqrecord import Dseqrecord
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import reverse_complement
from .dna_utils import compute_regex_site, dseqrecord_finditer

# We create a dictionary to map ambiguous bases to their consensus base
# For example, ambigous_base_dict['ACGT'] -> 'N'
ambiguous_base_dict = {}
for ambiguous, bases in ambiguous_dna_values.items():
    ambiguous_base_dict[''.join(sorted(bases))] = ambiguous

# To handle N values
ambiguous_base_dict['N'] = 'N'

# This is the original loxP sequence, here for reference
LOXP_SEQUENCE = 'ATAACTTCGTATAGCATACATTATACGAAGTTAT'

loxP_sequences = [
    # https://blog.addgene.org/plasmids-101-cre-lox
    # loxP
    'ATAACTTCGTATANNNTANNNTATACGAAGTTAT',
    # PMID:12202778
    # lox66
    'ATAACTTCGTATANNNTANNNTATACGAACGGTA',
    # lox71
    'TACCGTTCGTATANNNTANNNTATACGAAGTTAT',
]

loxP_consensus = ''

for pos in range(len(LOXP_SEQUENCE)):
    all_letters = set(seq[pos] for seq in loxP_sequences)
    key = ''.join(sorted(all_letters))
    loxP_consensus += ambiguous_base_dict[key]

# We compute the regex for the forward and reverse loxP sequences
loxP_regex = (compute_regex_site(loxP_consensus), compute_regex_site(reverse_complement(loxP_consensus)))


def cre_loxP_overlap(x: Dseqrecord, y: Dseqrecord, _l: None = None) -> list[tuple[int, int, int]]:
    """Find matching loxP sites between two sequences."""
    out = list()
    for pattern in loxP_regex:
        matches_x = dseqrecord_finditer(pattern, x)
        matches_y = dseqrecord_finditer(pattern, y)

        for match_x, match_y in product(matches_x, matches_y):
            value_x = match_x.group()
            value_y = match_y.group()
            if value_x[13:21] == value_y[13:21]:
                out.append((match_x.start() + 13, match_y.start() + 13, 8))
    # Unique values (keeping the order)
    unique_out = []
    for item in out:
        if item not in unique_out:
            unique_out.append(item)
    return unique_out
