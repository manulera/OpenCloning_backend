from itertools import product
from pydna.dseqrecord import Dseqrecord

from .dna_utils import compute_regex_site, dseqrecord_finditer

# This is the original loxP sequence, here for reference
LOXP_SEQUENCE = 'ATAACTTCGTATAGCATACATTATACGAAGTTAT'

# This is a consensus sequence, from this Addgene blog post: https://blog.addgene.org/plasmids-101-cre-lox
# IMPORTANT: Because it is palyndromic, we only look for it in the forward direction, if this was changed
# to a non-palindromic sequence, you would need to look for matches reversing it, like in Gateway cloning
LOXP_CONSENSUS = 'ATAACTTCGTATANNNTANNNTATACGAAGTTAT'


loxP_regex = compute_regex_site(LOXP_CONSENSUS)


def cre_loxP_overlap(x: Dseqrecord, y: Dseqrecord, _l: None = None) -> list[tuple[int, int, int]]:
    """Find matching loxP sites between two sequences."""
    out = list()
    matches_x = dseqrecord_finditer(loxP_regex, x)
    matches_y = dseqrecord_finditer(loxP_regex, y)

    for match_x, match_y in product(matches_x, matches_y):
        value_x = match_x.group()
        value_y = match_y.group()
        if value_x == value_y:
            out.append((match_x.start(), match_y.start(), len(value_x)))
    return out
