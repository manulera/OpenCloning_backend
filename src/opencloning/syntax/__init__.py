from pydna.dseqrecord import Dseqrecord
import networkx as nx
from .models import Syntax
from Bio.Seq import reverse_complement


def assign_plasmid_to_syntax_part(plasmid: Dseqrecord, syntax: Syntax) -> str:
    graph = syntax.to_edges_graph()
    assembly_enzyme = syntax.get_assembly_enzyme()
    result = []
    for fragment in plasmid.cut(assembly_enzyme):
        for rc in [True, False]:
            query = fragment.reverse_complement() if rc else fragment
            three_type, three_ovhg = query.seq.three_prime_end()
            five_type, five_ovhg = query.seq.five_prime_end()
            # It must only have 5' overhangs
            if three_type != five_type or five_type != "5'":
                continue
            # It must not contain the recognition site of the enzyme inside
            # since they are always in the backbone, not the part.
            # We use compsite, because the simple search method requires the
            # cutsite to be there, and not sure how behaviour will be querying
            # the overhangs.
            if assembly_enzyme.compsite.search(str(query.seq)) is not None:
                continue

            left_node = three_ovhg.upper()
            right_node = reverse_complement(five_ovhg).upper()

            if left_node in graph and right_node in graph and nx.has_path(graph, left_node, right_node):
                result.append(f"{left_node}-{right_node}")
    return result
