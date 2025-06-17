from pydna.dseqrecord import Dseqrecord
import opencloning.assembly2 as assembly


def test_insertion_edge_case():
    a = Dseqrecord('cccgaggggaatcgaa')
    b = Dseqrecord('Acccgagggggaatc')

    asm = assembly.Assembly([a, b], limit=5, use_all_fragments=True, use_fragment_order=False)

    product_seqs = set(str(prod.seq) for prod in asm.assemble_insertion())
    expected_seqs = {'cccgagggggaatcgaa', 'Acccgaggggaatc'}
    assert product_seqs == expected_seqs

    b_circ = b.looped()

    for shift in range(len(b_circ)):
        b_shifted = b_circ.shifted(shift)
        asm = assembly.Assembly([a, b_shifted], limit=5, use_all_fragments=True, use_fragment_order=False)

        product_seqs = [str(prod.seq) for prod in asm.assemble_insertion()]
        assert len(product_seqs) == 1
        assert product_seqs[0] == 'cccgagggggaatcgaa'
