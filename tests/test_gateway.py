from opencloning.gateway import gateway_overlap, find_gateway_sites
import pydna.assembly2 as assembly
import glob
from opencloning.dna_functions import custom_file_parser
from pydna.dseqrecord import Dseqrecord
from Bio.SeqFeature import SimpleLocation
import os

test_files = os.path.join(os.path.dirname(__file__), 'test_files')


def algoBP(x, y, _):
    return gateway_overlap(x, y, 'BP', greedy=False)


def algoLR(x, y, _):
    return gateway_overlap(x, y, 'LR', greedy=False)


def test_gateway_manual_cloning():

    with open(os.path.join(test_files, 'gateway_manual_cloning/pairing.tsv')) as f:
        for line in f:
            line = line.strip().split('\t')
            if len(line) < 5:
                backbone, pcr_product, entry_vector = line
                backbone_expression = None
                expression = None
            else:
                backbone, pcr_product, entry_vector, backbone_expression, expression = line

            with open(os.path.join(test_files, 'gateway_manual_cloning/' + backbone), 'rb') as f:
                backbone = custom_file_parser(f, 'snapgene')[0]
            with open(os.path.join(test_files, 'gateway_manual_cloning/' + pcr_product), 'rb') as f:
                pcr_product = custom_file_parser(f, 'snapgene')[0]
            with open(os.path.join(test_files, 'gateway_manual_cloning/' + entry_vector), 'rb') as f:
                entry_vector = custom_file_parser(f, 'snapgene')[0]

            # Works with the right reaction
            asm = assembly.Assembly([backbone, pcr_product], algorithm=algoBP)
            out = asm.assemble_circular()
            seguids = [seq.seguid() for seq in out]
            assert entry_vector.seguid() in seguids

            # Does not work with the wrong reaction
            asm = assembly.Assembly([backbone, pcr_product], algorithm=algoLR)
            out = asm.assemble_circular()
            assert len(out) == 0

            if backbone_expression is not None and expression is not None and len(out):
                with open(os.path.join(test_files, 'gateway_manual_cloning/' + backbone_expression), 'rb') as f:
                    backbone_expression = custom_file_parser(f, 'snapgene')[0]
                with open(os.path.join(test_files, 'gateway_manual_cloning/' + expression), 'rb') as f:
                    expression = custom_file_parser(f, 'snapgene')[0]

                # Works with the right reaction
                asm = assembly.Assembly([backbone_expression, entry_vector], algorithm=algoLR)

                out = asm.assemble_circular()
                seguids = [seq.seguid() for seq in out]
                assert expression.seguid() in seguids

                # Does not work with the wrong reaction
                asm = assembly.Assembly([backbone_expression, entry_vector], algorithm=algoBP)

                out = asm.assemble_circular()
                assert len(out) == 0

    example_valerie = glob.glob(os.path.join(test_files, 'gateway_manual_cloning/example_valerie/*.dna'))
    inputs = list()
    for file in example_valerie:
        with open(file, 'rb') as f:
            seq = custom_file_parser(f, 'snapgene')[0]
            seq.name = file.split('/')[-1]
            inputs.append(seq)

    asm = assembly.Assembly(inputs, algorithm=algoLR, use_all_fragments=True, use_fragment_order=False)

    out = asm.assemble_circular()


def test_find_gateway_sites():
    seq = Dseqrecord('CAACTTTGTATACAAAAGTTGaaaCAACTTTGTATAATAAAGTTG')
    assert find_gateway_sites(seq, greedy=False) == {
        'attB5': [SimpleLocation(0, 21, 1)],
        'attB3': [SimpleLocation(24, 45, 1)],
    }

    # Works on circular sequences
    seq = Dseqrecord('CAACTTTGTATACAAAAGTTGaaaCAACTTTGTATAATAAAGTTG', circular=True)
    for shift in range(len(seq)):
        seq_shifted = seq.shifted(shift)
        sites = find_gateway_sites(seq_shifted, greedy=False)
        assert str(sites['attB5'][0].extract(seq_shifted).seq) == 'CAACTTTGTATACAAAAGTTG'
        assert str(sites['attB3'][0].extract(seq_shifted).seq) == 'CAACTTTGTATAATAAAGTTG'

    # Finds several sites of the same type
    seq = Dseqrecord('CAACTTTGTATAATAAAGTTGaaaCAACTTTGTATAATAAAGTTG')
    assert find_gateway_sites(seq, greedy=False) == {'attB3': [SimpleLocation(0, 21, 1), SimpleLocation(24, 45, 1)]}
