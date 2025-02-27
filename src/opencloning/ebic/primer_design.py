from pydna.dseqrecord import Dseqrecord
from Bio.SeqFeature import SimpleLocation
from primer3 import bindings

from ..pydantic_models import PrimerModel
from .primer_design_settings import amanda_settings

padding = 1000

adapter_left_fwd = 'ataGGTCTCtGGAG'
adapter_left_rvs = 'ataGGTCTCtCATT'
adapter_right_fwd = 'ataGGTCTCtGCTT'
adapter_right_rvs = 'ataGGTCTCtAGCG'


def call_primer3(seq: str, seq_args: dict):
    report = bindings.design_primers(
        seq_args={
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': seq,
            **seq_args,
        },
        global_args=amanda_settings,
    )
    return report


def ebic_primers(
    input_seq: Dseqrecord,
    location: SimpleLocation,
    max_inside: int,
    max_outside: int,
) -> tuple[PrimerModel, PrimerModel, PrimerModel, PrimerModel]:
    """Design primers for EBIC"""

    # First, we keep only the part within the padding
    edge_left = location.start - padding
    edge_right = location.end + padding
    if edge_left < 0 or edge_right > len(input_seq):
        raise ValueError('The template is too short for the padding.')

    template_seq = str(input_seq.seq[edge_left:edge_right])
    inside_edge = padding + max_inside
    outside_edge = padding - max_outside

    left_template = template_seq[:inside_edge]
    right_template = template_seq[-inside_edge:]

    seq_args_left = {
        'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': f'0,{int(padding/2)},{outside_edge},{max_outside + max_inside}',
    }
    seq_args_right = {
        'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': f'0,{max_outside + max_inside},{len(right_template) - int(padding/2)},{int(padding/2)}',
    }

    report_left = call_primer3(left_template, seq_args_left)
    report_right = call_primer3(right_template, seq_args_right)

    primer_names = ['left_fwd', 'left_rvs', 'right_fwd', 'right_rvs']
    primer_seqs = [
        adapter_left_fwd + report_left['PRIMER_LEFT'][0]['SEQUENCE'],
        adapter_left_rvs + report_left['PRIMER_RIGHT'][0]['SEQUENCE'],
        adapter_right_fwd + report_right['PRIMER_LEFT'][0]['SEQUENCE'],
        adapter_right_rvs + report_right['PRIMER_RIGHT'][0]['SEQUENCE'],
    ]
    return [PrimerModel(id=0, name=primer_names[i], sequence=primer_seqs[i]) for i in range(4)]
