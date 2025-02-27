import asyncio
import os
from primer3 import bindings
import json
from fastapi import UploadFile, Response

from ...pydantic_models import (
    GenomeCoordinatesSource,
    TextFileSequence,
    PrimerModel,
    PCRSource,
    RestrictionAndLigationSource,
    BaseCloningStrategy,
    HomologousRecombinationSource,
)
from .primer_design_settings import amanda_settings
from ...endpoints.external_import import genome_coordinates, read_from_file
from ...endpoints.assembly import pcr, restriction_and_ligation, homologous_recombination
from ...dna_functions import read_dsrecord_from_json

# Settings for design
padding = 1000
max_inside = 20
max_outside = 70
outside_edge = padding - max_outside
inside_edge = padding + max_inside

adapter_left_fwd = 'ataGGTCTCtGGAG'
adapter_left_rvs = 'ataGGTCTCtCATT'
adapter_right_fwd = 'ataGGTCTCtGCTT'
adapter_right_rvs = 'ataGGTCTCtAGCG'


def design_primers_ebic(seq: str, seq_args: dict):

    report = bindings.design_primers(
        seq_args={
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': seq,
            **seq_args,
        },
        global_args=amanda_settings,
    )

    return report


async def main():
    cloning_strategy = BaseCloningStrategy(
        sequences=[],
        sources=[],
        primers=[],
        description='EBIC cloning strategy',
    )

    initial_source = GenomeCoordinatesSource(
        id=0,
        assembly_accession='GCF_002847425.1',
        sequence_accession='NZ_CP025534.1',
        start=2700777 - padding,
        end=2701295 + padding,
        strand=1,
    )

    if os.path.exists('cached_response.json'):
        with open('cached_response.json', 'r') as f:
            locus = json.load(f)
            locus['sequences'][0] = TextFileSequence.model_validate(locus['sequences'][0])
            locus['sources'][0] = GenomeCoordinatesSource.model_validate(locus['sources'][0])
    else:
        print('No cached response found, fetching from server')
        locus = await genome_coordinates(initial_source)
        with open('cached_response.json', 'w') as f:
            json.dump(
                {
                    'sequences': [locus['sequences'][0].model_dump()],
                    'sources': [locus['sources'][0].model_dump()],
                },
                f,
            )

    locus_seq: TextFileSequence = locus['sequences'][0]
    locus_source: GenomeCoordinatesSource = locus['sources'][0]
    locus_dseqr = read_dsrecord_from_json(locus_seq)

    cloning_strategy.add_source_and_sequence(locus_source, locus_seq)

    left_template = locus_dseqr[:inside_edge]
    right_template = locus_dseqr[-inside_edge:]

    seq_args_left = {
        'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': f'0,{int(padding/2)},{outside_edge},{max_outside + max_inside}',
    }
    seq_args_right = {
        'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': f'0,{max_outside + max_inside},{len(right_template) - int(padding/2)},{int(padding/2)}',
    }

    report_left = design_primers_ebic(str(left_template.seq), seq_args_left)
    report_right = design_primers_ebic(str(right_template.seq), seq_args_right)

    primer_names = ['left_fwd', 'left_rvs', 'right_fwd', 'right_rvs']
    primer_seqs = [
        adapter_left_fwd + report_left['PRIMER_LEFT'][0]['SEQUENCE'],
        adapter_left_rvs + report_left['PRIMER_RIGHT'][0]['SEQUENCE'],
        adapter_right_fwd + report_right['PRIMER_LEFT'][0]['SEQUENCE'],
        adapter_right_rvs + report_right['PRIMER_RIGHT'][0]['SEQUENCE'],
    ]
    primers = [PrimerModel(id=0, name=primer_names[i], sequence=primer_seqs[i]) for i in range(4)]
    for primer in primers:
        cloning_strategy.add_primer(primer)

    pcr_left_source = PCRSource(id=0, output_name='left_homology_arm')
    resp = await pcr(pcr_left_source, [locus_seq], [primers[0], primers[1]], 17, 0)
    pcr_product_left: TextFileSequence = resp['sequences'][0]
    pcr_left_source: PCRSource = resp['sources'][0]
    cloning_strategy.add_source_and_sequence(pcr_left_source, pcr_product_left)

    pcr_right_source = PCRSource(id=0, output_name='right_homology_arm')
    resp = await pcr(pcr_right_source, [locus_seq], [primers[2], primers[3]], 17, 0)
    pcr_product_right: TextFileSequence = resp['sequences'][0]
    pcr_right_source: PCRSource = resp['sources'][0]
    cloning_strategy.add_source_and_sequence(pcr_right_source, pcr_product_right)

    with open(os.path.join(os.path.dirname(__file__), 'barcode.gb'), 'rb') as f:
        dummy_resp = Response()
        resp = await read_from_file(dummy_resp, UploadFile(file=f, filename='barcode.gb'), None, None, True, 'barcode')

    barcode_source = resp['sources'][0]
    barcode_seq: TextFileSequence = resp['sequences'][0]
    cloning_strategy.add_source_and_sequence(barcode_source, barcode_seq)

    with open(os.path.join(os.path.dirname(__file__), 'common_plasmid.gb'), 'rb') as f:
        dummy_resp = Response()
        resp = await read_from_file(
            dummy_resp, UploadFile(file=f, filename='common_plasmid.gb'), None, None, True, 'common_plasmid'
        )

    common_plasmid_source = resp['sources'][0]
    common_plasmid_seq: TextFileSequence = resp['sequences'][0]
    cloning_strategy.add_source_and_sequence(common_plasmid_source, common_plasmid_seq)

    golgen_gate = RestrictionAndLigationSource(id=0, output_name='golgen_gate', restriction_enzymes=['BsaI'])
    resp = await restriction_and_ligation(
        golgen_gate, [pcr_product_left, pcr_product_right, barcode_seq, common_plasmid_seq], False, True
    )
    golgen_gate_product: TextFileSequence = resp['sequences'][0]
    golgen_gate_source: RestrictionAndLigationSource = resp['sources'][0]
    cloning_strategy.add_source_and_sequence(golgen_gate_source, golgen_gate_product)

    homologous_recombination_source = HomologousRecombinationSource(id=0, output_name='engineered_allele')
    resp = await homologous_recombination(homologous_recombination_source, [locus_seq, golgen_gate_product], 17)

    multi_site_sources = [
        i
        for i, s in enumerate(resp['sources'])
        if all(join.left_location != join.right_location for join in s.assembly)
    ]
    if len(multi_site_sources) > 1:
        raise ValueError('Multiple insertions possible')

    possible_sources = [resp['sources'][i] for i in multi_site_sources]
    possible_sequences = [resp['sequences'][i] for i in multi_site_sources]
    homologous_recombination_product: TextFileSequence = possible_sequences[0]
    homologous_recombination_source: HomologousRecombinationSource = possible_sources[0]
    cloning_strategy.add_source_and_sequence(homologous_recombination_source, homologous_recombination_product)

    # Save cloning strategy to JSON file
    with open('cloning_strategy.json', 'w') as f:
        json.dump(cloning_strategy.model_dump(), f, indent=2)
    return cloning_strategy


asyncio.run(main())
