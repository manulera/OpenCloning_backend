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
    UploadedFileSource,
    RestrictionAndLigationSource,
    BaseCloningStrategy,
    HomologousRecombinationSource,
)
from .primer_design_settings import amanda_settings
from ...endpoints.external_import import genome_coordinates, read_from_file
from ...endpoints.assembly import pcr, restriction_and_ligation, homologous_recombination
from ...dna_functions import read_dsrecord_from_json

padding = 0

example_left = GenomeCoordinatesSource(
    id=1,
    assembly_accession='GCF_002847425.1',
    sequence_accession='NZ_CP025534.1',
    start=2699639 - padding,
    end=2700602 + padding,
    strand=1,
)

example_right = GenomeCoordinatesSource(
    id=3,
    assembly_accession='GCF_002847425.1',
    sequence_accession='NZ_CP025534.1',
    start=2700775 - padding,
    end=2701674 + padding,
    strand=1,
)

adapter_left_fwd = 'ataGGTCTCtGGAG'
adapter_left_rvs = 'ataGGTCTCtCATT'
adapter_right_fwd = 'ataGGTCTCtGCTT'
adapter_right_rvs = 'ataGGTCTCtAGCG'


async def design_primers_ebic(seq: str):

    report = bindings.design_primers(
        seq_args={
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': seq,
        },
        global_args=amanda_settings,
    )

    return report


async def main(source_left: GenomeCoordinatesSource, source_right: GenomeCoordinatesSource):

    locus_left = await genome_coordinates(source_left)
    locus_left_seq: TextFileSequence = TextFileSequence.model_validate(locus_left['sequences'][0])
    locus_left_dseqr = read_dsrecord_from_json(locus_left_seq)
    locus_left_seq.id = 2
    locus_left_source: GenomeCoordinatesSource = GenomeCoordinatesSource.model_validate(locus_left['sources'][0])
    locus_left_source.output = 2

    locus_right = await genome_coordinates(source_right)
    locus_right_seq: TextFileSequence = TextFileSequence.model_validate(locus_right['sequences'][0])
    locus_right_dseqr = read_dsrecord_from_json(locus_right_seq)
    locus_right_seq.id = 4
    locus_right_source: GenomeCoordinatesSource = GenomeCoordinatesSource.model_validate(locus_right['sources'][0])
    locus_right_source.output = 4

    report_left = await design_primers_ebic(str(locus_left_dseqr.seq))
    report_right = await design_primers_ebic(str(locus_right_dseqr.seq))

    primer_names = ['left_fwd', 'left_rvs', 'right_fwd', 'right_rvs']
    primer_seqs = [
        adapter_left_fwd + report_left['PRIMER_LEFT'][0]['SEQUENCE'],
        adapter_left_rvs + report_left['PRIMER_RIGHT'][0]['SEQUENCE'],
        adapter_right_fwd + report_right['PRIMER_LEFT'][0]['SEQUENCE'],
        adapter_right_rvs + report_right['PRIMER_RIGHT'][0]['SEQUENCE'],
    ]
    primer_ids = [1, 2, 3, 4]
    primers = [PrimerModel(id=primer_ids[i], name=primer_names[i], sequence=primer_seqs[i]) for i in range(4)]

    pcr_left_source = PCRSource(id=5, output_name='left_homology_arm')
    resp = await pcr(pcr_left_source, [locus_left_seq], [primers[0], primers[1]], 17, 0)
    pcr_product_left: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    pcr_product_left.id = 6
    pcr_left_source: PCRSource = PCRSource.model_validate(resp['sources'][0])
    pcr_left_source.output = 6

    pcr_right_source = PCRSource(id=7, output_name='right_homology_arm')
    resp = await pcr(pcr_right_source, [locus_right_seq], [primers[2], primers[3]], 17, 0)
    pcr_product_right: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    pcr_product_right.id = 8
    pcr_right_source: PCRSource = PCRSource.model_validate(resp['sources'][0])
    pcr_right_source.output = 8

    with open(os.path.join(os.path.dirname(__file__), 'barcode.gb'), 'rb') as f:
        dummy_resp = Response()
        resp = await read_from_file(dummy_resp, UploadFile(file=f, filename='barcode.gb'), None, None, True, 'barcode')

    barcode_source = UploadedFileSource.model_validate(resp['sources'][0])
    barcode_source.id = 9
    barcode_source.output = 10
    barcode_seq: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    barcode_seq.id = 10

    with open(os.path.join(os.path.dirname(__file__), 'common_plasmid.gb'), 'rb') as f:
        dummy_resp = Response()
        resp = await read_from_file(
            dummy_resp, UploadFile(file=f, filename='common_plasmid.gb'), None, None, True, 'common_plasmid'
        )

    common_plasmid_source = UploadedFileSource.model_validate(resp['sources'][0])
    common_plasmid_source.id = 11
    common_plasmid_source.output = 12
    common_plasmid_seq: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    common_plasmid_seq.id = 12

    golgen_gate = RestrictionAndLigationSource(id=13, output_name='golgen_gate', restriction_enzymes=['BsaI'])
    resp = await restriction_and_ligation(
        golgen_gate, [pcr_product_left, pcr_product_right, barcode_seq, common_plasmid_seq], False, True
    )
    golgen_gate_product: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    golgen_gate_product.id = 14
    golgen_gate_source: RestrictionAndLigationSource = RestrictionAndLigationSource.model_validate(resp['sources'][0])
    golgen_gate_source.output = 14

    insertion_locus = GenomeCoordinatesSource(
        id=15,
        assembly_accession='GCF_002847425.1',
        sequence_accession='NZ_CP025534.1',
        start=locus_left_source.start,
        end=locus_right_source.end,
        strand=1,
    )
    resp = await genome_coordinates(insertion_locus)
    insertion_locus_seq: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    insertion_locus_seq.id = 16
    insertion_locus_source: GenomeCoordinatesSource = GenomeCoordinatesSource.model_validate(resp['sources'][0])
    insertion_locus_source.output = 16

    homologous_recombination_source = HomologousRecombinationSource(id=17, output_name='engineered_allele')
    resp = await homologous_recombination(
        homologous_recombination_source, [insertion_locus_seq, golgen_gate_product], 17
    )
    homologous_recombination_product: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    homologous_recombination_product.id = 18
    homologous_recombination_source: HomologousRecombinationSource = HomologousRecombinationSource.model_validate(
        resp['sources'][0]
    )
    homologous_recombination_source.output = 18

    cloning_strategy = BaseCloningStrategy(
        sequences=[
            locus_left_seq,
            locus_right_seq,
            pcr_product_left,
            pcr_product_right,
            barcode_seq,
            common_plasmid_seq,
            golgen_gate_product,
            insertion_locus_seq,
            homologous_recombination_product,
        ],
        sources=[
            locus_left_source,
            locus_right_source,
            pcr_left_source,
            pcr_right_source,
            barcode_source,
            common_plasmid_source,
            golgen_gate_source,
            insertion_locus_source,
            homologous_recombination_source,
        ],
        primers=primers,
        description='EBIC cloning strategy',
    )

    # Save cloning strategy to JSON file
    with open('cloning_strategy.json', 'w') as f:
        json.dump(cloning_strategy.model_dump(), f, indent=2)
    return cloning_strategy


asyncio.run(main(example_left, example_right))
