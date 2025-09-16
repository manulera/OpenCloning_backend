import os
from opencloning.endpoints.external_import import genome_coordinates, get_from_repository_id_addgene, read_from_file
from opencloning.endpoints.assembly import pcr, homologous_recombination
from opencloning.pydantic_models import (
    GenomeCoordinatesSource,
    TextFileSequence,
    AddgeneIdSource,
    PCRSource,
    PrimerModel,
    HomologousRecombinationSource,
    BaseCloningStrategy,
    UploadedFileSource,
)

from opencloning.ncbi_requests import get_annotations_from_query
import asyncio
import json
from Bio import SeqIO
from pydna.parsers import parse as pydna_parse


async def main(
    gene: str, assembly_accession: str, output_dir: str, plasmid: str | dict = '19343', padding: int = 1000
):
    print(f"\033[92mCloning {gene}\033[0m")
    # Parse primers =================================================================================
    primer_records = list(SeqIO.parse(os.path.join(output_dir, gene, 'primers.fa'), 'fasta'))
    checking_primers = list(SeqIO.parse(os.path.join(output_dir, 'checking_primers.fa'), 'fasta'))
    primer_records = primer_records[:3] + checking_primers[1:] + primer_records[3:] + checking_primers[:1]
    primers = []
    for primer in primer_records:
        primers.append(PrimerModel(sequence=str(primer.seq), id=0, name=primer.id))

    # Get genome region =====================================================================
    annotations = await get_annotations_from_query(gene, assembly_accession)
    if len(annotations) == 0:
        raise ValueError(f'No annotations found for {gene}')

    annotations = [a for a in annotations if gene.upper() in a['locus_tag'].upper()]
    if len(annotations) != 1:
        raise ValueError(f'No right annotation found for {gene}')

    annotation = annotations[0]

    gene_range = annotation['genomic_regions'][0]['gene_range']['range'][0]
    sequence_accession = annotation['genomic_regions'][0]['gene_range']['accession_version']
    locus_tag = annotation.get('locus_tag', None)
    gene_id = annotation.get('gene_id', None)
    start = int(gene_range['begin'])
    end = int(gene_range['end'])
    orientation = 1 if gene_range['orientation'] == 'plus' else -1

    source = GenomeCoordinatesSource(
        id=0,
        start=start - padding,
        end=end + padding,
        strand=orientation,
        assembly_accession=assembly_accession,
        sequence_accession=sequence_accession,
        locus_tag=locus_tag,
        gene_id=gene_id,
        output_name=gene,
    )
    locus = await genome_coordinates(source)

    cloning_strategy = BaseCloningStrategy(
        sequences=[],
        sources=[],
        primers=[],
        description=f'Cloning strategy for deleting the gene {gene} using PCR and homologous recombination',
    )
    for primer in primers:
        cloning_strategy.add_primer(primer)
    locus_seq: TextFileSequence = TextFileSequence.model_validate(locus['sequences'][0])
    locus_source: GenomeCoordinatesSource = GenomeCoordinatesSource.model_validate(locus['sources'][0])
    cloning_strategy.add_source_and_sequence(locus_source, locus_seq)

    # Get plasmid sequence =================s================================================================
    if not isinstance(plasmid, str):
        if plasmid.filename.endswith('.fa') or plasmid.filename.endswith('.fasta'):
            resp = await read_from_file(plasmid, None, None, True, None, None, None)
        else:
            resp = await read_from_file(plasmid, None, None, None, None, None, None)
        # Verify that plasmid is circular
        if not pydna_parse(resp['sequences'][0].file_content)[0].circular:
            raise ValueError('Plasmid is not circular')
        plasmid_source: UploadedFileSource = UploadedFileSource.model_validate(resp['sources'][0])
    else:
        addgene_source = AddgeneIdSource(
            id=0,
            repository_id=plasmid,
            repository_name='addgene',
        )
        resp = await get_from_repository_id_addgene(addgene_source)
        plasmid_source: AddgeneIdSource = AddgeneIdSource.model_validate(resp['sources'][0])

    plasmid_seq: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    cloning_strategy.add_source_and_sequence(plasmid_source, plasmid_seq)

    # PCR ================================================================================================
    pcr_source = PCRSource(id=0, output_name='amplified_marker')
    resp = await pcr(pcr_source, [plasmid_seq], [primers[0], primers[1]], 15, 0)

    pcr_product: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    pcr_source: PCRSource = PCRSource.model_validate(resp['sources'][0])
    cloning_strategy.add_source_and_sequence(pcr_source, pcr_product)

    # Homologous recombination ========================================================================
    hrec_source = HomologousRecombinationSource(id=0, output_name='deletion_allele')
    resp = await homologous_recombination(hrec_source, [locus_seq, pcr_product], 50)

    hrec_product: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    hrec_source: HomologousRecombinationSource = HomologousRecombinationSource.model_validate(resp['sources'][0])
    cloning_strategy.add_source_and_sequence(hrec_source, hrec_product)

    # Checking pcr 1 ======================================================================================
    check_pcr_source_left = PCRSource(id=0, output_name='check_pcr_left')
    resp = await pcr(check_pcr_source_left, [hrec_product], [primers[2], primers[3]], 15, 0)

    check_pcr_product_left: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    check_pcr_source_left: PCRSource = PCRSource.model_validate(resp['sources'][0])
    cloning_strategy.add_source_and_sequence(check_pcr_source_left, check_pcr_product_left)

    # Checking pcr 2 ======================================================================================
    check_pcr_source_right = PCRSource(id=0, output_name='check_pcr_right')
    resp = await pcr(check_pcr_source_right, [hrec_product], [primers[4], primers[5]], 15, 0)

    check_pcr_product_right: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    check_pcr_source_right: PCRSource = PCRSource.model_validate(resp['sources'][0])
    cloning_strategy.add_source_and_sequence(check_pcr_source_right, check_pcr_product_right)

    cloning_strategy.description = (
        f'Cloning strategy for deleting the gene {gene} using PCR and homologous recombination'
    )

    if not os.path.exists(os.path.join(output_dir, gene)):
        os.makedirs(os.path.join(output_dir, gene))

    with open(os.path.join(output_dir, gene, 'cloning_strategy.json'), 'w') as f:
        json.dump(cloning_strategy.model_dump(), f, indent=2)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='List of genes to delete from S. pombe')
    parser.add_argument(
        '--genes', type=str, required=True, help='Path to a file containing a list of genes, one per line'
    )
    args = parser.parse_args()

    parser.add_argument(
        '--assembly_accession',
        type=str,
        default='GCF_000002945.2',
        help='Assembly accession for S. pombe genome (default: GCF_000002945.2)',
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        default='batch_cloning_output',
        help='Directory to save the output files (default: batch_cloning_output)',
    )

    parser.add_argument(
        '--plasmid',
        type=str,
        default='19343',
        help='Addgene ID for the plasmid (default: 19343)',
    )

    args = parser.parse_args()
    assembly_accession = args.assembly_accession

    with open(args.genes, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    for gene in genes:
        asyncio.run(main(gene, assembly_accession, args.output_dir, args.plasmid))
