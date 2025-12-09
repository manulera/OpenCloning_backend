import os
from pydna.assembly2 import homologous_recombination_integration, pcr_assembly
from opencloning.dna_functions import request_from_addgene
from opencloning.ncbi_requests import get_annotations_from_query, get_genome_region_from_annotation
import asyncio
from Bio import SeqIO
from pydna.primer import Primer
from pydna.opencloning_models import CloningStrategy
from fastapi.datastructures import UploadFile
from pydna.parsers import parse as pydna_parse


async def main(
    gene: str,
    assembly_accession: str,
    output_dir: str,
    plasmid_input: UploadFile | str = '19343',
    padding: int = 1000,
):
    print(f"\033[92mCloning {gene}\033[0m")
    # Parse primers =================================================================================
    primers = [Primer(p) for p in SeqIO.parse(os.path.join(output_dir, gene, 'primers.fa'), 'fasta')]
    common_primers = [Primer(p) for p in SeqIO.parse(os.path.join(output_dir, 'checking_primers.fa'), 'fasta')]

    # Get plasmid sequence =================================================================================
    if isinstance(plasmid_input, UploadFile):
        file_content = (await plasmid_input.read()).decode()

        plasmid = pydna_parse(file_content)[0]
    else:
        plasmid = await request_from_addgene(plasmid_input)

    # Get genome region =====================================================================
    annotations = await get_annotations_from_query(gene, assembly_accession)
    if len(annotations) == 0:
        raise ValueError(f'No annotations found for {gene}')

    annotations = [a for a in annotations if gene.upper() in a['locus_tag'].upper()]
    if len(annotations) != 1:
        raise ValueError(f'No right annotation found for {gene}')

    locus = await get_genome_region_from_annotation(annotations[0], 1000, 1000)

    # PCR ================================================================================================
    pcr_products = pcr_assembly(plasmid, primers[0], primers[1], limit=14, mismatches=0)
    pcr_products[0].name = 'amplified_marker'
    alleles = homologous_recombination_integration(locus, [pcr_products[0]], 40)
    pcr_check1 = pcr_assembly(alleles[0], primers[2], common_primers[1], limit=14, mismatches=0)[0]
    pcr_check1.name = 'check_pcr_left'
    pcr_check2 = pcr_assembly(alleles[0], primers[3], common_primers[0], limit=14, mismatches=0)[0]
    pcr_check2.name = 'check_pcr_right'

    cs = CloningStrategy.from_dseqrecords([pcr_check1, pcr_check2])

    if not os.path.exists(os.path.join(output_dir, gene)):
        os.makedirs(os.path.join(output_dir, gene))

    with open(os.path.join(output_dir, gene, 'cloning_strategy.json'), 'w') as f:
        f.write(cs.model_dump_json(indent=2))


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
