from bs4 import BeautifulSoup
import asyncio
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ...http_client import get_http_client, Response

default_settings_primer_design = {
    'length': 80,
    'increment': 40,
    'add_seq': 400,
    'opt_len': 22,
    'min_len': 20,
    'max_len': 28,
    'opt_tm': 60.0,
    'min_tm': 57.0,
    'max_tm': 63.0,
    'min_gc': 30,
    'max_gc': 60,
    '.submit': 'Submit',
}


async def async_post(url, headers, data, params=None) -> Response:
    async with get_http_client() as client:
        return await client.post(url, headers=headers, data=data, params=params, timeout=20.0)


async def get_primers(gene):
    print(f"\033[92mGetting primers for {gene}\033[0m")
    # A first request to access the primers
    data = {
        'gene': gene,
        'length': 80,
        'plasmid': 'pFA6a',
        'increment': 40,
        '.submit': 'Submit',
        '.cgifields': 'plasmid',
    }

    url = 'http://bahlerweb.cs.ucl.ac.uk/cgi-bin/PPPP/pppp_deletion.pl'
    headers = None
    resp = await async_post(url, headers, data)
    # Parse with BeautifulSoup
    soup = BeautifulSoup(resp.text, 'html.parser')
    # select forward and reverse primers by default
    forward_primer = soup.find('input', {'name': 'for_sel', 'checked': True})
    reverse_primer = soup.find('input', {'name': 'rev_sel', 'checked': True})
    forward_primer_seq = re.sub(r'[^a-zA-Z]', '', forward_primer['value'])
    reverse_primer_seq = re.sub(r'[^a-zA-Z]', '', reverse_primer['value'])
    # Make a second request to get the checking primers
    data = {
        'gene': gene,
        'for_sel': forward_primer['value'],
        'rev_sel': reverse_primer['value'],
    }
    data.update(default_settings_primer_design)

    url = 'http://bahlerweb.cs.ucl.ac.uk/cgi-bin/PPPP/pppp_checking.pl'
    resp = await async_post(url, headers, data)
    # Parse with BeautifulSoup
    soup = BeautifulSoup(resp.text, 'html.parser')
    # Find a pre tag with the text "Left Primer:"
    left_check_primer = soup.find('pre', string=re.compile(r'Left Primer:\s+Sequence:')).get_text().strip()
    right_check_primer = soup.find('pre', string=re.compile(r'Right Primer:\s+Sequence:')).get_text().strip()

    pattern = r'Sequence:\s+(\S+)'
    left_check_primer_seq = re.search(pattern, left_check_primer).group(1)
    right_check_primer_seq = re.search(pattern, right_check_primer).group(1)
    return forward_primer_seq, reverse_primer_seq, left_check_primer_seq, right_check_primer_seq


async def main(gene: str, output_dir: str):
    gene_dir = os.path.join(output_dir, gene)
    os.makedirs(gene_dir, exist_ok=True)

    forward_primer, reverse_primer, forward_check, reverse_check = await get_primers(gene)

    primers = [
        SeqRecord(Seq(forward_primer), id=f"{gene}_fwd", description=''),
        SeqRecord(Seq(reverse_primer), id=f"{gene}_rvs", description=''),
        SeqRecord(Seq(forward_check), id=f"{gene}_fwd_check", description=''),
        SeqRecord(Seq(reverse_check), id=f"{gene}_rvs_check", description=''),
    ]

    with open(os.path.join(gene_dir, 'primers.fa'), 'w') as f:
        SeqIO.write(primers, f, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get primers for S. pombe genes')
    parser.add_argument(
        '--genes', type=str, required=True, help='Path to a file containing a list of genes, one per line'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        default='batch_cloning_output',
        help='Directory to save the output files (default: batch_cloning_output)',
    )

    args = parser.parse_args()

    with open(args.genes, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    for gene in genes:
        asyncio.run(main(gene, args.output_dir))
