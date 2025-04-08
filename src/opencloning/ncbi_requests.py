import requests
from fastapi import HTTPException
from pydna.parsers import parse as pydna_parse
from httpx import AsyncClient, Response
from pydna.dseqrecord import Dseqrecord
from .app_settings import settings

headers = None if settings.NCBI_API_KEY is None else {'api_key': settings.NCBI_API_KEY}


async def async_get(url, headers, params=None) -> Response:
    async with AsyncClient(timeout=20.0) as client:
        return await client.get(url, headers=headers, params=params)


# TODO: this does not return old assembly accessions, see https://github.com/ncbi/datasets/issues/380#issuecomment-2231142816
def get_assembly_accession_from_sequence_accession(sequence_accession: str) -> list[str]:
    """Get the assembly accession from a sequence accession"""

    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/sequence_accession/{sequence_accession}/sequence_assemblies'
    resp = requests.get(url, headers=headers)
    data = resp.json()
    if 'accessions' in data:
        return data['accessions']
    else:
        return []


async def get_sequence_accessions_from_assembly_accession(assembly_accession: str) -> list[str]:
    """Get the sequence accessions from an assembly accession"""
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/sequence_reports'
    resp = await async_get(url, headers=headers)
    data = resp.json()
    if 'reports' in data:
        refseq_accessions = [report['refseq_accession'] for report in data['reports'] if 'refseq_accession' in report]
        genbank_accessions = [
            report['genbank_accession'] for report in data['reports'] if 'genbank_accession' in report
        ]
        return refseq_accessions + genbank_accessions
    elif 'total_count' in data:
        raise HTTPException(400, f'No sequence accessions linked, see {url}')
    else:
        raise HTTPException(404, 'Wrong assembly accession number')


async def get_annotation_from_locus_tag(locus_tag: str, assembly_accession: str) -> dict:
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/annotation_report?search_text={locus_tag}'
    resp = await async_get(url, headers=headers)
    if resp.status_code == 404:
        raise HTTPException(404, 'wrong accession number')
    data = resp.json()
    if 'reports' not in data:
        raise HTTPException(404, 'wrong locus_tag')

    matching_annotations = list(a['annotation'] for a in data['reports'] if a['annotation']['locus_tag'] == locus_tag)

    if len(matching_annotations) == 0:
        raise HTTPException(404, 'wrong locus_tag')
    elif len(matching_annotations) > 1:
        # Not sure if this can ever happen, but just in case
        raise HTTPException(400, 'multiple matches for locus_tag')

    return matching_annotations[0]


async def get_annotations_from_query(query: str, assembly_accession: str) -> list[dict]:
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/annotation_report?search_text={query}'
    resp = await async_get(url, headers=headers)
    if resp.status_code == 404:
        raise HTTPException(404, 'wrong accession number')

    data = resp.json()
    if 'reports' not in data:
        raise HTTPException(404, f'query "{query}" gave no results')

    if len(data['reports']) > 1:
        raise HTTPException(400, 'multiple matches for query')

    return [r['annotation'] for r in data['reports']]


async def get_sequence_length_from_sequence_accession(sequence_accession: str) -> int:
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
    params = {'id': sequence_accession, 'db': 'nuccore', 'retmode': 'json'}
    if headers is not None:
        params['api_key'] = headers['api_key']
    resp = await async_get(url, headers=headers, params=params)
    data = resp.json()
    if 'result' not in data:
        raise HTTPException(503, 'NCBI returned an error (try again)')
    if len(data['result']['uids']) == 0:
        raise HTTPException(404, 'wrong sequence accession')
    sequence_id = data['result']['uids'][0]
    return data['result'][sequence_id]['slen']


async def get_genbank_sequence(sequence_accession, start=None, end=None, strand=None) -> Dseqrecord:
    gb_strand = 1 if strand == 1 or strand is None else 2
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    params = {
        'db': 'nuccore',
        'id': sequence_accession,
        'rettype': 'gbwithparts',
        'seq_start': start,
        'seq_stop': end,
        'strand': gb_strand,
        'retmode': 'text',
    }
    if headers is not None:
        params['api_key'] = headers['api_key']

    resp = await async_get(url, headers=headers, params=params)
    if resp.status_code == 200:
        return pydna_parse(resp.text)[0]
    elif resp.status_code == 400:
        raise HTTPException(404, 'wrong sequence accession')
    elif resp.status_code == 503:
        raise HTTPException(503, 'NCBI returned an error')
    else:
        raise HTTPException(500, 'NCBI returned an unexpected error')


def validate_coordinates_pre_request(start, end, strand):
    # TODO: move this to the class
    if strand not in [1, -1]:
        raise HTTPException(422, 'strand must be 1 or -1')
    if start >= end:
        raise HTTPException(422, 'start must be less than end')
    if start < 1:
        raise HTTPException(422, 'start must be greater than 0')
    if end - start > 100000:
        raise HTTPException(400, 'sequence is too long (max 100000 bp)')
