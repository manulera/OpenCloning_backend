"""
Some extra tests that are not covered by the main endpoint tests
"""

import opencloning.ncbi_requests as ncbi_requests
import pytest
import respx
from fastapi import HTTPException
import unittest


class NcbiAsyncRequestsTest(unittest.IsolatedAsyncioTestCase):

    @respx.mock
    async def test_get_genbank_sequence_subset(self):
        respx.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi').respond(503, text='')
        with pytest.raises(HTTPException) as e:
            await ncbi_requests.get_genbank_sequence('blah', 1, 10, 1)
        assert e.value.status_code == 503
        assert e.value.detail == 'NCBI returned an error'

        respx.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi').respond(500, text='')
        with pytest.raises(HTTPException) as e:
            await ncbi_requests.get_genbank_sequence('blah', 1, 10, 1)
        assert e.value.status_code == 500
        assert e.value.detail == 'NCBI is down, try again later'

        respx.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi').respond(504, text='')
        with pytest.raises(HTTPException) as e:
            await ncbi_requests.get_genbank_sequence('blah', 1, 10, 1)
        assert e.value.status_code == 500
        assert e.value.detail == 'NCBI returned an unexpected error'

    async def test_get_annotations_from_query(self):
        result = await ncbi_requests.get_annotations_from_query('SPAPB1A10.09', 'GCF_000002945.2')
        self.assertEqual(result[0]['symbol'], 'ase1')

    @respx.mock
    async def test_get_annotations_from_query_errors(self):
        def get_url(assembly_accession):
            return f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/annotation_report'

        respx.get(get_url('blah')).respond(404, text='')
        with pytest.raises(HTTPException) as e:
            await ncbi_requests.get_annotations_from_query('bluh', 'blah')
        assert e.value.status_code == 404
        assert e.value.detail == 'wrong accession number'

        respx.get(get_url('blah2')).respond(200, json={'dummy': 'data'})
        with pytest.raises(HTTPException) as e:
            await ncbi_requests.get_annotations_from_query('bluh2', 'blah2')
        assert e.value.status_code == 404
        assert e.value.detail == 'query "bluh2" gave no results'

    @respx.mock
    async def test_get_annotation_from_locus_tag_errors(self):
        def get_url(assembly_accession):
            return f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/annotation_report'

        respx.get(get_url('blah3')).respond(404, text='')
        with pytest.raises(HTTPException) as e:
            await ncbi_requests.get_annotation_from_locus_tag('bluh3', 'blah3')
        assert e.value.status_code == 404
        assert e.value.detail == 'wrong accession number'

        respx.get(get_url('blah4')).respond(200, json={'dummy': 'data'})
        with pytest.raises(HTTPException) as e:
            await ncbi_requests.get_annotation_from_locus_tag('bluh4', 'blah4')
        assert e.value.status_code == 404
        assert e.value.detail == 'wrong locus_tag'

        respx.get(get_url('blah5')).respond(
            200, json={'reports': [{'annotation': {'locus_tag': 'bluh5'}}, {'annotation': {'locus_tag': 'bluh5'}}]}
        )
        with pytest.raises(HTTPException) as e:
            await ncbi_requests.get_annotation_from_locus_tag('bluh5', 'blah5')
        assert e.value.status_code == 400
        assert e.value.detail == 'multiple matches for locus_tag'

    async def test_empty_get_assembly_accession_from_sequence_accession(self):
        # For accessions that are not linked to assemblies
        self.assertEqual([], await ncbi_requests.get_assembly_accession_from_sequence_accession('DQ208311.2'))

    @pytest.mark.xfail(reason='waiting on https://github.com/ncbi/datasets/issues/380#issuecomment-2231142816')
    async def test_get_assembly_accession_from_sequence_accession(self):
        # For accessions that are linked to assemblies
        self.assertEqual(
            ['GCF_000002945.2', 'GCF_000002945.1'],
            await ncbi_requests.get_assembly_accession_from_sequence_accession('NC_003424.3'),
        )
