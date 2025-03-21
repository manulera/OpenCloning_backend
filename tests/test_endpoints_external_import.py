from fastapi.testclient import TestClient
from pydna.dseqrecord import Dseqrecord
import unittest
import tempfile
import pytest
import os
import respx
import httpx

import opencloning.request_examples as request_examples
from opencloning.dna_functions import read_dsrecord_from_json
import opencloning.main as _main
from opencloning.pydantic_models import (
    RepositoryIdSource,
    TextFileSequence,
    UploadedFileSource,
    GenomeCoordinatesSource,
    AddgeneIdSource,
    BenchlingUrlSource,
    EuroscarfSource,
    SnapGenePlasmidSource,
    IGEMSource,
    WekWikGeneIdSource,
    SEVASource,
)


test_files = os.path.join(os.path.dirname(__file__), 'test_files')

client = TestClient(_main.app)


class ReadFileTest(unittest.TestCase):
    def test_read_files(self):
        """Test that uploading files with single and multiple sequences works."""

        examples = [
            {'file': f'{test_files}/pFA6a-hphMX6.gb', 'format': 'genbank', 'nb_sequences': 1},
            {'file': f'{test_files}/dummy_EcoRI.fasta', 'format': 'fasta', 'nb_sequences': 1},
            {'file': f'{test_files}/dummy_multi_fasta.fasta', 'format': 'fasta', 'nb_sequences': 2},
            {
                'file': f'{test_files}/addgene-plasmid-39296-sequence-49545.dna',
                'format': 'snapgene',
                'nb_sequences': 1,
            },
            {'file': f'{test_files}/ase1.embl', 'format': 'embl', 'nb_sequences': 1},
            # Ape files as of 2024-10-30 did not have a properly formatted LOCUS line
            {
                'file': f'{test_files}/example.ape',
                'format': 'genbank',
                'nb_sequences': 1,
                'warning': True,
                'circular': False,
            },
            # Euroscarf files as of 2024-10-30 did not have a properly formatted LOCUS line
            {
                'file': f'{test_files}/pKT128_euroscarf.gb',
                'format': 'genbank',
                'nb_sequences': 1,
                'warning': True,
                'circular': True,
            },
        ]

        for example in examples:
            with open(example['file'], 'rb') as f:
                response = client.post('/read_from_file', files={'file': f})
            self.assertEqual(response.status_code, 200)
            payload = response.json()

            resulting_sequences = [
                read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
            ]
            sources = [UploadedFileSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(sources), example['nb_sequences'])
            self.assertEqual(len(resulting_sequences), example['nb_sequences'])
            for seq in resulting_sequences:
                self.assertGreater(len(seq), 3)
            for source in sources:
                self.assertEqual(source.sequence_file_format, example['format'])

            if 'warning' in example and example['warning']:
                self.assertIn('x-warning', response.headers)
            else:
                self.assertNotIn('x-warning', response.headers)
            if 'circular' in example:
                for seq in resulting_sequences:
                    self.assertEqual(seq.circular, example['circular'])

        # Test naming
        with open(example['file'], 'rb') as f:
            response = client.post('/read_from_file?output_name=hello', files={'file': f})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        dseqr = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(dseqr.name, 'hello')

    def test_errors_read_files(self):
        # Create a temp empty file
        with tempfile.NamedTemporaryFile() as temp_empty_file:

            examples = [
                {
                    'file': './test_endpoints.py',
                    'format': None,
                    'error_message': 'could not guess',
                },
                {
                    'file': './test_endpoints.py',
                    'format': 'snapgene',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': './test_endpoints.py',
                    'format': 'genbank',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': './test_endpoints.py',
                    'format': 'embl',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': './test_endpoints.py',
                    'format': 'fasta',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': temp_empty_file.name,
                    'format': 'genbank',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': temp_empty_file.name,
                    'format': 'embl',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': temp_empty_file.name,
                    'format': 'fasta',
                    'error_message': 'Biopython cannot process',
                },
            ]
            for example in examples:
                if example['format'] is not None:
                    response = client.post(
                        f'/read_from_file?sequence_file_format={example["format"]}',
                        files={'file': temp_empty_file},
                    )
                else:
                    response = client.post('/read_from_file', files={'file': temp_empty_file})

                self.assertNotEqual(response.status_code, 200)
                self.assertIn(example['error_message'], response.json()['detail'])

    def test_file_index_known(self):
        """Test that if the index in file is specified it works."""

        with open(f'{test_files}/dummy_multi_fasta.fasta', 'rb') as f:
            response = client.post('/read_from_file?index_in_file=1', files={'file': f})

        self.assertEqual(response.status_code, 200)
        payload = response.json()

        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [UploadedFileSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(sources), 1)
        self.assertEqual(len(resulting_sequences), 1)

        # If the index is outside the range, it should raise an error
        with open(f'{test_files}/dummy_multi_fasta.fasta', 'rb') as f:
            response = client.post('/read_from_file?index_in_file=2', files={'file': f})
        self.assertEqual(response.status_code, 404)
        self.assertEqual(response.json()['detail'], 'The index_in_file is out of range.')

    def test_circularize_fasta_sequence(self):
        """Test that the circularize parameter works when reading from file"""
        file_paths = [
            f'{test_files}/dummy_EcoRI.fasta',
            f'{test_files}/example.ape',
            f'{test_files}/gateway_manual_cloning/pcr_product-attP1_1-attP2_1.dna',
        ]
        for file_path in file_paths:
            with open(file_path, 'rb') as f:
                response = client.post('/read_from_file?circularize=True', files={'file': f})

            self.assertEqual(response.status_code, 200)
            payload = response.json()

            resulting_sequences = [
                read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
            ]
            sources = [UploadedFileSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(sources), 1)
            self.assertEqual(len(resulting_sequences), 1)
            self.assertTrue(sources[0].circularize)
            self.assertTrue(resulting_sequences[0].circular)


class GenBankTest(unittest.TestCase):

    # TODO these tests will not work off-line, so the case where connection cannot be established should be handled in some way
    def test_request_gene(self):
        """Test whether the gene is requested from GenBank"""
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='NM_001018957.2',
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequence = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertIn('Ase1', sequence.description)

    def test_request_wrong_id(self):
        """Test a wrong Genbank id"""
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='wrong_id',
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 404)

    # Extremely unlikely case where the request that checks for the length
    # succeeds, but the request that gets the sequence fails
    @respx.mock
    def test_request_wrong_id2(self):
        respx.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi').mock(
            return_value=httpx.Response(200, json={'result': {'uids': ['1'], '1': {'slen': 1000}}})
        )
        # 400 is the error code for a wrong sequence accession :_)
        respx.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi').respond(400, text='')
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='wrong_id',
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 404)

    @respx.mock
    def test_eutils_down(self):
        """Test that the request fails if the NCBI is down"""

        # First request fails
        respx.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi').mock(
            side_effect=httpx.ConnectError('Connection error')
        )
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='NM_001018957.2',
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 504)

        # Second request fails
        respx.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi').mock(
            return_value=httpx.Response(200, json={'result': {'uids': ['1'], '1': {'slen': 1000}}})
        )
        respx.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi').mock(
            side_effect=httpx.ConnectError('Connection error')
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 504)

    def test_redirect(self):
        """The repository_id endpoint should redirect based on repository_name value"""
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='NM_001018957.2',
        )
        response = client.post('/repository_id', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequence: Dseqrecord = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertIn('Ase1', sequence.description)

    def test_rename(self):
        """If passing output_name, it renames the output"""
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='NM_001018957.2',
            output_name='hello',
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequence = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(sequence.name, 'hello')

    def test_long_sequence(self):
        """Test that a long sequence raises an error"""
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='CU329670.1',
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 400)
        self.assertIn('sequence is too long', response.json()['detail'])


class AddgeneTest(unittest.TestCase):
    def test_request_plasmid(self):
        """Test whether the gene is requested from Addgene and returns the right info"""
        examples = [
            {
                'id': '39282',
                'url': 'https://media.addgene.org/snapgene-media/v2.0.0/sequences/240599/4936a6ae-6b4d-4d24-b7ac-2339fad5755d/addgene-plasmid-39282-sequence-240599.gbk',
                'type': 'addgene-full',
                'name': 'pFA6a-kanMX6-P81nmt1',
            },
            {
                'id': '39289',
                'url': 'https://media.addgene.org/snapgene-media/v2.0.0/sequences/49640/4fa9f18b-d5ca-4ac6-a50c-76a6cd15cbab/addgene-plasmid-39289-sequence-49640.gbk',
                'type': 'depositor-full',
                'name': 'pFA6a-kanMX6-P3nmt1-GFP',
            },
        ]
        for example in examples:
            source = AddgeneIdSource(
                id=1,
                repository_name='addgene',
                repository_id=example['id'],
            )

            response = client.post('/repository_id/addgene', json=source.model_dump())
            self.assertEqual(response.status_code, 200)
            payload = response.json()
            resulting_sequences = [
                read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
            ]
            sources = [AddgeneIdSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)
            self.assertEqual(sources[0].addgene_sequence_type, example['type'])
            self.assertEqual(sources[0].sequence_file_url, example['url'])
            self.assertEqual(resulting_sequences[0].name, example['name'])

            # We get the same response when making the response with the url
            response2 = client.post('/repository_id/addgene', json=payload['sources'][0])
            self.assertEqual(response.json(), response2.json())

    @pytest.mark.xfail(reason='This file was removed from Addgene, not sure what the ideal behavior should be')
    def test_old_url(self):
        """Works for an Addgene url that has now been replaced by a newer one"""
        source = AddgeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='65109',
            addgene_sequence_type='addgene-full',
            sequence_file_url='https://media.addgene.org/snapgene-media/v1.7.9-0-g88a3305/sequences/110162/c1c98803-c8ba-44a6-95b8-d6a94097e36f/addgene-plasmid-65109-sequence-110162.gbk',
        )
        response = client.post('/repository_id/addgene', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [AddgeneIdSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)
        self.assertEqual(sources[0].addgene_sequence_type, 'addgene-full')
        self.assertEqual(
            sources[0].sequence_file_url,
            'https://media.addgene.org/snapgene-media/v1.7.9-0-g88a3305/sequences/110162/c1c98803-c8ba-44a6-95b8-d6a94097e36f/addgene-plasmid-65109-sequence-110162.gbk',
        )
        self.assertEqual(resulting_sequences[0].name, 'pYTK002')

    def test_missing_sequences(self):
        # Non-existing id
        source = AddgeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='DUMMYTEST',
        )

        response = client.post('/repository_id/addgene', json=source.model_dump())
        self.assertEqual(response.status_code, 404)
        self.assertIn('wrong addgene id', response.json()['detail'])

        # Id that has no full-sequences
        source = AddgeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='39291',
        )
        response = client.post('/repository_id/addgene', json=source.model_dump())
        self.assertEqual(response.status_code, 404)
        self.assertIn('The requested plasmid does not have full sequences', response.json()['detail'])

        # url does not exist
        source = AddgeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='39282',
            sequence_file_url='https://media.addgene.org/snapgene-media/wrongggggggg.gbk',
        )
        response = client.post('/repository_id/addgene', json=source.model_dump())
        self.assertEqual(response.status_code, 404)

        # TODO url exists but does not match id, or does not match

    def test_redirect(self):
        """Test repository_id endpoint should redirect based on repository_name value"""
        source = AddgeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='39282',
        )
        response = client.post('/repository_id', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequence: Dseqrecord = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertIn('synthetic circular DNA', sequence.description)

    @respx.mock
    def test_addgene_down(self):
        respx.get('https://www.addgene.org/39282/sequences/').mock(side_effect=httpx.ConnectError('Connection error'))
        source = AddgeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='39282',
        )
        response = client.post('/repository_id/addgene', json=source.model_dump())
        self.assertEqual(response.status_code, 504)
        self.assertIn('unable to connect to Addgene', response.json()['detail'])


class WekWikGeneSourceTest(unittest.TestCase):

    def test_valid_id(self):
        source = WekWikGeneIdSource(
            id=1,
            repository_name='wekwikgene',
            repository_id='0000304',
        )
        response = client.post('/repository_id/wekwikgene', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)
        sequence: Dseqrecord = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertTrue(sequence.circular)
        self.assertIn('RPL15', sequence.name)

    def test_invalid_id(self):
        source = WekWikGeneIdSource(
            id=1,
            repository_name='wekwikgene',
            repository_id='999999999999999999999999999999',  # Non-existent ID
        )
        response = client.post('/repository_id/wekwikgene', json=source.model_dump())
        self.assertEqual(response.status_code, 404)
        self.assertIn('invalid wekwikgene id', response.json()['detail'])

    @respx.mock
    def test_wekwikgene_down(self):
        source = WekWikGeneIdSource(
            id=1,
            repository_name='wekwikgene',
            repository_id='0000304',
        )

        respx.get('https://wekwikgene.wllsb.edu.cn/plasmids/0000304').mock(
            side_effect=httpx.ConnectError('Connection error')
        )
        response = client.post('/repository_id/wekwikgene', json=source.model_dump())
        self.assertEqual(response.status_code, 504)
        self.assertIn('unable to connect to WekWikGene', response.json()['detail'])

    def test_redirect(self):
        source = WekWikGeneIdSource(
            id=1,
            repository_name='wekwikgene',
            repository_id='0000304',
        )
        response = client.post('/repository_id', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequence: Dseqrecord = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertIn('RPL15', sequence.name)


class BenchlingUrlSourceTest(unittest.TestCase):

    def test_valid_url(self):
        url = 'https://benchling.com/siverson/f/lib_B94YxDHhQh-cidar-moclo-library/seq_dh1FrJTc-b0015_dh.gb'
        source = BenchlingUrlSource(id=0, repository_id=url, repository_name='benchling')
        response = client.post('/repository_id/benchling', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)
        self.assertEqual(payload['sources'][0], source.model_dump())

    def test_invalid_url(self):
        # We have to initialize the object with a valid url
        url = 'https://benchling.com/siverson/f/lib_B94YxDHhQh-cidar-moclo-library/seq_dh1FrJTc-b0015_dh.gb'
        source_object = BenchlingUrlSource(id=0, repository_id=url, repository_name='benchling')

        # In the dict, we can then edit
        source_dict = source_object.model_dump()

        # url missing the .gb
        source_dict['repository_id'] = (
            'https://benchling.com/siverson/f/lib_B94YxDHhQh-cidar-moclo-library/seq_dh1FrJTc-b0015_dh'
        )
        response = client.post('/repository_id/benchling', json=source_dict)
        self.assertEqual(response.status_code, 422)

        # The /edit url
        source_dict['repository_id'] = (
            'https://benchling.com/siverson/f/lib_B94YxDHhQh-cidar-moclo-library/seq_dh1FrJTc-b0015_dh/edit'
        )
        response = client.post('/repository_id/benchling', json=source_dict)
        self.assertEqual(response.status_code, 422)

        # One that matches the pattern but does not exist
        url = 'https://benchling.com/bluh/blah.gb'
        source = BenchlingUrlSource(id=0, repository_id=url, repository_name='benchling')
        response = client.post('/repository_id/benchling', json=source.model_dump())
        self.assertEqual(response.status_code, 404)
        self.assertIn('file requested from url not found', response.json()['detail'])


class SnapGenePlasmidSourceTest(unittest.TestCase):

    def test_valid_url(self):

        source = SnapGenePlasmidSource(
            id=0, repository_id='basic_cloning_vectors/pEASY-T1_(linearized)', repository_name='snapgene'
        )
        response = client.post('/repository_id/snapgene', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)
        out_source = payload['sources'][0]
        # It should have been renamed
        self.assertEqual(out_source['output_name'], 'pEASY-T1_(linearized)')
        # Remove that and compare with original source
        out_source['output_name'] = None
        self.assertEqual(out_source, source.model_dump())
        seq = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(seq.name, 'pEASY-T1_(linearized)')

        # We can also provide a name
        source2 = SnapGenePlasmidSource(
            id=0,
            repository_id='basic_cloning_vectors/pEASY-T1_(linearized)',
            repository_name='snapgene',
            output_name='my_name',
        )
        response = client.post('/repository_id/snapgene', json=source2.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload['sources'][0]['output_name'], 'my_name')
        seq = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(seq.name, 'my_name')

    def test_invalid_url(self):
        # Invalid url
        source = SnapGenePlasmidSource(id=0, repository_id='hello/world', repository_name='snapgene')
        response = client.post('/repository_id/snapgene', json=source.model_dump())
        self.assertEqual(response.status_code, 404)

        # Wrongly formatted url
        source_dict = source.model_dump()
        source_dict['repository_id'] = 'hello'
        response = client.post('/repository_id/snapgene', json=source_dict)
        self.assertEqual(response.status_code, 422)


class EuroscarfSourceTest(unittest.TestCase):

    def test_valid_url(self):
        source = EuroscarfSource(id=0, repository_id='P30174', repository_name='euroscarf')
        response = client.post('/repository_id/euroscarf', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)
        out_source = payload['sources'][0]
        self.assertEqual(out_source, source.model_dump())
        sequence = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(sequence.name, 'pKT128')
        self.assertEqual(len(sequence), 4738)
        self.assertTrue(any('yEGFP' in f.qualifiers['gene'] for f in sequence.features))

        # Ensure that linear files are circularised
        source = EuroscarfSource(id=0, repository_id='P30555', repository_name='euroscarf')
        response = client.post('/repository_id/euroscarf', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        sequence = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertTrue(sequence.circular)

    def test_invalid_url(self):
        # Compatible with regex, but does not exist
        source = EuroscarfSource(id=0, repository_id='P99999999999999', repository_name='euroscarf')
        response = client.post('/repository_id/euroscarf', json=source.model_dump())
        self.assertEqual(response.status_code, 404)

        # Not compatible with regex
        source_dict = source.model_dump()
        source_dict['repository_id'] = 'hello'
        response = client.post('/repository_id/euroscarf', json=source_dict)
        self.assertEqual(response.status_code, 422)

    @respx.mock
    def test_circularize_plasmid(self):
        # We mock a request in which we would get a linear plasmid
        source = EuroscarfSource(id=0, repository_id='P9999999999999', repository_name='euroscarf')
        respx.get('http://www.euroscarf.de/plasmid_details.php').respond(
            200, text='<html><body><a href="files/dna/test.gb">Download</a></body></html>'
        )
        with open(f'{test_files}/ase1.gb', 'r') as f:
            str_content = f.read()

        respx.get('http://www.euroscarf.de/files/dna/test.gb').respond(200, text=str_content)
        response = client.post('/repository_id/euroscarf', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        seq = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(seq.circular, True)


class IGEMSourceTest(unittest.TestCase):
    good_url = 'https://raw.githubusercontent.com/manulera/annotated-igem-distribution/master/results/plasmids/1.gb'
    no_gb_url = 'https://blah.com/1.txt'
    wrong_url = (
        'https://raw.githubusercontent.com/manulera/annotated-igem-distribution/master/results/plasmids/dummy.gb'
    )

    @pytest.mark.flaky(reruns=3, reruns_delay=2)
    def test_igem(self):
        source = IGEMSource(
            id=0, repository_name='igem', repository_id='BBa_C0062-dummy', sequence_file_url=self.good_url
        )
        response = client.post('/repository_id/igem', json=source.model_dump())
        payload = response.json()
        self.assertEqual(response.status_code, 200)
        self.assertEqual(payload['sources'][0]['repository_id'], 'BBa_C0062-dummy')
        self.assertEqual(payload['sources'][0]['repository_name'], 'igem')

    def test_errors(self):

        # The repository_id does not start with the part_name, even if url is valid
        source = IGEMSource(
            id=0, repository_name='igem', repository_id='BBa_C0062-dummy', sequence_file_url=self.good_url
        )

        # The url is not a GenBank file
        source_json = source.model_dump()
        source_json['sequence_file_url'] = self.no_gb_url
        response = client.post('/repository_id/igem', json=source_json)
        self.assertEqual(response.status_code, 422)

        # The url does not exist
        source = IGEMSource(id=0, repository_name='igem', repository_id='dummy-test', sequence_file_url=self.wrong_url)
        response = client.post('/repository_id/igem', json=source.model_dump())
        self.assertEqual(response.status_code, 404)


class GenomeRegionTest(unittest.TestCase):
    def assertStatusCode(self, response_code: int, expected: int, msg: str = ''):
        if response_code == 503:
            self.skipTest('NCBI not available')
        else:
            self.assertEqual(response_code, expected, msg)

    @pytest.mark.flaky(reruns=3, reruns_delay=2)
    def test_examples(self):
        for example_name in request_examples.genome_region_examples:
            example = request_examples.genome_region_examples[example_name]
            response = client.post('/genome_coordinates', json=example['value'])
            msg = f'Error in example {example_name}'
            self.assertStatusCode(response.status_code, 200, msg)
            payload = response.json()
            try:
                sources = [GenomeCoordinatesSource.model_validate(s) for s in payload['sources']]
            except Exception:
                self.fail(f'Cannot parse the sources for example {example_name}')
            response_source = sources[0]
            request_source = GenomeCoordinatesSource.model_validate(example['value'])
            if example_name == 'full' or example_name == 'full_with_genbank_accession':
                self.assertEqual(response_source, request_source, msg)
            elif example_name == 'id_omitted':
                request_source.gene_id = 2543372
                self.assertEqual(response_source, request_source, msg)
            elif example_name == 'assembly_accession_omitted':
                self.assertEqual(response_source, request_source, msg)
            elif example_name == 'viral_sequence':
                self.assertEqual(response_source, request_source, msg)

    @pytest.mark.flaky(reruns=3, reruns_delay=2)
    def test_exceptions(self):
        # Load first example
        correct_source = GenomeCoordinatesSource.model_validate(
            request_examples.genome_region_examples['full']['value']
        )

        # Ommit assembly accession
        s = correct_source.model_copy()
        s.assembly_accession = None
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 422)

        # Ommit locus_tag keeping gene id is now supported
        # Before it was not, but see https://github.com/ncbi/datasets/issues/397
        s = correct_source.model_copy()
        s.locus_tag = None
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 200)

        # Wrong gene_id (not matching that of the locus_tag)
        s = correct_source.model_copy()
        s.gene_id = 123
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 400)

        # Wrong assembly accession
        s = correct_source.model_copy()
        s.assembly_accession = 'blah'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 404)

        # Wrong locus_tag
        s = correct_source.model_copy()
        s.locus_tag = 'blah'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 404)

        # Wrong coordinates
        s = correct_source.model_copy()
        s.start = 1
        s.end = 10
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 400)

        # Wrong assembly accession
        s = correct_source.model_copy()
        s.locus_tag = None
        s.gene_id = None
        s.assembly_accession = 'blah'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 404)
        self.assertIn('Wrong assembly accession', response.json()['detail'])

        # Assembly accession not linked to any sequence record
        s = correct_source.model_copy()
        s.locus_tag = None
        s.gene_id = None
        # It used to be the case, but it changed so now we just mock the response
        s.assembly_accession = 'GCF_000146045.1'

        with respx.mock:
            respx.get(
                'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000146045.1/sequence_reports'
            ).mock(return_value=httpx.Response(200, json={'total_count': 0}))

            response = client.post('/genome_coordinates', json=s.model_dump())
            self.assertStatusCode(response.status_code, 400)
            self.assertIn('No sequence accessions linked', response.json()['detail'])

        # Assembly accession not linked to that sequence accession
        s = correct_source.model_copy()
        s.locus_tag = None
        s.gene_id = None
        s.assembly_accession = 'GCF_000146045.2'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 400)
        self.assertIn('not contained in assembly accession', response.json()['detail'])

        # Wrong sequence accession
        s = correct_source.model_copy()
        s.locus_tag = None
        s.gene_id = None
        s.assembly_accession = None
        s.sequence_accession = 'blah'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 404)

        # Coordinates malformatted
        viral_source = GenomeCoordinatesSource.model_validate(
            request_examples.genome_region_examples['viral_sequence']['value']
        )
        viral_source.start = 10
        viral_source.end = 1
        response = client.post('/genome_coordinates', json=viral_source.model_dump())
        self.assertStatusCode(response.status_code, 422)

        viral_source.start = 0
        viral_source.end = 20
        response = client.post('/genome_coordinates', json=viral_source.model_dump())
        self.assertStatusCode(response.status_code, 422)

        viral_source.start = 1
        viral_source.end = 20
        viral_source.strand = 0
        response = client.post('/genome_coordinates', json=viral_source.model_dump())
        self.assertStatusCode(response.status_code, 422)

        # Coordinates outside of the sequence
        viral_source.start = 1
        # the length is 2151
        viral_source.end = 2152
        viral_source.strand = 1
        response = client.post('/genome_coordinates', json=viral_source.model_dump())
        self.assertStatusCode(response.status_code, 400)

        # Coordinates too long
        viral_source.start = 1
        viral_source.end = 100004
        response = client.post('/genome_coordinates', json=viral_source.model_dump())
        self.assertStatusCode(response.status_code, 400)


class SEVASourceTest(unittest.TestCase):
    def test_seva_url(self):
        source = SEVASource(
            id=0,
            repository_id='pSEVA261',
            repository_name='seva',
            sequence_file_url='https://seva-plasmids.com/maps-canonical/maps-plasmids-SEVAs-canonical-versions-web-1-3-gbk/pSEVA261.gbk',
        )
        response = client.post('/repository_id/seva', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)
        out_source = payload['sources'][0]
        self.assertEqual(out_source, source.model_dump())
        seq = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(seq.name, 'pSEVA261')

    def test_ncbi_url(self):
        source = SEVASource(
            id=0,
            repository_id='pSEVA2214',
            repository_name='seva',
            sequence_file_url='https://www.ncbi.nlm.nih.gov/nuccore/MH650998',
        )
        response = client.post('/repository_id/seva', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)
        out_source = payload['sources'][0]
        self.assertEqual(out_source, source.model_dump())
        seq = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(seq.name, 'pSEVA2214')

    def test_errors(self):
        source = SEVASource(
            id=0,
            repository_id='pSEVA261',
            repository_name='seva',
            sequence_file_url='https://seva-plasmids.com/dummy.gbk',
        )

        response = client.post('/repository_id/seva', json=source.model_dump())
        self.assertEqual(response.status_code, 404)

        source_dict = source.model_dump()
        source_dict['repository_id'] = 'hello'
        response = client.post('/repository_id/seva', json=source_dict)
        self.assertEqual(response.status_code, 422)

        source_dict = source.model_dump()
        source_dict['sequence_file_url'] = 'hello'
        response = client.post('/repository_id/seva', json=source_dict)
        self.assertEqual(response.status_code, 422)

        source_dict = source.model_dump()
        source_dict['sequence_file_url'] = 'https://hello.com'
        response = client.post('/repository_id/seva', json=source_dict)
        self.assertEqual(response.status_code, 404)
        payload = response.json()
        self.assertIn('invalid SEVA url', payload['detail'])

        # Mock connection error
        with respx.mock:
            respx.get('https://seva-plasmids.com/dummy.gbk').mock(side_effect=httpx.ConnectError)
            response = client.post('/repository_id/seva', json=source.model_dump())
            self.assertEqual(response.status_code, 504)
            payload = response.json()
            self.assertEqual(payload['detail'], 'unable to connect to SEVA')

        # Mock incorrect file
        with respx.mock:
            respx.get('https://seva-plasmids.com/dummy.gbk').mock(return_value=httpx.Response(200, text='dummy test'))
            response = client.post('/repository_id/seva', json=source.model_dump())
            self.assertEqual(response.status_code, 400)
            payload = response.json()
            self.assertIn('Error parsing file', payload['detail'])
            self.assertIn('No sequences found in SEVA file', payload['detail'])

        with respx.mock:
            with open(f'{test_files}/ase1_body_error.gb', 'r') as f:
                mock_seq = f.read()
            respx.get('https://seva-plasmids.com/dummy.gbk').mock(return_value=httpx.Response(200, text=mock_seq))
            response = client.post('/repository_id/seva', json=source.model_dump())
            self.assertEqual(response.status_code, 400)
            payload = response.json()
            self.assertIn('Error parsing file', payload['detail'])
            self.assertIn('Premature end of line', payload['detail'])

    def test_redirect(self):
        source = SEVASource(
            id=0,
            repository_id='pSEVA261',
            repository_name='seva',
            sequence_file_url='https://seva-plasmids.com/maps-canonical/maps-plasmids-SEVAs-canonical-versions-web-1-3-gbk/pSEVA261.gbk',
        )
        with open(f'{test_files}/ase1.gb', 'r') as f:
            mock_seq = f.read()
        # We mock to avoid extra requests
        with respx.mock:
            respx.get(source.sequence_file_url).mock(return_value=httpx.Response(200, text=mock_seq))
            response = client.post('/repository_id', json=source.model_dump())
            self.assertEqual(response.status_code, 200)

    def test_circularize(self):
        source = SEVASource(
            id=0,
            repository_id='pSEVA261',
            repository_name='seva',
            sequence_file_url='https://seva-plasmids.com/maps-canonical/maps-plasmids-SEVAs-canonical-versions-web-1-3-gbk/pSEVA261.gbk',
        )
        with open(f'{test_files}/ase1.gb', 'r') as f:
            mock_seq = f.read()
        # We mock to avoid extra requests
        with respx.mock:
            respx.get(source.sequence_file_url).mock(return_value=httpx.Response(200, text=mock_seq))
            response = client.post('/repository_id', json=source.model_dump())
            self.assertEqual(response.status_code, 200)
            payload = response.json()
            seq = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
            self.assertEqual(seq.circular, True)
