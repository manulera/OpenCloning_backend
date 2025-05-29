import os
import unittest
from fastapi.testclient import TestClient
import json
from pydna.dseqrecord import Dseqrecord
from pydna.parsers import parse
import opencloning.main as _main
from opencloning.dna_functions import format_sequence_genbank, read_dsrecord_from_json
from opencloning.pydantic_models import TextFileSequence, BaseCloningStrategy
from opencloning_linkml._version import __version__ as schema_version


test_files = os.path.join(os.path.dirname(__file__), 'test_files')

client = TestClient(_main.app)


class VersionTestWithFiles(unittest.TestCase):

    def setUp(self):
        with open('./version.txt', 'w') as f:
            f.write('1.2.3')
        with open('./commit_sha.txt', 'w') as f:
            f.write('1234567890')

    def tearDown(self):
        os.remove('./version.txt')
        os.remove('./commit_sha.txt')

    def test_version_file(self):
        response = client.get('/version')
        self.assertEqual(response.status_code, 200)
        resp = response.json()
        self.assertEqual(resp['version'], '1.2.3')
        self.assertEqual(resp['commit_sha'], '1234567890')


class VersionTestWithoutFiles(unittest.TestCase):
    def test_version_empty(self):
        response = client.get('/version')
        self.assertEqual(response.status_code, 200)
        resp = response.json()
        self.assertIsNone(resp['version'])
        self.assertIsNone(resp['commit_sha'])


class ValidatorTest(unittest.TestCase):
    def test_validator(self):
        with open(f'{test_files}/homologous_recombination.json') as ins:
            # Read it to json
            data = json.load(ins)
        BaseCloningStrategy.model_validate(data)


class ValidateEndPointTest(unittest.TestCase):

    def test_validate(self):
        # Valid file
        with open(f'{test_files}/homologous_recombination.json') as ins:
            # Read it to json
            data = json.load(ins)
        response = client.post('/validate', json=data)
        self.assertEqual(response.status_code, 200)
        self.assertNotIn('x-warning', response.headers)

        # Invalid file
        data['dummy'] = 'dummy'
        response = client.post('/validate', json=data)
        self.assertEqual(response.status_code, 422)

        # Completely wrong file
        data = {'dummy': 'dummy'}
        response = client.post('/validate', json=data)
        self.assertEqual(response.status_code, 422)

        # Old format file
        with open(f'{test_files}/homologous_recombination_old_format.json') as ins:
            # Read it to json
            data = json.load(ins)
        response = client.post('/validate', json=data)
        self.assertEqual(response.status_code, 200)
        self.assertIn('x-warning', response.headers)
        self.assertNotIn('contained an error', response.headers['x-warning'])
        self.assertIn(
            'The cloning strategy is in a previous version of the model and has been migrated to the latest version.',
            response.headers['x-warning'],
        )
        # The data has been migrated to the latest version
        cs = BaseCloningStrategy.model_validate(response.json())
        self.assertEqual(cs.schema_version, schema_version)

        # File with errors
        with open(f'{test_files}/bug_fixing/digestion_spanning_origin.json') as ins:
            # Read it to json
            data = json.load(ins)
        response = client.post('/validate', json=data)
        self.assertEqual(response.status_code, 200)
        self.assertIn('x-warning', response.headers)
        self.assertIn('previous version', response.headers['x-warning'])
        self.assertIn('contained an error', response.headers['x-warning'])

        # The right source has been turned into a template
        data = response.json()
        seq = next(s for s in data['sequences'] if s['id'] == 6)
        self.assertEqual(seq['type'], 'TemplateSequence')


class RenameSequenceTest(unittest.TestCase):

    def test_rename(self):
        dseqr = Dseqrecord('ACGT')
        dseqr.name = 'original'
        json_seq = format_sequence_genbank(dseqr)
        json_seq.id = 0
        response = client.post('/rename_sequence?name=hello', json=json_seq.model_dump())
        self.assertEqual(response.status_code, 200)

        payload = response.json()
        dseq_resp = read_dsrecord_from_json(TextFileSequence.model_validate(payload))
        self.assertEqual(dseq_resp.name, 'hello')

    def test_error(self):
        """Does not allow spaces"""
        dseqr = Dseqrecord('ACGT')
        dseqr.name = 'original'
        json_seq = format_sequence_genbank(dseqr)
        json_seq.id = 0
        response = client.post('/rename_sequence?name=hello world', json=json_seq.model_dump())
        self.assertEqual(response.status_code, 422)


class RestrictionEnzymeListTest(unittest.TestCase):

    def test_restriction_enzyme_list(self):
        response = client.get('/restriction_enzyme_list')
        assert response.status_code == 200
        assert 'EcoRI' in response.json()['enzyme_names']


class AlignSangerTest(unittest.TestCase):

    def test_align_sanger(self):
        seq = parse(os.path.join(test_files, 'GIN11M86.gb'))[0].looped()
        json_seq = format_sequence_genbank(seq)
        json_seq.id = 0
        trace = 'ttgcagcattttgtctttctataaaaatgtgtcgttcctttttttcattttttggcgcgtcgcctcggggtcgtatagaatatg'
        response = client.post('/align_sanger', json={'sequence': json_seq.model_dump(), 'traces': [trace]})
        assert response.status_code == 200
        assert len(response.json()) == 2

    def test_errors(self):
        seq = Dseqrecord('ACGT', circular=True)
        json_seq = format_sequence_genbank(seq)
        json_seq.id = 0
        response = client.post('/align_sanger', json={'sequence': json_seq.model_dump(), 'traces': ['ACGT']})
        assert response.status_code == 400
