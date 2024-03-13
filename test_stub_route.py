from fastapi.testclient import TestClient
import unittest
import shutil
import os
from pydantic_models import ManuallyTypedSource, RestrictionEnzymeDigestionSource
from pytest import MonkeyPatch
from importlib import reload
from dna_functions import format_sequence_genbank
from pydna.dseqrecord import Dseqrecord


class StubRouteTest(unittest.TestCase):
    # DO this before each test
    def setUp(self):
        # Has to be imported here to get the right environment variable
        MonkeyPatch().setenv('RECORD_STUBS', '1')
        import main

        reload(main)
        client = TestClient(main.app)
        self.client = client
        # remove the stubs folder
        shutil.rmtree('stubs', ignore_errors=True)

    # Remove the stubs folder after each test
    def tearDown(self):
        shutil.rmtree('stubs', ignore_errors=True)

    def test_stub_route(self):
        source = ManuallyTypedSource(
            user_input='ATGC',
        )

        response = self.client.post('/manually_typed', json=source.model_dump())
        self.assertEqual(response.status_code, 200)

        self.assertTrue(os.path.exists('stubs/manually_typed/'))
        self.assertTrue(len(os.listdir('stubs/manually_typed/')), 1)

        # Also works for 422 response
        source.user_input = 'io'
        response = self.client.post('/manually_typed', json=source.model_dump())
        self.assertEqual(response.status_code, 422)
        self.assertEqual(len(os.listdir('stubs/manually_typed/')), 2)

        # Works for 404
        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', circular=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # One enzyme
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['helloworld'],
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = self.client.post('/restriction', json=data)
        self.assertEqual(response.status_code, 404)
        self.assertEqual(len(os.listdir('stubs/restriction/')), 1)


if __name__ == '__main__':
    unittest.main()
