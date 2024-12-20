from fastapi.testclient import TestClient
import unittest
import shutil
import os
from shareyourcloning.pydantic_models import ManuallyTypedSource, RestrictionEnzymeDigestionSource
from pytest import MonkeyPatch
from importlib import reload
from pydna.dseqrecord import Dseqrecord

from shareyourcloning.dna_functions import format_sequence_genbank
import shareyourcloning.get_router as get_router
import shareyourcloning.endpoints.no_input as no_input_endpoints
import shareyourcloning.endpoints.no_assembly as no_assembly_endpoints
import shareyourcloning.main as main
import shareyourcloning.app_settings as app_settings


class StubRouteTest(unittest.TestCase):
    # DO this before each test
    def setUp(self):
        # Has to be imported here to get the right environment variable
        MonkeyPatch().setenv('RECORD_STUBS', '1')
        reload(app_settings)
        reload(get_router)
        reload(no_input_endpoints)
        reload(no_assembly_endpoints)
        reload(main)

        client = TestClient(main.app)
        self.client = client
        # remove the stubs folder
        shutil.rmtree('stubs', ignore_errors=True)

    # Remove the stubs folder after each test
    def tearDown(self):
        shutil.rmtree('stubs', ignore_errors=True)
        MonkeyPatch().setenv('RECORD_STUBS', '0')

        reload(app_settings)
        reload(get_router)
        reload(no_input_endpoints)
        reload(no_assembly_endpoints)
        reload(main)

    def test_stub_route(self):
        source = ManuallyTypedSource(
            id=0,
            user_input='ATGC',
        )

        response = self.client.post('/manually_typed', json=source.model_dump())
        self.assertEqual(response.status_code, 200)

        self.assertTrue(os.path.exists('stubs/manually_typed/'))
        self.assertTrue(len(os.listdir('stubs/manually_typed/')), 1)

        # Also works for 422 response
        source_dict = source.model_dump()
        source_dict['user_input'] = 'io'
        response = self.client.post('/manually_typed', json=source_dict)
        self.assertEqual(response.status_code, 422)
        self.assertEqual(len(os.listdir('stubs/manually_typed/')), 2)

        # Works for 404
        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', circular=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # One enzyme
        source = RestrictionEnzymeDigestionSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = self.client.post('/restriction', json=data, params={'restriction_enzymes': ['helloworld']})
        self.assertEqual(response.status_code, 404)
        self.assertEqual(len(os.listdir('stubs/restriction/')), 1)


if __name__ == '__main__':
    unittest.main()
