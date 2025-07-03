from fastapi.testclient import TestClient
import unittest
import os

from opencloning.dna_functions import read_dsrecord_from_json
import opencloning.main as _main
from opencloning.pydantic_models import (
    TextFileSequence,
)


test_files = os.path.join(os.path.dirname(__file__), 'test_files')

client = TestClient(_main.app)


class ZiqiangEtAl2024Test(unittest.TestCase):
    def test_ziqiang_et_al_2024(self):
        # Default call should produce the desired product
        protospacers = [
            'GCTGGCTAACCGTGAGGGGA',
            'CCGTGTACTGTAGTTACAGT',
            'TGTGGTTCCCCGGCCGTCTT',
            'ATACTCTAGTCCTCAACGCC',
        ]
        response = client.post('/batch_cloning/ziqiang_et_al2024', json=protospacers)
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['primers']), 14)
        # Find the LR source
        lr_source = next(s for s in payload['sources'] if s['type'] == 'GatewaySource' and s['reaction_type'] == 'LR')
        self.assertIsNotNone(lr_source)
        # Get the output sequence
        seq_id = lr_source['id']
        seq = next(s for s in payload['sequences'] if s['id'] == seq_id)
        self.assertIsNotNone(seq)
        dseq = read_dsrecord_from_json(TextFileSequence.model_validate(seq))
        self.assertEqual(dseq.name, 'expression_clone')
        self.assertEqual(dseq.seguid(), 'cdseguid=kzzhO-94Zw1KwHwS1vRcggoXxqU')

        # If we stop after BP, we should get the same sequence
        response = client.post('/batch_cloning/ziqiang_et_al2024', json=protospacers, params={'until_bp': True})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['primers']), 14)
        # No LR should exist
        lr_source = next(
            (s for s in payload['sources'] if s['type'] == 'GatewaySource' and s['reaction_type'] == 'LR'), None
        )
        self.assertIsNone(lr_source)
        bp_source = next(
            (s for s in payload['sources'] if s['type'] == 'GatewaySource' and s['reaction_type'] == 'BP'), None
        )
        self.assertIsNotNone(bp_source)
        seq = next(s for s in payload['sequences'] if s['id'] == bp_source['id'])
        self.assertIsNotNone(seq)
        dseq = read_dsrecord_from_json(TextFileSequence.model_validate(seq))
        self.assertEqual(dseq.name, 'entry_clone')
        self.assertEqual(dseq.seguid(), 'cdseguid=UC2hIJs2Ba3M6he3oUOoLf-I1to')

    def test_protospacer_validation(self):

        # Test that the protospacers are valid
        response = client.post('/batch_cloning/ziqiang_et_al2024', json=[])
        self.assertEqual(response.status_code, 422)
        response = client.post('/batch_cloning/ziqiang_et_al2024', json=[''])
        self.assertEqual(response.status_code, 400)
        response = client.post('/batch_cloning/ziqiang_et_al2024', json=['A'])
        self.assertEqual(response.status_code, 400)
        response = client.post('/batch_cloning/ziqiang_et_al2024', json=['A' * 21])
        self.assertEqual(response.status_code, 400)
        response = client.post('/batch_cloning/ziqiang_et_al2024', json=['b' * 20])
        self.assertEqual(response.status_code, 400)
        response = client.post('/batch_cloning/ziqiang_et_al2024', json=['A' * 8 + 'AACA' + 'A' * 8])
        self.assertEqual(response.status_code, 400)
        response = client.post('/batch_cloning/ziqiang_et_al2024', json=['A' * 8 + 'GCTT' + 'A' * 8])
        self.assertEqual(response.status_code, 400)
        response = client.post(
            '/batch_cloning/ziqiang_et_al2024', json=['A' * 8 + 'ACCA' + 'A' * 8, 'T' * 8 + 'ACCA' + 'T' * 8]
        )
        self.assertEqual(response.status_code, 400)
