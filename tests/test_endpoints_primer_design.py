from fastapi.testclient import TestClient
from pydna.dseqrecord import Dseqrecord
import unittest
import copy
from Bio.Seq import reverse_complement
import os
from primer3 import bindings
from itertools import product

from opencloning.dna_functions import format_sequence_genbank
import opencloning.main as _main
from opencloning.pydantic_models import (
    PrimerModel,
    SimpleSequenceLocation as PydanticSimpleLocation,
)
from opencloning.endpoints.primer_design import (
    PrimerDetailsResponse,
    ThermodynamicResult,
)

test_files = os.path.join(os.path.dirname(__file__), 'test_files')

client = TestClient(_main.app)


class PrimerDesignTest(unittest.TestCase):

    def test_homologous_recombination(self):
        pcr_seq = format_sequence_genbank(Dseqrecord('AATGATGGATGACATTCAAAGCACTGATTCTATTGCTGAAAAAGATAAT'))
        pcr_seq.id = 1
        pcr_loc = PydanticSimpleLocation(start=4, end=44)
        hr_seq = format_sequence_genbank(Dseqrecord('AAACGTTT'))
        hr_seq.id = 2
        hr_loc_replace = PydanticSimpleLocation(start=3, end=5)

        homology_length = 3
        minimal_hybridization_length = 10
        insert_forward = True
        target_tm = 30

        # First we replace the CG
        data = {
            'pcr_template': {
                'sequence': pcr_seq.model_dump(),
                'location': pcr_loc.model_dump(),
            },
            'homologous_recombination_target': {
                'sequence': hr_seq.model_dump(),
                'location': hr_loc_replace.model_dump(),
            },
        }
        params = {
            'homology_length': homology_length,
            'minimal_hybridization_length': minimal_hybridization_length,
            'insert_forward': insert_forward,
            'target_tm': target_tm,
        }
        response = client.post('/primer_design/homologous_recombination', json=data, params=params)
        payload = response.json()

        self.assertEqual(response.status_code, 200)
        self.assertEqual(payload['primers'][0]['sequence'], 'aaaATGGATGACATT')
        self.assertEqual(payload['primers'][1]['sequence'], 'aaaCTTTTTCAGCAA')

        # Raise valuerror
        params['homology_length'] = 10
        response = client.post('/primer_design/homologous_recombination', json=data, params=params)
        self.assertEqual(response.status_code, 400)
        payload = response.json()
        self.assertEqual(payload['detail'], 'Forward homology region is out of bounds.')

        # Test an insertion with spacers and reversed insert
        params['homology_length'] = 3
        data['pcr_template']['forward_orientation'] = False
        data['homologous_recombination_target']['location'] = PydanticSimpleLocation(start=3, end=3).model_dump()
        data['spacers'] = ['attt', 'cggg']
        response = client.post('/primer_design/homologous_recombination', json=data, params=params)
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload['primers'][0]['sequence'], 'aaaatttCTTTTTCAGCAA')
        self.assertEqual(payload['primers'][1]['sequence'], 'acgcccgATGGATGACATT')

        # Raise error if the number of spacers is incorrect
        data['spacers'] = ['attt', 'cggg', 'tttt']
        response = client.post('/primer_design/homologous_recombination', json=data, params=params)
        self.assertEqual(response.status_code, 422)
        self.assertIn('The number of spacers must be', response.json()['detail'])

        # Raise error if the spacer is not DNA
        data['spacers'] = ['zzz', 'cggg']
        response = client.post('/primer_design/homologous_recombination', json=data, params=params)
        self.assertEqual(response.status_code, 422)
        self.assertIn('Spacer can only contain ACGT bases', response.json()['detail'])

    def test_gibson_assembly(self):
        # Test case for gibson_assembly_primers endpoint
        templates = [
            Dseqrecord('AAACAGTAATACGTTCCTTTTTTATGATGATGGATGACATTCAAAGCACTGATTCTAT'),
            Dseqrecord('GTTTACAACGGCAATGAACGTTCCTTTTTTATGATATGCCCAGCTTCATGAAATGGAA'),
            Dseqrecord('AAGGACAACGTTCCTTTTTTATGATATATATGGCACAGTATGATCAAAAGTTAAGTAC'),
        ]

        queries = []
        for i, template in enumerate(templates):
            json_seq = format_sequence_genbank(template)
            json_seq.id = i
            queries.append(
                {
                    'sequence': json_seq.model_dump(),
                    'location': PydanticSimpleLocation(start=0, end=len(template)).model_dump(),
                    'forward_orientation': True,
                }
            )

        params = {'homology_length': 20, 'minimal_hybridization_length': 15, 'target_tm': 55, 'circular': True}

        response = client.post(
            '/primer_design/gibson_assembly', json={'pcr_templates': queries, 'spacers': None}, params=params
        )

        self.assertEqual(response.status_code, 200)

        payload = response.json()
        self.assertEqual(len(payload['primers']), 6)  # 2 primers per template
        for i, p in enumerate(payload['primers']):
            p = PrimerModel.model_validate(p)
            # Check that the name is correct
            if i % 2 == 0:
                self.assertEqual(p.name, f'seq_{i//2}_fwd')
            else:
                self.assertEqual(p.name, f'seq_{i//2}_rvs')

        # Primer naming also work for named sequences
        for i, t in enumerate(templates):
            t.name = f'template_{i}'

        queries = []
        for i, template in enumerate(templates):
            json_seq = format_sequence_genbank(template)
            json_seq.id = i
            queries.append(
                {
                    'sequence': json_seq.model_dump(),
                    'location': PydanticSimpleLocation(start=0, end=len(template)).model_dump(),
                    'forward_orientation': True,
                }
            )

        params = {'homology_length': 20, 'minimal_hybridization_length': 15, 'target_tm': 55, 'circular': True}

        response = client.post(
            '/primer_design/gibson_assembly', json={'pcr_templates': queries, 'spacers': None}, params=params
        )

        self.assertEqual(response.status_code, 200)

        payload = response.json()
        self.assertEqual(len(payload['primers']), 6)  # 2 primers per template
        for i, p in enumerate(payload['primers']):
            p = PrimerModel.model_validate(p)
            # Check that the name is correct
            if i % 2 == 0:
                self.assertEqual(p.name, f'template_{i//2}_fwd')
            else:
                self.assertEqual(p.name, f'template_{i//2}_rvs')

        # Test error case with invalid parameters
        params['minimal_hybridization_length'] = 100  # Too long
        response = client.post(
            '/primer_design/gibson_assembly', json={'pcr_templates': queries, 'spacers': None}, params=params
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn('Primers could not be designed', response.json()['detail'])

        # Test case with spacers
        params['homology_length'] = 20
        params['circular'] = True
        params['minimal_hybridization_length'] = 15
        spacers = ['aaaa', 'tttt', 'cccc']
        response = client.post(
            '/primer_design/gibson_assembly', json={'pcr_templates': queries, 'spacers': spacers}, params=params
        )
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['primers']), 6)  # 2 primers per template
        primers = [PrimerModel.model_validate(p) for p in payload['primers']]
        self.assertTrue(primers[0].sequence.startswith('TTAAGTACccccAAACAGTA'))
        self.assertTrue(reverse_complement(primers[-1].sequence).endswith('TTAAGTACccccAAACAGTA'))
        self.assertTrue(reverse_complement(primers[1].sequence).endswith('GATTCTATaaaaGTTTACAA'))
        self.assertTrue(primers[2].sequence.startswith('GATTCTATaaaaGTTTACAA'))
        self.assertTrue(reverse_complement(primers[3].sequence).endswith('AAATGGAAttttAAGGACAA'))
        self.assertTrue(primers[4].sequence.startswith('AAATGGAAttttAAGGACAA'))

        # Test that wrong number of spacers fails
        response = client.post(
            '/primer_design/gibson_assembly',
            json={'pcr_templates': queries, 'spacers': ['aaaa', 'tttt']},
            params=params,
        )
        self.assertEqual(response.status_code, 422)

        # Test that non-DNA spacers fails
        response = client.post(
            '/primer_design/gibson_assembly',
            json={'pcr_templates': queries, 'spacers': ['hello', 'TTTT']},
            params=params,
        )
        self.assertEqual(response.status_code, 422)

    def test_simple_pair(self):
        from Bio.Restriction import EcoRI, BamHI

        # Create a test sequence
        dseqr = Dseqrecord('aaaGGCTTCACCAAGTCCTTGGAACAGccc')
        dseqr.name = 'test_sequence'
        json_seq = format_sequence_genbank(dseqr)
        json_seq.id = 0

        query = {
            'sequence': json_seq.model_dump(),
            'location': PydanticSimpleLocation(start=3, end=27).model_dump(),
            'forward_orientation': True,
        }

        params = {
            'minimal_hybridization_length': 10,
            'target_tm': 30,
            'left_enzyme': 'EcoRI',
            'right_enzyme': 'BamHI',
            'filler_bases': 'GC',
        }

        response = client.post('/primer_design/simple_pair', json={'pcr_template': query}, params=params)
        self.assertEqual(response.status_code, 200)

        payload = response.json()
        self.assertEqual(len(payload['primers']), 2)  # 2 primers (forward and reverse)

        fwd_primer = PrimerModel.model_validate(payload['primers'][0])
        rvs_primer = PrimerModel.model_validate(payload['primers'][1])

        self.assertEqual(fwd_primer.name, 'test_sequence_EcoRI_fwd')
        self.assertEqual(rvs_primer.name, 'test_sequence_BamHI_rvs')

        self.assertTrue(fwd_primer.sequence.startswith('GC' + str(EcoRI.site)))
        self.assertTrue(rvs_primer.sequence.startswith('GC' + str(BamHI.site)))

        # Same without enzymes
        params_no_enzymes = copy.deepcopy(params)
        params_no_enzymes.pop('left_enzyme')
        params_no_enzymes.pop('right_enzyme')
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query}, params=params_no_enzymes)
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['primers']), 2)  # 2 primers (forward and reverse)
        fwd_primer = PrimerModel.model_validate(payload['primers'][0])
        rvs_primer = PrimerModel.model_validate(payload['primers'][1])
        self.assertEqual(fwd_primer.name, 'test_sequence_fwd')
        self.assertEqual(rvs_primer.name, 'test_sequence_rvs')
        self.assertTrue(fwd_primer.sequence.startswith('GGCTT'))
        self.assertTrue(rvs_primer.sequence.startswith('CTGTT'))

        # Same, now inverted
        query2 = copy.deepcopy(query)
        query2['forward_orientation'] = False
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query2}, params=params_no_enzymes)
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['primers']), 2)  # 2 primers (forward and reverse)
        fwd_primer = PrimerModel.model_validate(payload['primers'][0])
        rvs_primer = PrimerModel.model_validate(payload['primers'][1])
        self.assertEqual(fwd_primer.name, 'test_sequence_fwd')
        self.assertEqual(rvs_primer.name, 'test_sequence_rvs')
        self.assertTrue(fwd_primer.sequence.startswith('CTGTT'))
        self.assertTrue(rvs_primer.sequence.startswith('GGCTT'))

        # Test primer name when sequence name is unset
        query2 = copy.deepcopy(query)
        dseqr2 = copy.deepcopy(dseqr)
        dseqr2.name = 'name'
        json_seq2 = format_sequence_genbank(dseqr2)
        json_seq2.id = 0
        query2['sequence'] = json_seq2.model_dump()
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query2}, params=params)
        payload = response.json()
        fwd_primer = PrimerModel.model_validate(payload['primers'][0])
        rvs_primer = PrimerModel.model_validate(payload['primers'][1])

        self.assertEqual(fwd_primer.name, 'seq_0_EcoRI_fwd')
        self.assertEqual(rvs_primer.name, 'seq_0_BamHI_rvs')

        # Test with spacers
        response = client.post(
            '/primer_design/simple_pair',
            json={'pcr_template': query, 'spacers': ['ATTT', 'CCAG']},
            params=params,
        )
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['primers']), 2)

        fwd_primer = PrimerModel.model_validate(payload['primers'][0])
        rvs_primer = PrimerModel.model_validate(payload['primers'][1])

        self.assertEqual(fwd_primer.name, 'test_sequence_EcoRI_fwd')
        self.assertEqual(rvs_primer.name, 'test_sequence_BamHI_rvs')

        self.assertTrue(fwd_primer.sequence.startswith('GC' + str(EcoRI.site) + 'ATTT'))
        self.assertTrue(rvs_primer.sequence.startswith('GC' + str(BamHI.site) + 'CTGG'))

        # Test error case with invalid parameters
        params['minimal_hybridization_length'] = 100  # Too long
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query}, params=params)
        self.assertEqual(response.status_code, 400)
        self.assertIn('Primers could not be designed', response.json()['detail'])

        # Test error case with invalid enzyme
        params['minimal_hybridization_length'] = 10  # Reset to valid value
        params['left_enzyme'] = 'InvalidEnzyme'
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query}, params=params)
        self.assertEqual(response.status_code, 404)
        self.assertIn('These enzymes do not exist', response.json()['detail'])

        # Test error case with wrong filler bases
        params['left_enzyme'] = 'EcoRI'
        params['filler_bases'] = 'zAA'
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query}, params=params)
        self.assertEqual(response.status_code, 400)
        self.assertIn('Filler bases can only contain ACGT bases.', response.json()['detail'])

        # Test error case with wrong spacer number
        params['filler_bases'] = 'GC'
        response = client.post(
            '/primer_design/simple_pair', json={'pcr_template': query, 'spacers': ['ATGC']}, params=params
        )
        self.assertEqual(response.status_code, 422)
        self.assertIn('The number of spacers must be', response.json()['detail'])

    def test_primer_details(self):
        # Works with short sequences
        SEQUENCE = 'ATGCATGCATGCATGC'
        response = client.get('/primer_details', params={'sequence': SEQUENCE})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resp = PrimerDetailsResponse.model_validate(payload)
        self.assertEqual(resp.gc_content, 0.5)
        self.assertEqual(resp.melting_temperature, bindings.calc_tm(SEQUENCE))
        homodimerResult = bindings.calc_homodimer(SEQUENCE, output_structure=True)
        self.assertEqual(resp.homodimer.melting_temperature, homodimerResult.tm)
        self.assertEqual(resp.homodimer.deltaG, homodimerResult.dg)
        self.assertEqual(resp.homodimer.figure, '\n'.join(homodimerResult.ascii_structure_lines))
        hairpinResult = bindings.calc_hairpin(SEQUENCE, output_structure=True)
        self.assertEqual(resp.hairpin.melting_temperature, hairpinResult.tm)
        self.assertEqual(resp.hairpin.deltaG, hairpinResult.dg)
        self.assertEqual(resp.hairpin.figure, '\n'.join(hairpinResult.ascii_structure_lines))

        # Splits long sequence
        LONG_SEQUENCE = 'GGAAAAGCATTTTTCTAAAATTGAAAGGCTTCACCAAGTCCTTGGAACAGATGGAGACAATTCATCATTA'
        response = client.get('/primer_details', params={'sequence': LONG_SEQUENCE})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resp = PrimerDetailsResponse.model_validate(payload)
        # Use the homodimer result with the lowest deltaG
        min_deltaG_homodimer = min(
            bindings.calc_homodimer(seq, output_structure=True).dg for seq in [LONG_SEQUENCE[:60], LONG_SEQUENCE[60:]]
        )
        self.assertEqual(resp.homodimer.deltaG, min_deltaG_homodimer)

        # Same with hairpin
        min_deltaG_hairpin = min(
            bindings.calc_hairpin(seq, output_structure=True).dg for seq in [LONG_SEQUENCE[:60], LONG_SEQUENCE[60:]]
        )
        self.assertEqual(resp.hairpin.deltaG, min_deltaG_hairpin)

        # Response is the same for upper and lower case
        response_lower = client.get('/primer_details', params={'sequence': LONG_SEQUENCE.lower()})
        self.assertEqual(response_lower.json()['melting_temperature'], response.json()['melting_temperature'])

        # Error 422 if sequence is not DNA
        response = client.get('/primer_details', params={'sequence': 'ATGCATGCATGCATGCX'})
        self.assertEqual(response.status_code, 422)

        # Handles the case where there is no hairpin / homodimer
        response = client.get('/primer_details', params={'sequence': 'AAAA'})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resp = PrimerDetailsResponse.model_validate(payload)
        self.assertIsNone(resp.homodimer)
        self.assertIsNone(resp.hairpin)

    def test_primer_heterodimer(self):
        # Works with short sequences
        SEQUENCE1 = 'ATGCATGCATGCATGC'
        SEQUENCE2 = 'CTAAAATTGAAAGGCTTCACCAAGT'
        response = client.get('/primer_heterodimer', params={'sequence1': SEQUENCE1, 'sequence2': SEQUENCE2})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resp = ThermodynamicResult.model_validate(payload)
        heterodimer = bindings.calc_heterodimer(SEQUENCE1, SEQUENCE2, output_structure=True)
        self.assertEqual(resp.melting_temperature, heterodimer.tm)
        self.assertEqual(resp.deltaG, heterodimer.dg)
        self.assertEqual(resp.figure, '\n'.join(heterodimer.ascii_structure_lines))

        # Same with lower case
        response_lower = client.get(
            '/primer_heterodimer', params={'sequence1': SEQUENCE1.lower(), 'sequence2': SEQUENCE2.lower()}
        )
        self.assertEqual(response_lower.json()['melting_temperature'], response.json()['melting_temperature'])

        # Error 422 if sequence is not DNA
        response = client.get('/primer_heterodimer', params={'sequence1': 'ATGCATGCATGCATGCX', 'sequence2': SEQUENCE2})
        self.assertEqual(response.status_code, 422)

        # Handles long sequences
        LONG_SEQUENCE1 = 'GGAAAAGCATTTTTCTAAAATTGAAAGGCTTCACCAAGTCCTTGGAACAGATGGAGACAATTCATCATTA'
        LONG_SEQUENCE2 = 'GGAAAAGCATTTTTCTAAAATTGAAAGGCTTCGACTATTCCTTGGAACAGATGGAGACAATTCATCATTA'
        response = client.get('/primer_heterodimer', params={'sequence1': LONG_SEQUENCE1, 'sequence2': LONG_SEQUENCE2})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resp = ThermodynamicResult.model_validate(payload)
        combs = list(product([LONG_SEQUENCE1[:60], LONG_SEQUENCE1[60:]], [LONG_SEQUENCE2[:60], LONG_SEQUENCE2[60:]]))
        min_deltaG = min(bindings.calc_heterodimer(seq1, seq2, output_structure=True).dg for seq1, seq2 in combs)
        self.assertEqual(resp.deltaG, min_deltaG)

        # Handles case where there is no heterodimer
        response = client.get('/primer_heterodimer', params={'sequence1': 'AAA', 'sequence2': 'AAA'})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertIsNone(payload)
