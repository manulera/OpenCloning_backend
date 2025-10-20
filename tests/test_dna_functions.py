from urllib.error import HTTPError
import unittest
import os
import respx
import httpx
from pydna.dseq import Dseq
from Bio.Seq import reverse_complement

from opencloning.dna_functions import (
    custom_file_parser,
    correct_name,
    MyGenBankScanner,
    get_sequence_from_euroscarf_url,
    oligonucleotide_hybridization_overhangs,
)

test_files = os.path.join(os.path.dirname(__file__), 'test_files')


class PermisiveParserWithApeTest(unittest.TestCase):
    def test_permisive_parser_with_ape_circular(self):
        with open(f'{test_files}/P2RP3.ape', 'r') as f:
            plasmid = custom_file_parser(f, 'genbank')[0]
            # Since APE files are not correctly gb formatted (as of 2024-11-27)
            # the Bio.SeqIO.parse may not recognize the topology of the plasmid
            # Our custom permissive parser should be then used and the topology
            # parameter properly recognized
            self.assertEqual(plasmid.circular, True)

    def test_permisive_parser_with_ape_linear(self):
        with open(f'{test_files}/P2RP3_linear.ape', 'r') as f:
            # I manually changed the topology of the plasmid to linear
            plasmid = custom_file_parser(f, 'genbank')[0]
            self.assertEqual(plasmid.circular, False)

    def test_permisive_parser_no_topology(self):
        with open(f'{test_files}/ase1_no_topology.gb', 'r') as f:
            plasmid = custom_file_parser(f, 'genbank')[0]
            self.assertEqual(plasmid.circular, False)

    def test_custom_file_parser_body_error(self):
        with open(f'{test_files}/ase1_body_error.gb', 'r') as f:
            with self.assertRaises(ValueError):
                custom_file_parser(f, 'genbank')


class PermissiveParserOtherTest(unittest.TestCase):
    def test_permissive_parser_other(self):
        with open(f'{test_files}/pSEVA427.gbk', 'r') as f:
            plasmid = custom_file_parser(f, 'genbank')[0]
            self.assertEqual(plasmid.circular, True)


class MinorFunctionsTest(unittest.TestCase):
    def test_correct_name(self):
        file = f'{test_files}/addgene-plasmid-39296-sequence-49545.gbk'
        with open(file, 'r') as f:
            dseq = custom_file_parser(f, 'genbank')[0]
        correct_name(dseq)
        self.assertEqual(dseq.name, 'pFA6a-kanMX6')

    def test_error_on_genbank_scanner(self):

        with self.assertRaises(ValueError):
            MyGenBankScanner(debug=0)._feed_first_line(None, 'LOCUS hello bye')


class MinorFunctionsAsyncTest(unittest.IsolatedAsyncioTestCase):
    @respx.mock
    async def test_error_euroscarf(self):

        # Connection error
        respx.get('http://www.euroscarf.de/plasmid_details.php').mock(
            side_effect=httpx.ConnectError('Connection error')
        )
        with self.assertRaises(HTTPError) as e:
            await get_sequence_from_euroscarf_url('blah')
        self.assertEqual(e.exception.code, 504)
        self.assertIn('could not connect to euroscarf', str(e.exception))

        # As far as I can tell, this never happens (it always returns a 200 even if the page is missing)
        respx.get('http://www.euroscarf.de/plasmid_details.php').respond(503, text='')
        with self.assertRaises(HTTPError) as e:
            await get_sequence_from_euroscarf_url('blah')
        self.assertEqual(e.exception.code, 503)
        self.assertIn('could not connect to euroscarf', str(e.exception))

        # If the format of the page would change, these errors should be raised
        respx.get('http://www.euroscarf.de/plasmid_details.php').respond(200, text='')
        with self.assertRaises(HTTPError) as e:
            await get_sequence_from_euroscarf_url('blah')
        self.assertEqual(e.exception.code, 503)
        self.assertIn('Could not retrieve plasmid details', str(e.exception))

        respx.get('http://www.euroscarf.de/plasmid_details.php').respond(200, text='<body>missing other</body>')
        with self.assertRaises(HTTPError) as e:
            await get_sequence_from_euroscarf_url('blah')
        self.assertEqual(e.exception.code, 503)
        self.assertIn('Could not retrieve plasmid details', str(e.exception))


class OligonucleotideHybridizationTest(unittest.TestCase):

    def test_oligonucleotide_hybridization_overhangs(self):
        self.assertEqual(
            [-4], oligonucleotide_hybridization_overhangs('CTCGatcggtgtgaaaagtcagtatccagtcgtgtag', 'tttcacaccgat', 12)
        )
        self.assertEqual(
            [21], oligonucleotide_hybridization_overhangs('tttcacaccgat', 'CTCGatcggtgtgaaaagtcagtatccagtcgtgtag', 12)
        )

        for ovhgs in [[3, 0], [0, 2], [3, 2], [-3, -2], [3, -2]]:
            dseq = Dseq.from_full_sequence_and_overhangs('GGACAATATATGGCAC', *ovhgs)
            self.assertEqual([ovhgs[0]], oligonucleotide_hybridization_overhangs(dseq.watson, dseq.crick, 10))
            self.assertEqual([ovhgs[1]], oligonucleotide_hybridization_overhangs(dseq.crick, dseq.watson, 10))

        seq1 = 'GGACAATATATGGCAC'
        seq2 = 'a' + reverse_complement(seq1)
        seq1 += 'c'

        self.assertRaises(ValueError, oligonucleotide_hybridization_overhangs, seq1, seq2, 10)
        self.assertRaises(
            ValueError, oligonucleotide_hybridization_overhangs, reverse_complement(seq1), reverse_complement(seq2), 10
        )
