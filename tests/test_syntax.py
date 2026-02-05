import os
from opencloning.syntax.models import Syntax, Part
from pydantic import ValidationError
import unittest
from opencloning.syntax import assign_plasmid_to_syntax_part
from pydna.parsers import parse as pydna_parse
from Bio.Restriction import BsaI

test_files = os.path.join(os.path.dirname(__file__), 'test_files')

moclo_syntax = Syntax.model_validate_json(open(os.path.join(test_files, 'syntax', 'moclo_syntax.json')).read())


class TestSyntax(unittest.TestCase):

    def _get_valid_syntax_dict(self):
        """Helper method to create a minimal valid syntax dictionary."""
        return {
            'syntaxName': 'Test Syntax',
            'assemblyEnzyme': 'BsaI',
            'domesticationEnzyme': 'BsmBI',
            'relatedDois': ['10.1000/xyz123'],
            'submitters': ['0000-0000-0000-0000'],
            'overhangNames': {'ACGT': 'test_overhang'},
            'parts': [
                {
                    'id': 1,
                    'name': 'part1',
                    'left_overhang': 'ACGT',
                    'right_overhang': 'CGTA',
                }
            ],
        }

    def test_validate_submitters(self):
        syntax_dict = self._get_valid_syntax_dict()
        Syntax.model_validate(syntax_dict)

        # Invalid ORCID formats
        self.assertRaises(
            ValidationError,
            Syntax.model_validate,
            {**syntax_dict, 'submitters': ['0000-0000-0000']},
        )
        self.assertRaises(
            ValidationError,
            Syntax.model_validate,
            {**syntax_dict, 'submitters': ['invalid']},
        )

    def test_validate_overhang_names(self):
        syntax_dict = self._get_valid_syntax_dict()
        Syntax.model_validate(syntax_dict)

        # Invalid overhang names (not 4 characters)
        self.assertRaises(
            ValidationError,
            Syntax.model_validate,
            {**syntax_dict, 'overhangNames': {'ACG': 'test'}},
        )
        self.assertRaises(
            ValidationError,
            Syntax.model_validate,
            {**syntax_dict, 'overhangNames': {'ACGTG': 'test'}},
        )

    def test_validate_parts_unique_ids(self):
        syntax_dict = self._get_valid_syntax_dict()
        Syntax.model_validate(syntax_dict)

        # Invalid parts with duplicate IDs
        self.assertRaises(
            ValidationError,
            Syntax.model_validate,
            {
                **syntax_dict,
                'parts': [
                    {
                        'id': 1,
                        'name': 'part1',
                        'left_overhang': 'ACGT',
                        'right_overhang': 'CGTA',
                    },
                    {
                        'id': 1,
                        'name': 'part2',
                        'left_overhang': 'TTTT',
                        'right_overhang': 'AAAA',
                    },
                ],
            },
        )

    def test_get_assembly_enzyme(self):
        self.assertEqual(moclo_syntax.get_assembly_enzyme(), BsaI)


class TestPart(unittest.TestCase):
    def test_part_validation(self):
        part_dict = {
            'id': 1,
            'name': '1',
            'info': 'Assembly connector',
            'glyph': 'three-prime-sticky-restriction-site',
            'left_overhang': 'CCCT',
            'right_overhang': 'AACG',
            'left_inside': 'A',
            'right_inside': 'ATTTTTTTT',
            'left_codon_start': 0,
            'right_codon_start': 0,
            'color': '#84c5de',
        }
        Part.model_validate(part_dict)

        # Invalid overhangs or inside
        for side in ['left', 'right']:
            self.assertRaises(ValidationError, Part.model_validate, {**part_dict, f"{side}_overhang": 'CCCTT'})
            self.assertRaises(ValidationError, Part.model_validate, {**part_dict, f"{side}_overhang": 'NNNN'})
            self.assertRaises(ValidationError, Part.model_validate, {**part_dict, f"{side}_overhang": ''})
            self.assertRaises(ValidationError, Part.model_validate, {**part_dict, f"{side}_inside": 'NNNN'})

            Part.model_validate({**part_dict, f"{side}_inside": 'ATTTTTTTT'})
            Part.model_validate({**part_dict, f"{side}_inside": ''})
            Part.model_validate({**part_dict, f"{side}_codon_start": 0})
            Part.model_validate({**part_dict, f"{side}_codon_start": 1000})

            # Invalid codon start
            self.assertRaises(ValidationError, Part.model_validate, {**part_dict, f"{side}_codon_start": -1})

        # Invalid color
        self.assertRaises(ValidationError, Part.model_validate, {**part_dict, 'color': 'invalid'})

    def test_default_values(self):
        part = Part(
            id=1,
            name='1',
            left_overhang='CCCT',
            right_overhang='AACG',
        )
        self.assertEqual(part.info, '')
        self.assertEqual(part.glyph, '')
        self.assertEqual(part.left_inside, '')
        self.assertEqual(part.right_inside, '')
        self.assertEqual(part.left_codon_start, 0)
        self.assertEqual(part.right_codon_start, 0)
        self.assertEqual(part.color, '')


class TestAssignPlasmidToSyntaxPart(unittest.TestCase):
    def test_assign_plasmid_to_syntax_part(self):

        # Has a part that is specifically named
        plasmid1 = pydna_parse('tests/test_files/syntax/pYTK002.gb')[0]
        # Spans multiple parts
        plasmid2 = pydna_parse('tests/test_files/syntax/pYTK095.gb')[0]
        # Has 2 parts
        plasmid3 = pydna_parse('tests/test_files/syntax/moclo_ytk_multi_part.gb')[0]
        # Has no parts that match the syntax
        plasmid4 = pydna_parse('tests/test_files/pAG25.gb')[0]

        self.assertEqual(assign_plasmid_to_syntax_part(plasmid1, moclo_syntax), ['CCCT-AACG'])
        self.assertEqual(assign_plasmid_to_syntax_part(plasmid2, moclo_syntax), ['TACA-CCCT'])
        self.assertEqual(assign_plasmid_to_syntax_part(plasmid3, moclo_syntax), ['ATCC-TGGC', 'CCCT-AACG'])
        self.assertEqual(assign_plasmid_to_syntax_part(plasmid4, moclo_syntax), [])
