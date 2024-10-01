from primer_design import homologous_recombination_primers, gibson_assembly_primers, restriction_enzyme_primers
from Bio.SeqFeature import SimpleLocation, SeqFeature
from unittest import TestCase
from pydna.dseqrecord import Dseqrecord
from pydantic_models import PrimerModel
from pydna.amplify import pcr
from assembly2 import Assembly, gibson_overlap
import pytest
from Bio.Data.IUPACData import ambiguous_dna_values


class TestHomologousRecombinationPrimers(TestCase):

    def test_normal_examples(self):
        """
        Test the homologous_recombination_primers function.
        """

        #                      >>>>>>>>>>> <<<<<<<<<<
        pcr_seq = Dseqrecord('GAAATGGAACAGTGCCAGAAATTTTT')
        pcr_loc = SimpleLocation(1, 23)
        hr_seq = Dseqrecord('AAACGTTT')
        homology_length = 3
        minimal_hybridization_length = 10
        insert_forward = True
        target_tm = 30

        # First we replace the CG
        hr_loc_replace = SimpleLocation(3, 5)
        primers = homologous_recombination_primers(
            pcr_seq,
            pcr_loc,
            hr_seq,
            hr_loc_replace,
            homology_length,
            minimal_hybridization_length,
            insert_forward,
            target_tm,
        )
        self.assertEqual(primers, ('aaaAAATGGAACAG', 'aaaAATTTCTGGC'))

        # An insertion
        hr_loc_insert = SimpleLocation(3, 3)
        primers = homologous_recombination_primers(
            pcr_seq,
            pcr_loc,
            hr_seq,
            hr_loc_insert,
            homology_length,
            minimal_hybridization_length,
            insert_forward,
            target_tm,
        )
        self.assertEqual(primers, ('aaaAAATGGAACAG', 'acgAATTTCTGGC'))

        # Same, circular and we loop
        pcr_seq = pcr_seq.looped()
        hr_seq = hr_seq.looped()
        pcr_seq.features.append(SeqFeature(pcr_loc))
        hr_seq.features.append(SeqFeature(hr_loc_replace))
        hr_seq.features.append(SeqFeature(hr_loc_insert))

        for shift_pcr in range(len(pcr_seq)):
            pcr_shifted = pcr_seq.shifted(shift_pcr)
            pcr_loc = pcr_shifted.features[0].location

            hr_shifted = hr_seq
            hr_loc_replace = hr_shifted.features[0].location
            hr_loc_insert = hr_shifted.features[1].location
            solutions = (
                ('aaaAAATGGAACAG', 'aaaAATTTCTGGC'),
                ('aaaAAATGGAACAG', 'acgAATTTCTGGC'),
            )
            for shift_hr in range(len(hr_seq)):
                hr_shifted = hr_seq.shifted(shift_hr)
                hr_loc_replace = hr_shifted.features[0].location
                hr_loc_insert = hr_shifted.features[1].location
                solutions = (
                    ('aaaAAATGGAACAG', 'aaaAATTTCTGGC'),
                    ('aaaAAATGGAACAG', 'acgAATTTCTGGC'),
                )
                for hr_loc, solution in zip([hr_loc_replace, hr_loc_insert], solutions):
                    for insert_forward in (True, False):
                        if not insert_forward:
                            solution = (
                                solution[0][:homology_length] + solution[1][homology_length:],
                                solution[1][:homology_length] + solution[0][homology_length:],
                            )

                        primers = homologous_recombination_primers(
                            pcr_shifted,
                            pcr_loc,
                            hr_shifted,
                            hr_loc,
                            homology_length,
                            minimal_hybridization_length,
                            insert_forward,
                            target_tm,
                        )
                        self.assertEqual(primers, solution)
        # With spacers
        spacers = ['attt', 'cggg']
        primers = homologous_recombination_primers(
            pcr_shifted,
            pcr_loc,
            hr_shifted,
            hr_loc,
            homology_length,
            minimal_hybridization_length,
            True,  # insert_forward
            target_tm,
            spacers,
        )

        solution = ('aaaatttAAATGGAACAG', 'acgcccgAATTTCTGGC')
        self.assertEqual(primers, solution)

        # The spacer is defined with respect to the locus, not the insert
        # so if we insert inverse, it should look like this:

        primers = homologous_recombination_primers(
            pcr_shifted,
            pcr_loc,
            hr_shifted,
            hr_loc,
            homology_length,
            minimal_hybridization_length,
            False,  # insert_forward
            target_tm,
            spacers,
        )

        solution = ('aaaatttAATTTCTGGC', 'acgcccgAAATGGAACAG')
        self.assertEqual(primers, solution)

    def test_clashing_homology(self):

        pcr_seq = Dseqrecord('GAAATGGAACAGTGCCAGAAATTTTT')
        pcr_loc = SimpleLocation(1, 23)
        hr_seq = Dseqrecord('AAACGTTT')
        homology_length = 5
        minimal_hybridization_length = 10
        insert_forward = True
        target_tm = 30
        hr_loc_insert = SimpleLocation(3, 3)

        try:
            homologous_recombination_primers(
                pcr_seq,
                pcr_loc,
                hr_seq,
                hr_loc_insert,
                homology_length,
                minimal_hybridization_length,
                insert_forward,
                target_tm,
            )
            self.fail('Expected ValueError.')

        except ValueError as e:
            self.assertEqual(str(e), 'Forward homology region is out of bounds.')

        hr_loc_insert = SimpleLocation(5, 5)

        try:
            homologous_recombination_primers(
                pcr_seq,
                pcr_loc,
                hr_seq,
                hr_loc_insert,
                homology_length,
                minimal_hybridization_length,
                insert_forward,
                target_tm,
            )
            self.fail('Expected ValueError.')
        except ValueError as e:
            self.assertEqual(str(e), 'Reverse homology region is out of bounds.')

        # circular case with clashing homology
        hr_seq = hr_seq.looped()
        homology_length = 4
        # exact match
        homologous_recombination_primers(
            pcr_seq,
            pcr_loc,
            hr_seq,
            hr_loc_insert,
            homology_length,
            minimal_hybridization_length,
            insert_forward,
            target_tm,
        )
        # clashing homology
        homology_length = 5
        try:
            homologous_recombination_primers(
                pcr_seq,
                pcr_loc,
                hr_seq,
                hr_loc_insert,
                homology_length,
                minimal_hybridization_length,
                insert_forward,
                target_tm,
            )
            self.fail('Expected ValueError.')
        except ValueError as e:
            self.assertEqual(str(e), 'Homology arms overlap.')

    def test_errors(self):
        pcr_seq = Dseqrecord('GAAATGGAACAGTGCCAGAAATTTTT')
        pcr_loc = SimpleLocation(1, 23)
        hr_seq = Dseqrecord('AAACGTTT')
        homology_length = 3
        minimal_hybridization_length = 10
        insert_forward = True
        target_tm = 30
        hr_loc_insert = SimpleLocation(5, 5)

        # Wrong number of spacers
        try:
            homologous_recombination_primers(
                pcr_seq,
                pcr_loc,
                hr_seq,
                hr_loc_insert,
                homology_length,
                minimal_hybridization_length,
                insert_forward,
                target_tm,
                spacers=['ATTT'],
            )
            self.fail('Expected ValueError.')

        except ValueError as e:
            self.assertIn("The 'spacers' list", str(e))

        # Primers can't be designed

        try:
            homologous_recombination_primers(
                pcr_seq,
                pcr_loc,
                hr_seq,
                hr_loc_insert,
                homology_length,
                100,
                insert_forward,
                target_tm,
            )
            self.fail('Expected ValueError.')
        except ValueError as e:
            self.assertIn('Primers could not be designed', str(e))


class TestGibsonAssemblyPrimers(TestCase):

    def test_normal_examples(self):
        """
        Test the gibson_assembly_primers function.
        """

        # Test case for gibson_assembly_primers function
        templates = [
            Dseqrecord('AAACAGTAATACGTTCCTTTTTTATGATGATGGATGACATTCAAAGCACTGATTCTAT'),
            Dseqrecord('GTTTACAACGGCAATGAACGTTCCTTTTTTATGATATGCCCAGCTTCATGAAATGGAA'),
            Dseqrecord('AAGGACAACGTTCCTTTTTTATGATATATATGGCACAGTATGATCAAAAGTTAAGTAC'),
        ]
        # Set the name to 'name' to check formatting of name
        for i, template in enumerate(templates):
            template.name = 'name'
            template.id = f'{i}'

        homology_length = 20
        minimal_hybridization_length = 15
        target_tm = 55
        circular = True

        primers = gibson_assembly_primers(
            templates, homology_length, minimal_hybridization_length, target_tm, circular
        )

        # Check if the correct number of primers is returned
        self.assertEqual(len(primers), 6)  # 2 primers per template

        # Check if all primers are instances of PrimerModel
        for primer in primers:
            self.assertIsInstance(primer, PrimerModel)

        # Check if primer names are correctly formatted
        for i, primer in enumerate(primers):
            template_index = i // 2
            primer_type = 'fwd' if i % 2 == 0 else 'rvs'
            expected_name = f'seq_{templates[template_index].id}_{primer_type}'
            self.assertEqual(primer.name, expected_name)

        # Test that it yields the right sequence
        amplicons = list()
        for i in range(len(templates)):
            amplicons.append(pcr(primers[i * 2].to_pydna_primer(), primers[i * 2 + 1].to_pydna_primer(), templates[i]))

        asm = Assembly(
            amplicons,
            algorithm=gibson_overlap,
            use_fragment_order=False,
            use_all_fragments=True,
            limit=homology_length,
        )

        circular_assemblies = asm.assemble_circular()
        self.assertEqual(len(circular_assemblies), 1)
        self.assertEqual(circular_assemblies[0].seq.seguid(), sum(templates, Dseqrecord('')).seq.looped().seguid())

        # Test with circular=False
        circular = False
        primers = gibson_assembly_primers(
            templates, homology_length, minimal_hybridization_length, target_tm, circular
        )

        # Check if the correct number of primers is returned for linear assembly
        self.assertEqual(len(primers), 6)  # 2 primers per template

        # Test that it yields the right sequence
        amplicons = list()
        for i in range(len(templates)):
            amplicons.append(pcr(primers[i * 2].to_pydna_primer(), primers[i * 2 + 1].to_pydna_primer(), templates[i]))

        asm = Assembly(
            amplicons,
            algorithm=gibson_overlap,
            use_fragment_order=False,
            use_all_fragments=True,
            limit=homology_length,
        )

        linear_assemblies = asm.assemble_linear()
        self.assertEqual(len(linear_assemblies), 1)
        self.assertEqual(linear_assemblies[0].seq, sum(templates, Dseqrecord('')).seq)

    @pytest.mark.xfail(
        reason='Waiting on https://github.com/BjornFJohansson/pydna/issues/265 and https://github.com/BjornFJohansson/pydna/issues/264'
    )
    def test_primer_errors(self):
        """
        Test the gibson_assembly_primers function.
        """

        templates = [
            Dseqrecord('AAACAGTAATACGTTCCTTTTTTATGATGATGGATGACATTCAAAGCACTGATTCTAT'),
            Dseqrecord('GTTTACTGGAA'),
            Dseqrecord('AAGGACAACGTTCCTTTTTTATGATATATATGGCACAGTATGATCAAAAGTTAAGTAC'),
        ]
        homology_length = 20
        minimal_hybridization_length = 15
        target_tm = 55
        circular = True

        try:
            gibson_assembly_primers(templates, homology_length, minimal_hybridization_length, target_tm, circular)
            self.fail('Expected ValueError.')
        except ValueError as e:
            self.assertEqual(str(e), 'Primers could not be designed for template 2, try changing settings.')

        templates = [templates[0], templates[2]]

        # Now ask for too long minimal_hybridization_length
        minimal_hybridization_length = 100
        try:
            gibson_assembly_primers(templates, homology_length, minimal_hybridization_length, target_tm, circular)
            self.fail('Expected ValueError.')
        except ValueError as e:
            self.assertEqual(str(e), 'Primers could not be designed for template 1, 2, try changing settings.')

        # Too long homology_length
        homology_length = 200
        minimal_hybridization_length = 15
        try:
            gibson_assembly_primers(templates, homology_length, minimal_hybridization_length, target_tm, circular)
            self.fail('Expected ValueError.')
        except ValueError as e:
            self.assertEqual(str(e), 'Primers could not be designed for template 1, 2, try changing settings.')


class TestRestrictionEnzymePrimers(TestCase):
    def test_restriction_enzyme_primers(self):
        """
        Test the restriction_enzyme_primers function.
        """
        from Bio.Restriction import EcoRI, BamHI, AflIII

        template = Dseqrecord('ATGCATGCATGCATGCATGCATGC')
        minimal_hybridization_length = 10
        target_tm = 55
        left_enzyme = EcoRI
        right_enzyme = BamHI
        filler_bases = 'GC'

        fwd, rvs = restriction_enzyme_primers(
            template, minimal_hybridization_length, target_tm, left_enzyme, right_enzyme, filler_bases
        )

        # Check that primers contain the correct restriction sites
        self.assertTrue(fwd.sequence.startswith('GC' + str(EcoRI.site)))
        self.assertTrue(rvs.sequence.startswith('GC' + str(BamHI.site)))

        # Test with only left enzyme
        fwd, rvs = restriction_enzyme_primers(
            template, minimal_hybridization_length, target_tm, left_enzyme, None, filler_bases
        )

        self.assertTrue(fwd.sequence.startswith('GC' + str(EcoRI.site)))
        self.assertFalse(rvs.sequence.startswith('GC' + str(BamHI.site)))

        # Test with only right enzyme
        fwd, rvs = restriction_enzyme_primers(
            template, minimal_hybridization_length, target_tm, None, right_enzyme, filler_bases
        )

        self.assertFalse(fwd.sequence.startswith('GC' + str(EcoRI.site)))
        self.assertTrue(rvs.sequence.startswith('GC' + str(BamHI.site)))

        # Test with no enzymes
        fwd, rvs = restriction_enzyme_primers(
            template, minimal_hybridization_length, target_tm, None, None, filler_bases
        )

        self.assertFalse(fwd.sequence.startswith('GC' + str(EcoRI.site)))
        self.assertFalse(rvs.sequence.startswith('GC' + str(BamHI.site)))

        # Test with enzyme that has ambiguous bases
        fwd, rvs = restriction_enzyme_primers(
            template, minimal_hybridization_length, target_tm, AflIII, AflIII, filler_bases
        )
        actual_site = (
            str(AflIII.site).replace('R', ambiguous_dna_values['R'][0]).replace('Y', ambiguous_dna_values['Y'][0])
        )

        self.assertTrue(fwd.sequence.startswith('GC' + actual_site))
        self.assertTrue(rvs.sequence.startswith('GC' + actual_site))

        # Test spacers on one side only
        spacers = ['ACGT' * 10, '']
        fwd, rvs = restriction_enzyme_primers(
            template, minimal_hybridization_length, target_tm, left_enzyme, right_enzyme, filler_bases, spacers=spacers
        )

        self.assertTrue('GC' + str(left_enzyme.site) + 'ACGT' * 10 in fwd.sequence)
        self.assertTrue('GC' + str(right_enzyme.site) in rvs.sequence)

        # Both spacers
        spacers = ['ACGT' * 10, 'ACGT' * 10]
        fwd, rvs = restriction_enzyme_primers(
            template, minimal_hybridization_length, target_tm, left_enzyme, right_enzyme, filler_bases, spacers=spacers
        )

        self.assertTrue('GC' + str(left_enzyme.site) + 'ACGT' * 10 in fwd.sequence)
        self.assertTrue('GC' + str(right_enzyme.site) + 'ACGT' * 10 in rvs.sequence)
