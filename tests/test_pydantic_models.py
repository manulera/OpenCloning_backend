from unittest import TestCase
from opencloning_linkml.datamodel import AssemblySource, AssemblyFragment
from pydna.opencloning_models import SequenceLocationStr
from Bio.SeqFeature import SimpleLocation


class DummyFragment:
    def __init__(self, id):
        self.id = id


class AssemblySourceTest(TestCase):

    def test_get_assembly_plan(self):
        # Linear assembly
        assembly = (
            (1, 2, SimpleLocation(0, 10), SimpleLocation(10, 20)),
            (2, 3, SimpleLocation(0, 10), SimpleLocation(10, 20)),
        )

        fragments = [DummyFragment(4), DummyFragment(5), DummyFragment(6)]

        assembly_source = AssemblySource.from_assembly(assembly=assembly, fragments=fragments, id=0, circular=False)
        assembly_plan = assembly_source.get_assembly_plan(fragments)
        self.assertEqual(assembly_plan, assembly)

        # Circular assembly
        assembly = (
            (1, 2, SimpleLocation(0, 10), SimpleLocation(10, 20)),
            (2, 3, SimpleLocation(0, 10), SimpleLocation(10, 20)),
            (3, 1, SimpleLocation(0, 10), SimpleLocation(10, 20)),
        )

        fragments = [DummyFragment(4), DummyFragment(5), DummyFragment(6)]

        assembly_source = AssemblySource.from_assembly(assembly=assembly, fragments=fragments, id=0, circular=True)
        assembly_plan = assembly_source.get_assembly_plan(fragments)
        self.assertEqual(assembly_plan, assembly)


class AssemblyFragmentTest(TestCase):
    def test_field_validator(self):
        # No locations
        AssemblyFragment(sequence=1, reverse_complemented=False, left_location=None, right_location=None)
        location_str = SequenceLocationStr.from_start_and_end(start=0, end=10)
        # SequenceLocationStr
        AssemblyFragment(
            sequence=1, reverse_complemented=False, left_location=location_str, right_location=location_str
        )
        # Normal string
        AssemblyFragment(
            sequence=1, reverse_complemented=False, left_location=str(location_str), right_location=str(location_str)
        )
