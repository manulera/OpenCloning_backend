from unittest import TestCase
from opencloning.pydantic_models import AssemblySource, SequenceLocationStr, AssemblyFragment
from pydna.assembly2 import edge_representation2subfragment_representation
from Bio.SeqFeature import SimpleLocation


class DummyFragment:
    def __init__(self, id):
        self.id = id


class AssemblySourceTest(TestCase):

    def test_from_assembly(self):
        assemblies = [
            [
                (1, 2, SimpleLocation(0, 10), SimpleLocation(10, 20)),
                (2, 3, SimpleLocation(0, 10), SimpleLocation(10, 20)),
            ],
            [
                (1, -2, SimpleLocation(0, 10), SimpleLocation(10, 20)),
                (-2, 3, SimpleLocation(0, 10), SimpleLocation(10, 20)),
            ],
        ]

        for i, assembly in enumerate(assemblies):

            fragments = [DummyFragment(4), DummyFragment(5), DummyFragment(6)]
            assembly_source = AssemblySource.from_assembly(
                assembly=assembly, fragments=fragments, id=0, circular=False
            )
            fragment_assembly = edge_representation2subfragment_representation(assembly, False)

            if i == 0:
                # Check first fragment
                self.assertEqual(assembly_source.input[0].sequence, 4)
                self.assertEqual(assembly_source.input[0].reverse_complemented, False)
                self.assertEqual(assembly_source.input[0].left_location, None)
                self.assertEqual(
                    assembly_source.input[0].right_location,
                    SequenceLocationStr.from_start_and_end(start=0, end=10),
                )

                # Check second fragment
                self.assertEqual(assembly_source.input[1].sequence, 5)
                self.assertEqual(assembly_source.input[1].reverse_complemented, False)
                self.assertEqual(
                    assembly_source.input[1].left_location,
                    SequenceLocationStr.from_start_and_end(start=10, end=20),
                )
                self.assertEqual(
                    assembly_source.input[1].right_location,
                    SequenceLocationStr.from_start_and_end(start=0, end=10),
                )

            for obj, tup in zip(assembly_source.input, fragment_assembly):
                self.assertEqual(obj.to_fragment_tuple(fragments), tup)

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
