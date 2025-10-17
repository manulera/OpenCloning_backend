from opencloning_linkml.datamodel import AssemblySource
from Bio.SeqFeature import Location


def is_assembly_complete(source: AssemblySource) -> bool:
    return any(f.type == 'AssemblyFragment' for f in source.input)


def minimal_assembly_overlap(source: AssemblySource) -> int:
    all_overlaps = list()
    for f in source.input:
        if f.left_location is not None:
            all_overlaps.append(len(Location.fromstring(f.left_location)))
        if f.right_location is not None:
            all_overlaps.append(len(Location.fromstring(f.right_location)))
    return min(all_overlaps)
