from opencloning_linkml.datamodel.models import TextFileSequence, Source
from opencloning_db.models import SequenceType
from opencloning.dna_functions import read_dsrecord_from_json


def guess_sequence_type(sequence: TextFileSequence, source: Source) -> SequenceType:
    seqrecord = read_dsrecord_from_json(sequence)

    if seqrecord.circular:
        return SequenceType.plasmid

    source_type = source.type
    if source_type == 'PCRSource':
        return SequenceType.pcr_product
    elif source_type == 'GenomeCoordinatesSource':
        return SequenceType.locus
    elif source_type in ['HomologousRecombinationSource', 'CRISPRSource', 'RecombinaseSource']:
        return SequenceType.allele
    elif source_type == 'RestrictionEnzymeDigestionSource':
        return SequenceType.restriction_fragment
    else:
        return SequenceType.linear_dna
