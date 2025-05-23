from fastapi import Body, Query, HTTPException
from pydantic import create_model
import re
from Bio.Restriction import RestrictionBatch
from Bio.SeqUtils import gc_fraction
from primer3 import bindings
from itertools import product

from ..dna_functions import get_invalid_enzyme_names
from ..pydantic_models import PrimerModel, PrimerDesignQuery
from ..dna_functions import read_dsrecord_from_json
from ..primer_design import (
    homologous_recombination_primers,
    gibson_assembly_primers,
    simple_pair_primers,
)
from ..get_router import get_router
from ..ebic.primer_design import ebic_primers
from pydantic import BaseModel

router = get_router()

PrimerDesignResponse = create_model('PrimerDesignResponse', primers=(list[PrimerModel], ...))


def validate_spacers(spacers: list[str] | None, nb_templates: int, circular: bool):
    if spacers is None:
        return

    if circular and len(spacers) != nb_templates:
        raise HTTPException(
            422, 'The number of spacers must be the same as the number of templates when the assembly is circular.'
        )
    if not circular and len(spacers) != (nb_templates + 1):
        raise HTTPException(
            422, 'The number of spacers must be one more than the number of templates when the assembly is linear.'
        )
    for spacer in spacers:
        # If it's not only ACGt
        if not re.match(r'^[ACGT]*$', spacer.upper()):
            raise HTTPException(422, 'Spacer can only contain ACGT bases.')


@router.post('/primer_design/homologous_recombination', response_model=PrimerDesignResponse)
async def primer_design_homologous_recombination(
    pcr_template: PrimerDesignQuery,
    homologous_recombination_target: PrimerDesignQuery,
    spacers: list[str] | None = Body(
        None,
        description='Spacers to add at the left and right side of the insertion.',
    ),
    homology_length: int = Query(..., description='The length of the homology region in bps.'),
    minimal_hybridization_length: int = Query(
        ..., description='The minimal length of the hybridization region in bps.'
    ),
    target_tm: float = Query(
        ..., description='The desired melting temperature for the hybridization part of the primer.'
    ),
):
    """Design primers for homologous recombination"""

    validate_spacers(spacers, 1, False)

    pcr_seq = read_dsrecord_from_json(pcr_template.sequence)
    pcr_loc = pcr_template.location.to_biopython_location()

    hr_seq = read_dsrecord_from_json(homologous_recombination_target.sequence)
    hr_loc = homologous_recombination_target.location.to_biopython_location()

    insert_forward = pcr_template.forward_orientation

    try:
        forward_primer, reverse_primer = homologous_recombination_primers(
            pcr_seq,
            pcr_loc,
            hr_seq,
            hr_loc,
            homology_length,
            minimal_hybridization_length,
            insert_forward,
            target_tm,
            spacers,
        )
    except ValueError as e:
        raise HTTPException(400, *e.args)

    forward_primer = PrimerModel(id=0, sequence=forward_primer, name='forward')
    reverse_primer = PrimerModel(id=0, sequence=reverse_primer, name='reverse')

    return {'primers': [forward_primer, reverse_primer]}


@router.post('/primer_design/gibson_assembly', response_model=PrimerDesignResponse)
async def primer_design_gibson_assembly(
    pcr_templates: list[PrimerDesignQuery],
    spacers: list[str] | None = Body(
        None,
        description='Spacers to add between the restriction site and the 5\' end of the primer footprint (the part that binds the DNA).',
    ),
    homology_length: int = Query(..., description='The length of the homology region in bps.'),
    minimal_hybridization_length: int = Query(
        ..., description='The minimal length of the hybridization region in bps.'
    ),
    target_tm: float = Query(
        ..., description='The desired melting temperature for the hybridization part of the primer.'
    ),
    circular: bool = Query(False, description='Whether the assembly is circular.'),
):
    """Design primers for Gibson assembly"""
    # Validate the spacers
    validate_spacers(spacers, len(pcr_templates), circular)
    templates = list()
    for query in pcr_templates:
        dseqr = read_dsrecord_from_json(query.sequence)
        location = query.location.to_biopython_location()
        template = location.extract(dseqr)
        if not query.forward_orientation:
            template = template.reverse_complement()
        # For naming the primers
        template.name = dseqr.name
        template.id = dseqr.id
        templates.append(template)
    try:
        primers = gibson_assembly_primers(
            templates, homology_length, minimal_hybridization_length, target_tm, circular, spacers
        )
    except ValueError as e:
        raise HTTPException(400, *e.args)

    return {'primers': primers}


@router.post('/primer_design/simple_pair', response_model=PrimerDesignResponse)
async def primer_design_simple_pair(
    pcr_template: PrimerDesignQuery,
    spacers: list[str] | None = Body(
        None,
        description='Spacers to add between the restriction site and the 5\' end of the primer footprint (the part that binds the DNA).',
    ),
    minimal_hybridization_length: int = Query(
        ..., description='The minimal length of the hybridization region in bps.'
    ),
    target_tm: float = Query(
        ..., description='The desired melting temperature for the hybridization part of the primer.'
    ),
    left_enzyme: str | None = Query(None, description='The restriction enzyme for the left side of the sequence.'),
    right_enzyme: str | None = Query(None, description='The restriction enzyme for the right side of the sequence.'),
    left_enzyme_inverted: bool = Query(False, description='Whether the left enzyme site is inverted.'),
    right_enzyme_inverted: bool = Query(False, description='Whether the right enzyme site is inverted.'),
    filler_bases: str = Query(
        'TTT',
        description='These bases are added to the 5\' end of the primer to ensure proper restriction enzyme digestion.',
    ),
):
    """
    Design primers for restriction ligation. If the restriction site contains ambiguous bases, the primer bases will
    be chosen by order of appearance in the dictionary `ambiguous_dna_values` from `Bio.Data.IUPACData`.
    """

    if not re.match(r'^[ACGT]*$', filler_bases.upper()):
        raise HTTPException(400, 'Filler bases can only contain ACGT bases.')

    invalid_enzymes = get_invalid_enzyme_names([left_enzyme, right_enzyme])
    if len(invalid_enzymes):
        raise HTTPException(404, 'These enzymes do not exist: ' + ', '.join(invalid_enzymes))

    validate_spacers(spacers, 1, False)

    dseqr = read_dsrecord_from_json(pcr_template.sequence)
    location = pcr_template.location.to_biopython_location()
    template = location.extract(dseqr)
    if not pcr_template.forward_orientation:
        template = template.reverse_complement()
    template.name = dseqr.name
    template.id = dseqr.id
    # This is to my knowledge the only way to get the enzymes
    rb = RestrictionBatch()
    try:
        fwd, rvs = simple_pair_primers(
            template,
            minimal_hybridization_length,
            target_tm,
            rb.format(left_enzyme) if left_enzyme is not None else None,
            rb.format(right_enzyme) if right_enzyme is not None else None,
            filler_bases,
            spacers,
            left_enzyme_inverted,
            right_enzyme_inverted,
        )
    except ValueError as e:
        raise HTTPException(400, *e.args)

    return {'primers': [fwd, rvs]}


@router.post('/primer_design/ebic', response_model=PrimerDesignResponse)
async def primer_design_ebic(
    template: PrimerDesignQuery,
    max_inside: int = Query(..., description='The maximum length of the inside edge of the EBIC primer.'),
    max_outside: int = Query(..., description='The maximum length of the outside edge of the EBIC primer.'),
):
    """Design primers for EBIC"""
    dseqr = read_dsrecord_from_json(template.sequence)
    location = template.location.to_biopython_location()
    return {'primers': ebic_primers(dseqr, location, max_inside, max_outside)}


# @router.post('/primer_design/gateway_attB', response_model=PrimerDesignResponse)
# async def primer_design_gateway_attB(
#     template: PrimerDesignQuery,
#     left_site: str = Query(..., description='The left attB site to recombine.', pattern=r'^attB[1-5]$'),
#     right_site: str = Query(..., description='The right attB site to recombine.', pattern=r'^attB[1-5]$'),
#     spacers: list[str] | None = Body(None, description='Spacers to add between the attB site and the primer.'),
#     filler_bases: str = Query(
#         'GGGG',
#         description='These bases are added to the 5\' end of the primer to ensure proper restriction enzyme digestion.',
#     ),
#     minimal_hybridization_length: int = Query(
#         ..., description='The minimal length of the hybridization region in bps.'
#     ),
#     target_tm: float = Query(
#         ..., description='The desired melting temperature for the hybridization part of the primer.'
#     ),
# ):
#     """Design primers for Gateway attB"""
#     dseqr = read_dsrecord_from_json(template.sequence)
#     location = template.location.to_biopython_location(dseqr.circular, len(dseqr))
#     template = location.extract(dseqr)
#     if not template.forward_orientation:
#         template = template.reverse_complement()
#     template.name = dseqr.name
#     template.id = dseqr.id
#     try:
#         primers = gateway_attB_primers(
#             template, minimal_hybridization_length, target_tm, (left_site, right_site), spacers, filler_bases
#         )
#     except ValueError as e:
#         raise HTTPException(400, *e.args)

#     return {'primers': primers}


class ThermodynamicResult(BaseModel):
    melting_temperature: float
    deltaG: float
    figure: str | None

    @classmethod
    def from_binding(cls, result):
        return cls(
            melting_temperature=result.tm,
            deltaG=result.dg,
            figure='\n'.join(result.ascii_structure_lines),
        )


def get_sequence_thermodynamic_result(sequences: list[str], method: callable) -> ThermodynamicResult | None:
    """Get the thermodynamic result for a sequence, if the sequence is longer than primer3 60bp limit, it will be split into two
    and the result with the lowest deltaG will be returned."""
    results = [method(seq, output_structure=True) for seq in sequences]
    results = [r for r in results if r.structure_found]
    if len(results) == 0:
        return None

    result = min(results, key=lambda r: r.dg)
    return ThermodynamicResult.from_binding(result)


class PrimerDetailsResponse(BaseModel):
    melting_temperature: float
    gc_content: float
    homodimer: ThermodynamicResult | None
    hairpin: ThermodynamicResult | None


@router.get('/primer_details', response_model=PrimerDetailsResponse)
async def primer_details(
    sequence: str = Query(..., description='Primer sequence', pattern=r'^[ACGTacgt]+$'),
):
    """Get information about a primer"""
    sequence = sequence.upper()
    tm = bindings.calc_tm(sequence)
    gc_content = gc_fraction(sequence)

    thermodynamic_sequences = [sequence]
    if len(sequence) > 60:
        thermodynamic_sequences = [sequence[:60], sequence[-60:]]
    homodimer = get_sequence_thermodynamic_result(thermodynamic_sequences, bindings.calc_homodimer)
    hairpin = get_sequence_thermodynamic_result(thermodynamic_sequences, bindings.calc_hairpin)

    return {
        'melting_temperature': tm,
        'gc_content': gc_content,
        'homodimer': homodimer,
        'hairpin': hairpin,
    }


@router.get('/primer_heterodimer', response_model=ThermodynamicResult | None)
async def primer_heterodimer(
    sequence1: str = Query(..., description='First primer sequence', pattern=r'^[ACGTacgt]+$'),
    sequence2: str = Query(..., description='Second primer sequence', pattern=r'^[ACGTacgt]+$'),
):
    """Get information about a primer pair"""

    if len(sequence1) <= 60 or len(sequence2) <= 60:
        sequence_pairs = [(sequence1, sequence2)]
    else:
        sequence_pairs = list(product((sequence1[:60], sequence1[-60:]), (sequence2[:60], sequence2[-60:])))

    results = [bindings.calc_heterodimer(s1, s2, output_structure=True) for s1, s2 in sequence_pairs]
    results = [r for r in results if r.structure_found]
    if len(results) == 0:
        return None

    result = min(results, key=lambda r: r.dg)
    return ThermodynamicResult.from_binding(result)
