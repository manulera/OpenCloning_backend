from fastapi import Query, HTTPException
from typing import Union, Literal, Callable
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer as PydnaPrimer
from pydna.crispr import cas9
from pydantic import conlist, create_model
from Bio.Restriction.Restriction import RestrictionBatch
from opencloning.cre_lox import cre_loxP_overlap, annotate_loxP_sites
from ..dna_functions import (
    get_invalid_enzyme_names,
    format_sequence_genbank,
    read_dsrecord_from_json,
)
from ..pydantic_models import (
    PCRSource,
    PrimerModel,
    TextFileSequence,
    LigationSource,
    HomologousRecombinationSource,
    CRISPRSource,
    GibsonAssemblySource,
    InFusionSource,
    RestrictionAndLigationSource,
    AssemblySource,
    OverlapExtensionPCRLigationSource,
    GatewaySource,
    CreLoxRecombinationSource,
    InVivoAssemblySource,
)
from ..assembly2 import (
    Assembly,
    assemble,
    sticky_end_sub_strings,
    PCRAssembly,
    gibson_overlap,
    filter_linear_subassemblies,
    restriction_ligation_overlap,
    SingleFragmentAssembly,
    blunt_overlap,
    combine_algorithms,
    annotate_primer_binding_sites,
)

from ..gateway import gateway_overlap, find_gateway_sites, annotate_gateway_sites
from ..get_router import get_router

router = get_router()


def format_known_assembly_response(
    source: AssemblySource,
    out_sources: list[AssemblySource],
    fragments: list[Dseqrecord],
    product_callback: Callable[[Dseqrecord], Dseqrecord] = lambda x: x,
):
    """Common function for assembly sources, when assembly is known"""
    # If a specific assembly is requested
    assembly_plan = source.get_assembly_plan(fragments)
    for s in out_sources:
        # TODO: it seems that assemble() is not getting is_insertion ever
        other_assembly_plan = s.get_assembly_plan(fragments)
        if assembly_plan == other_assembly_plan:
            return {
                'sequences': [
                    format_sequence_genbank(product_callback(assemble(fragments, assembly_plan)), s.output_name)
                ],
                'sources': [s],
            }
    raise HTTPException(400, 'The provided assembly is not valid.')


@router.post(
    '/crispr',
    response_model=create_model(
        'CrisprResponse',
        sources=(list[CRISPRSource], ...),
        sequences=(list[TextFileSequence], ...),
    ),
)
async def crispr(
    source: CRISPRSource,
    guides: list[PrimerModel],
    sequences: conlist(TextFileSequence, min_length=2, max_length=2),
    minimal_homology: int = Query(40, description='The minimum homology between the template and the insert.'),
):
    """Return the sequence after performing CRISPR editing by Homology directed repair
    TODO: Support repair through NHEJ
    TODO: Check support for circular DNA targets
    """
    template, insert = [read_dsrecord_from_json(seq) for seq in sequences]

    # TODO: check input method for guide (currently as a primer)
    # TODO: support user input PAM

    # Check cutsites from guide provided by user
    guide_cuts = []
    for guide in guides:
        enzyme = cas9(guide.sequence)
        possible_cuts = template.seq.get_cutsites(enzyme)
        if len(possible_cuts) == 0:
            raise HTTPException(
                400, f'Could not find Cas9 cutsite in the target sequence using the guide: {guide.name}'
            )
        guide_cuts.append(possible_cuts)

    # Check if homologous recombination is possible
    fragments = [template, insert]
    asm = Assembly(fragments, minimal_homology, use_all_fragments=True)
    try:
        possible_assemblies = [a for a in asm.get_insertion_assemblies() if a[0][0] == 1]
    except ValueError as e:
        raise HTTPException(400, *e.args)

    if not possible_assemblies:
        raise HTTPException(400, 'Repair fragment cannot be inserted in the target sequence through homology')

    valid_assemblies = []
    # Check if Cas9 cut is within the homologous recombination region
    for a in possible_assemblies:
        hr_start = int(a[0][2].start)
        hr_end = int(a[1][3].end)

        for cuts in guide_cuts:
            reparable_cuts = [c for c in cuts if c[0][0] > hr_start and c[0][0] <= hr_end]
            if len(reparable_cuts):
                valid_assemblies.append(a)
            if len(reparable_cuts) != len(cuts):
                # TODO: warning a cutsite falls outside
                pass

    if len(valid_assemblies) == 0:
        raise HTTPException(
            400, 'A Cas9 cutsite was found, and a homologous recombination region, but they do not overlap.'
        )
    # elif len(valid_assemblies) != len(possible_assemblies):
    #     # TODO: warning that some assemblies were discarded
    #     pass

    # TODO: double check that this works for circular DNA -> for now get_insertion_assemblies() is only
    # meant for linear DNA

    out_sources = [
        CRISPRSource.from_assembly(id=source.id, assembly=a, guides=source.guides, fragments=fragments)
        for a in valid_assemblies
    ]

    # If a specific assembly is requested
    if len(source.assembly):
        return format_known_assembly_response(source, out_sources, [template, insert])

    out_sequences = [
        format_sequence_genbank(assemble([template, insert], a), source.output_name) for a in valid_assemblies
    ]
    return {'sources': out_sources, 'sequences': out_sequences}


def generate_assemblies(
    source: AssemblySource,
    create_source: Callable[[list, bool], AssemblySource],
    fragments: list[TextFileSequence],
    circular_only: bool,
    algo: Callable,
    allow_insertion_assemblies: bool,
    assembly_kwargs: dict | None = None,
    product_callback: Callable[[Dseqrecord], Dseqrecord] = lambda x: x,
    recombination_mode: bool = False,
) -> dict[Literal['sources', 'sequences'], list[AssemblySource] | list[TextFileSequence]]:
    if assembly_kwargs is None:
        assembly_kwargs = {}
    try:
        out_sources = []
        if len(fragments) > 1:
            asm = Assembly(
                fragments,
                algorithm=algo,
                use_all_fragments=True,
                use_fragment_order=False,
                **assembly_kwargs,
            )
            circular_assemblies = asm.get_circular_assemblies()
            out_sources += [create_source(a, True) for a in circular_assemblies]
            if not circular_only:
                if not recombination_mode:
                    out_sources += [
                        create_source(a, False)
                        for a in filter_linear_subassemblies(
                            asm.get_linear_assemblies(), circular_assemblies, fragments
                        )
                    ]
                else:
                    out_sources += [create_source(a, False) for a in asm.get_insertion_assemblies()]
        else:
            asm = SingleFragmentAssembly(fragments, algorithm=algo, **assembly_kwargs)
            out_sources.extend(create_source(a, True) for a in asm.get_circular_assemblies())
            if not circular_only and allow_insertion_assemblies:
                out_sources.extend(create_source(a, False) for a in asm.get_insertion_assemblies())

    except ValueError as e:
        raise HTTPException(400, *e.args)

    # If a specific assembly is requested
    if len(source.assembly):
        return format_known_assembly_response(source, out_sources, fragments, product_callback)

    out_sequences = [
        format_sequence_genbank(
            product_callback(assemble(fragments, s.get_assembly_plan(fragments))), source.output_name
        )
        for s in out_sources
    ]

    return {'sources': out_sources, 'sequences': out_sequences}


@router.post(
    '/ligation',
    response_model=create_model(
        'LigationResponse', sources=(list[LigationSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def ligation(
    source: LigationSource,
    sequences: conlist(TextFileSequence, min_length=1),
    blunt: bool = Query(False, description='Use blunt ligation as well as sticky ends.'),
    allow_partial_overlap: bool = Query(False, description='Allow for partially overlapping sticky ends.'),
    circular_only: bool = Query(False, description='Only return circular assemblies.'),
):

    fragments = [read_dsrecord_from_json(seq) for seq in sequences]

    # Lambda function for code clarity
    def create_source(a, is_circular):
        return LigationSource.from_assembly(assembly=a, circular=is_circular, id=source.id, fragments=fragments)

    # If the assembly is known, the blunt parameter is ignored, and we set the algorithm type from the assembly
    # (blunt ligations have features without length)
    if len(source.assembly):
        asm = source.get_assembly_plan(fragments)
        blunt = len(asm[0][2]) == 0

    algo = combine_algorithms(blunt_overlap, sticky_end_sub_strings) if blunt else sticky_end_sub_strings
    resp = generate_assemblies(
        source, create_source, fragments, circular_only, algo, False, {'limit': allow_partial_overlap}
    )
    if len(resp['sources']) == 0:
        raise HTTPException(400, 'No ligations were found.')

    return resp


@router.post(
    '/pcr',
    response_model=create_model(
        'PCRResponse', sources=(list[PCRSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def pcr(
    source: PCRSource,
    sequences: conlist(TextFileSequence, min_length=1, max_length=1),
    primers: conlist(PrimerModel, min_length=1, max_length=2),
    minimal_annealing: int = Query(20, description='The minimal annealing length for each primer.'),
    allowed_mismatches: int = Query(0, description='The number of mismatches allowed'),
):
    if len(primers) != len(sequences) * 2:
        raise HTTPException(400, 'The number of primers should be twice the number of sequences.')

    pydna_sequences = [read_dsrecord_from_json(s) for s in sequences]
    pydna_primers = [PydnaPrimer(p.sequence, id=str(p.id), name=p.name) for p in primers]

    # TODO: This may have to be re-written if we allow mismatches
    # If an assembly is provided, we ignore minimal_annealing
    # What happens if annealing is zero? That would mean
    # mismatch in the 3' of the primer, which maybe should
    # not be allowed.
    if len(source.assembly):
        minimal_annealing = source.minimal_overlap()
        # Only the ones that match are included in the output assembly
        # location, so the submitted assembly should be returned without
        # allowed mistmatches
        # TODO: tests for this
        allowed_mismatches = 0

    # Arrange the fragments in the order primer, sequence, primer
    fragments = list()
    while len(pydna_primers):
        fragments.append(pydna_primers.pop(0))
        fragments.append(pydna_sequences.pop(0))
        fragments.append(pydna_primers.pop(0))

    asm = PCRAssembly(fragments, limit=minimal_annealing, mismatches=allowed_mismatches)
    try:
        possible_assemblies = asm.get_linear_assemblies()
    except ValueError as e:
        raise HTTPException(400, *e.args)

    # Edge case: where both primers are identical, remove
    # duplicate assemblies that represent just reverse complement
    if len(sequences) == 1 and primers[0].id == primers[1].id:
        possible_assemblies = [a for a in possible_assemblies if (a[0][0] == 1 and a[0][1] == 2)]

    out_sources = [
        PCRSource.from_assembly(
            id=source.id,
            assembly=a,
            circular=False,
            fragments=fragments,
            add_primer_features=source.add_primer_features,
        )
        for a in possible_assemblies
    ]

    # If a specific assembly is requested
    if len(source.assembly):

        def callback(x):
            if source.add_primer_features:
                return annotate_primer_binding_sites(x, fragments, source.get_assembly_plan(fragments))
            else:
                return x

        return format_known_assembly_response(source, out_sources, fragments, callback)

    if len(possible_assemblies) == 0:
        raise HTTPException(400, 'No pair of annealing primers was found. Try changing the annealing settings.')

    def callback(fragments, a):
        out_seq = assemble(fragments, a)
        if source.add_primer_features:
            return annotate_primer_binding_sites(out_seq, fragments, possible_assemblies)
        else:
            return out_seq

    out_sequences = [
        format_sequence_genbank(callback(fragments, a), source.output_name)
        for s, a in zip(out_sources, possible_assemblies)
    ]

    return {'sources': out_sources, 'sequences': out_sequences}


@router.post(
    '/homologous_recombination',
    response_model=create_model(
        'HomologousRecombinationResponse',
        sources=(list[HomologousRecombinationSource], ...),
        sequences=(list[TextFileSequence], ...),
    ),
)
async def homologous_recombination(
    source: HomologousRecombinationSource,
    sequences: conlist(TextFileSequence, min_length=2, max_length=2),
    minimal_homology: int = Query(40, description='The minimum homology between the template and the insert.'),
):

    template, insert = [read_dsrecord_from_json(seq) for seq in sequences]

    # If an assembly is provided, we ignore minimal_homology
    if len(source.assembly):
        minimal_homology = source.minimal_overlap()

    asm = Assembly((template, insert), limit=minimal_homology, use_all_fragments=True)

    # The condition is that the first and last fragments are the template
    try:
        if not template.circular:
            possible_assemblies = [a for a in asm.get_insertion_assemblies() if a[0][0] == 1]
        else:
            possible_assemblies = [a for a in asm.get_circular_assemblies()]

    except ValueError as e:
        raise HTTPException(400, *e.args)

    if len(possible_assemblies) == 0:
        raise HTTPException(400, 'No homologous recombination was found.')

    out_sources = [
        HomologousRecombinationSource.from_assembly(
            id=source.id, assembly=a, circular=False, fragments=[template, insert]
        )
        for a in possible_assemblies
    ]

    # If a specific assembly is requested
    if len(source.assembly):
        return format_known_assembly_response(source, out_sources, [template, insert])

    out_sequences = [
        format_sequence_genbank(assemble([template, insert], a), source.output_name) for a in possible_assemblies
    ]

    return {'sources': out_sources, 'sequences': out_sequences}


@router.post(
    '/gibson_assembly',
    response_model=create_model(
        'GibsonAssemblyResponse',
        sources=(
            list[Union[GibsonAssemblySource, OverlapExtensionPCRLigationSource, InFusionSource, InVivoAssemblySource]],
            ...,
        ),
        sequences=(list[TextFileSequence], ...),
    ),
)
async def gibson_assembly(
    sequences: conlist(TextFileSequence, min_length=1),
    source: Union[GibsonAssemblySource, OverlapExtensionPCRLigationSource, InFusionSource, InVivoAssemblySource],
    minimal_homology: int = Query(
        40, description='The minimum homology between consecutive fragments in the assembly.'
    ),
    circular_only: bool = Query(False, description='Only return circular assemblies.'),
):

    fragments = [read_dsrecord_from_json(seq) for seq in sequences]

    # Lambda function for code clarity
    def create_source(a, is_circular):
        return source.__class__.from_assembly(assembly=a, circular=is_circular, id=source.id, fragments=fragments)

    resp = generate_assemblies(
        source, create_source, fragments, circular_only, gibson_overlap, False, {'limit': minimal_homology}
    )

    if len(resp['sources']) == 0:
        raise HTTPException(
            400,
            f'No {"circular " if circular_only else ""}assembly with at least {minimal_homology} bps of homology was found.',
        )

    return resp


@router.post(
    '/restriction_and_ligation',
    response_model=create_model(
        'RestrictionAndLigationResponse',
        sources=(list[RestrictionAndLigationSource], ...),
        sequences=(list[TextFileSequence], ...),
    ),
    summary='Restriction and ligation in a single step. Can also be used for Golden Gate assembly.',
)
async def restriction_and_ligation(
    source: RestrictionAndLigationSource,
    sequences: conlist(TextFileSequence, min_length=1),
    allow_partial_overlap: bool = Query(False, description='Allow for partially overlapping sticky ends.'),
    circular_only: bool = Query(False, description='Only return circular assemblies.'),
):

    fragments = [read_dsrecord_from_json(seq) for seq in sequences]
    invalid_enzymes = get_invalid_enzyme_names(source.restriction_enzymes)
    if len(invalid_enzymes):
        raise HTTPException(404, 'These enzymes do not exist: ' + ', '.join(invalid_enzymes))
    enzymes = RestrictionBatch(first=[e for e in source.restriction_enzymes if e is not None])

    # Lambda function for code clarity
    def create_source(a, is_circular):
        return RestrictionAndLigationSource.from_assembly(
            assembly=a,
            circular=is_circular,
            id=source.id,
            restriction_enzymes=source.restriction_enzymes,
            fragments=fragments,
        )

    # Algorithm used by assembly class
    def algo(x, y, _l):
        # By default, we allow blunt ends
        return restriction_ligation_overlap(x, y, enzymes, allow_partial_overlap, True)

    resp = generate_assemblies(source, create_source, fragments, circular_only, algo, True)

    if len(resp['sources']) == 0:
        raise HTTPException(400, 'No compatible restriction-ligation was found.')

    return resp


@router.post(
    '/gateway',
    response_model=create_model(
        'GatewayResponse', sources=(list[GatewaySource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def gateway(
    source: GatewaySource,
    sequences: conlist(TextFileSequence, min_length=1),
    circular_only: bool = Query(False, description='Only return circular assemblies.'),
    only_multi_site: bool = Query(
        False, description='Only return assemblies where more than one site per sequence recombined.'
    ),
):

    fragments = [read_dsrecord_from_json(seq) for seq in sequences]
    greedy = source.greedy

    # Lambda function for code clarity
    def create_source(a, is_circular):
        return GatewaySource.from_assembly(
            assembly=a,
            circular=is_circular,
            id=source.id,
            reaction_type=source.reaction_type,
            fragments=fragments,
        )

    # Algorithm used by assembly class
    def algo(x, y, _l):
        # By default, we allow blunt ends
        return gateway_overlap(x, y, source.reaction_type, greedy)

    def annotate(x):
        return annotate_gateway_sites(x, greedy)

    resp = generate_assemblies(source, create_source, fragments, circular_only, algo, False, product_callback=annotate)

    if len(resp['sources']) == 0:
        # Build a list of all the sites in the fragments
        sites_in_fragments = list()
        for frag in fragments:
            sites_in_fragments.append(list(find_gateway_sites(frag, greedy).keys()))
        formatted_strings = [f'fragment {i + 1}: {", ".join(sites)}' for i, sites in enumerate(sites_in_fragments)]
        raise HTTPException(
            400,
            f'Inputs are not compatible for {source.reaction_type} reaction.\n\n' + '\n'.join(formatted_strings),
        )

    if only_multi_site:
        multi_site_sources = [
            i
            for i, s in enumerate(resp['sources'])
            if all(join.left_location != join.right_location for join in s.assembly)
        ]
        sources = [resp['sources'][i] for i in multi_site_sources]
        sequences = [resp['sequences'][i] for i in multi_site_sources]
        return {'sources': sources, 'sequences': sequences}

    return resp


@router.post(
    '/cre_lox_recombination',
    response_model=create_model(
        'CreLoxRecombinationResponse',
        sources=(list[CreLoxRecombinationSource], ...),
        sequences=(list[TextFileSequence], ...),
    ),
)
async def cre_lox_recombination(source: CreLoxRecombinationSource, sequences: conlist(TextFileSequence, min_length=1)):
    fragments = [read_dsrecord_from_json(seq) for seq in sequences]

    # Lambda function for code clarity
    def create_source(a, is_circular):
        return CreLoxRecombinationSource.from_assembly(
            assembly=a, circular=is_circular, id=source.id, fragments=fragments
        )

    resp = generate_assemblies(
        source,
        create_source,
        fragments,
        False,
        cre_loxP_overlap,
        True,
        recombination_mode=True,
        product_callback=annotate_loxP_sites,
    )

    if len(resp['sources']) == 0:
        raise HTTPException(400, 'No compatible Cre/Lox recombination was found.')

    return resp
