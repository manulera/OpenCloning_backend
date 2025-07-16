from fastapi import Query, HTTPException
from pydna.dseqrecord import Dseqrecord
from pydantic import create_model, Field
from typing import Annotated
from Bio.Restriction import RestrictionBatch

from ..dna_functions import (
    format_sequence_genbank,
    read_dsrecord_from_json,
    get_invalid_enzyme_names,
)
from ..pydantic_models import (
    RestrictionEnzymeDigestionSource,
    TextFileSequence,
    PolymeraseExtensionSource,
    ReverseComplementSource,
)
from ..get_router import get_router

router = get_router()


@router.post(
    '/restriction',
    response_model=create_model(
        'RestrictionEnzymeDigestionResponse',
        sources=(list[RestrictionEnzymeDigestionSource], ...),
        sequences=(list[TextFileSequence], ...),
    ),
)
async def restriction(
    source: RestrictionEnzymeDigestionSource,
    sequences: Annotated[list[TextFileSequence], Field(min_length=1, max_length=1)],
    restriction_enzymes: Annotated[list[str], Query(default_factory=list)],
):
    # There should be 1 or 2 enzymes in the request if the source does not have cuts
    if source.left_edge is None and source.right_edge is None:
        if len(restriction_enzymes) < 1 or len(restriction_enzymes) > 2:
            raise HTTPException(422, 'There should be 1 or 2 restriction enzymes in the request.')
    else:
        if len(restriction_enzymes) != 0:
            raise HTTPException(422, 'There should be no restriction enzymes in the request if source is populated.')
        restriction_enzymes = source.get_enzymes()

    # TODO: this could be moved to the class
    invalid_enzymes = get_invalid_enzyme_names(restriction_enzymes)
    if len(invalid_enzymes):
        raise HTTPException(404, 'These enzymes do not exist: ' + ', '.join(invalid_enzymes))
    enzymes = RestrictionBatch(first=[e for e in restriction_enzymes if e is not None])

    seqr = read_dsrecord_from_json(sequences[0])
    # TODO: return error if the id of the sequence does not correspond

    cutsites = seqr.seq.get_cutsites(*enzymes)
    cutsite_pairs = seqr.seq.get_cutsite_pairs(cutsites)
    sources = [
        RestrictionEnzymeDigestionSource.from_cutsites(*p, [{'sequence': sequences[0].id}], source.id)
        for p in cutsite_pairs
    ]

    all_enzymes = set(enzyme for s in sources for enzyme in s.get_enzymes())
    enzymes_not_cutting = set(restriction_enzymes) - set(all_enzymes)
    if len(enzymes_not_cutting):
        raise HTTPException(400, 'These enzymes do not cut: ' + ', '.join(enzymes_not_cutting))

    try:
        # If the output is known
        if source.left_edge is not None or source.right_edge is not None:

            for i, s in enumerate(sources):
                if s == source:
                    return {
                        'sequences': [format_sequence_genbank(seqr.apply_cut(*cutsite_pairs[i]), source.output_name)],
                        'sources': [s],
                    }

            raise HTTPException(400, 'Invalid restriction enzyme pair.')

        products = [format_sequence_genbank(seqr.apply_cut(*p), source.output_name) for p in cutsite_pairs]

        return {'sequences': products, 'sources': sources}
    except ValueError as e:
        raise HTTPException(400, str(e))


@router.post(
    '/polymerase_extension',
    response_model=create_model(
        'PolymeraseExtensionResponse',
        sources=(list[PolymeraseExtensionSource], ...),
        sequences=(list[TextFileSequence], ...),
    ),
)
async def polymerase_extension(
    source: PolymeraseExtensionSource,
    sequences: Annotated[list[TextFileSequence], Field(min_length=1, max_length=1)],
):
    """Return the sequence from a polymerase extension reaction"""

    dseq = read_dsrecord_from_json(sequences[0])

    if dseq.circular:
        raise HTTPException(400, 'The sequence must be linear.')

    if dseq.seq.ovhg == dseq.seq.watson_ovhg() == 0:
        raise HTTPException(400, 'The sequence must have an overhang.')

    out_sequence = Dseqrecord(dseq.seq.fill_in(), features=dseq.features)

    return {'sequences': [format_sequence_genbank(out_sequence, source.output_name)], 'sources': [source]}


@router.post(
    '/reverse_complement',
    response_model=create_model(
        'ReverseComplementResponse',
        sources=(list[ReverseComplementSource], ...),
        sequences=(list[TextFileSequence], ...),
    ),
)
async def reverse_complement(
    source: ReverseComplementSource,
    sequences: Annotated[list[TextFileSequence], Field(min_length=1, max_length=1)],
):
    dseq = read_dsrecord_from_json(sequences[0])
    out_sequence = dseq.reverse_complement()
    seq_name = source.output_name if source.output_name is not None else dseq.name + '_rc'
    return {'sequences': [format_sequence_genbank(out_sequence, seq_name)], 'sources': [source]}
