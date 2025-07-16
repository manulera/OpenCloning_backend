from fastapi import Query, HTTPException
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydantic import create_model, Field
from typing import Annotated

from ..dna_functions import (
    format_sequence_genbank,
    oligonucleotide_hybridization_overhangs,
)
from ..pydantic_models import (
    PrimerModel,
    TextFileSequence,
    ManuallyTypedSource,
    OligoHybridizationSource,
    SourceInput,
)

from .. import request_examples
from ..get_router import get_router

router = get_router()


@router.post(
    '/manually_typed',
    response_model=create_model(
        'ManuallyTypedResponse', sources=(list[ManuallyTypedSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def manually_typed(source: ManuallyTypedSource):
    """Return the sequence from a manually typed sequence"""
    if source.circular:
        seq = Dseqrecord(source.user_input, circular=source.circular)
    else:
        seq = Dseqrecord(
            Dseq.from_full_sequence_and_overhangs(
                source.user_input, source.overhang_crick_3prime, source.overhang_watson_3prime
            )
        )
    return {'sequences': [format_sequence_genbank(seq, source.output_name)], 'sources': [source]}


@router.post(
    '/oligonucleotide_hybridization',
    response_model=create_model(
        'OligoHybridizationResponse',
        sources=(list[OligoHybridizationSource], ...),
        sequences=(list[TextFileSequence], ...),
    ),
    openapi_extra={
        'requestBody': {
            'content': {'application/json': {'examples': request_examples.oligonucleotide_hybridization_examples}}
        }
    },
)
async def oligonucleotide_hybridization(
    source: OligoHybridizationSource,
    primers: Annotated[list[PrimerModel], Field(min_length=1, max_length=2)],
    minimal_annealing: int = Query(20, description='The minimal annealing length for each primer.'),
):
    if len(source.input):
        watson_seq = next((p.sequence for p in primers if p.id == source.input[0].sequence), None)
        crick_seq = next((p.sequence for p in primers if p.id == source.input[1].sequence), None)
    else:
        watson_seq = primers[0].sequence
        crick_seq = primers[1].sequence if len(primers) > 1 else watson_seq
        source.input = [SourceInput(sequence=primers[0].id), SourceInput(sequence=primers[1].id)]

    if watson_seq is None or crick_seq is None:
        raise HTTPException(404, 'Invalid oligo id.')

    # The overhang is provided
    if source.overhang_crick_3prime is not None:
        ovhg_watson = len(watson_seq) - len(crick_seq) + source.overhang_crick_3prime
        minimal_annealing = len(watson_seq)
        if source.overhang_crick_3prime < 0:
            minimal_annealing += source.overhang_crick_3prime
        if ovhg_watson > 0:
            minimal_annealing -= ovhg_watson

    try:
        possible_overhangs = oligonucleotide_hybridization_overhangs(watson_seq, crick_seq, minimal_annealing)
    except ValueError as e:
        raise HTTPException(400, *e.args)

    if len(possible_overhangs) == 0:
        raise HTTPException(400, 'No pair of annealing oligos was found. Try changing the annealing settings.')

    if source.overhang_crick_3prime is not None:
        if source.overhang_crick_3prime not in possible_overhangs:
            raise HTTPException(400, 'The provided overhang is not compatible with the primers.')

        return {
            'sources': [source],
            'sequences': [
                format_sequence_genbank(
                    Dseqrecord(Dseq(watson_seq, crick_seq, source.overhang_crick_3prime)), source.output_name
                )
            ],
        }

    out_sources = list()
    out_sequences = list()
    for overhang in possible_overhangs:
        new_source = source.model_copy()
        new_source.overhang_crick_3prime = overhang
        out_sources.append(new_source)
        out_sequences.append(
            format_sequence_genbank(
                Dseqrecord(Dseq(watson_seq, crick_seq, new_source.overhang_crick_3prime)), source.output_name
            )
        )

    return {'sources': out_sources, 'sequences': out_sequences}
