from fastapi import Query, HTTPException, Response
from Bio.Restriction.Restriction_Dictionary import rest_dict
from pydantic import ValidationError
from opencloning_linkml.migrations import migrate
from opencloning_linkml._version import __version__ as schema_version

from ..bug_fixing.backend_v0_3 import fix_backend_v0_3

from ..dna_functions import (
    format_sequence_genbank,
    read_dsrecord_from_json,
)
from ..dna_utils import align_sanger_traces
from ..pydantic_models import (
    TextFileSequence,
    BaseCloningStrategy,
)
from ..get_router import get_router
from .._version import __version__ as backend_version


router = get_router()


@router.get('/version')
async def get_version():
    return {'backend_version': backend_version, 'schema_version': schema_version}


@router.get('/restriction_enzyme_list', response_model=dict[str, list[str]])
async def get_restriction_enzyme_list():
    """Return the dictionary of restriction enzymes"""
    return {'enzyme_names': list(rest_dict.keys())}


@router.post(
    '/validate',
    summary='Validate a cloning strategy',
    responses={
        200: {
            'description': 'The cloning strategy is valid',
            'headers': {
                'x-warning': {
                    'description': 'A warning returned if the file either contains errors or is in a previous version of the model',
                    'schema': {'type': 'string'},
                },
            },
        },
        422: {
            'description': 'The cloning strategy is invalid',
        },
    },
)
async def cloning_strategy_is_valid(data: dict, response: Response):
    """Validate a cloning strategy"""
    warnings = []
    if any(key not in data for key in ['primers', 'sources', 'sequences']):
        raise HTTPException(status_code=422, detail='The cloning strategy is invalid')

    try:
        migrated_data = migrate(data)
        if migrated_data is None:
            BaseCloningStrategy.model_validate(data)
            return None

        data = migrated_data
        warnings.append(
            'The cloning strategy is in a previous version of the model and has been migrated to the latest version.'
        )

        fixed_data = fix_backend_v0_3(data)
        if fixed_data is not None:
            data = fixed_data
            warnings.append('The cloning strategy contained an error and has been turned into a template.')
        cs = BaseCloningStrategy.model_validate(data)
        if len(warnings) > 0:
            response.headers['x-warning'] = ';'.join(warnings)
            return cs
        return None

    except ValidationError:
        raise HTTPException(status_code=422, detail='The cloning strategy is invalid')


@router.post('/rename_sequence', response_model=TextFileSequence)
async def rename_sequence(
    sequence: TextFileSequence,
    name: str = Query(..., description='The new name of the sequence.', pattern=r'^[^\s]+$'),
):
    """Rename a sequence"""
    dseqr = read_dsrecord_from_json(sequence)
    return format_sequence_genbank(dseqr, name)


@router.post('/align_sanger', response_model=list[str])
async def align_sanger(
    sequence: TextFileSequence,
    traces: list[str],
):
    """Align a list of sanger traces to a sequence"""

    dseqr = read_dsrecord_from_json(sequence)
    try:
        return align_sanger_traces(dseqr, traces)
    except RuntimeError as e:
        raise HTTPException(status_code=400, detail=str(e))
