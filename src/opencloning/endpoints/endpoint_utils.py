from fastapi import HTTPException
from pydna.dseqrecord import Dseqrecord
from opencloning_linkml.datamodel import Source, TextFileSequence
from typing import Literal
from opencloning.dna_functions import format_sequence_genbank
from pydna.opencloning_models import id_mode


def format_products(
    products: list[Dseqrecord],
    chosen_source: Source | None,
    output_name: str,
    no_products_error_message: str = 'No products were found.',
) -> dict[Literal['sources', 'sequences'], list[Source] | list[TextFileSequence]]:

    if chosen_source is not None:
        this_source_dict = chosen_source.model_dump()
        for prod in products:
            with id_mode(use_python_internal_id=False):
                if prod.source.to_pydantic_model(0).model_dump() == this_source_dict:
                    return {
                        'sources': [this_source_dict],
                        'sequences': [format_sequence_genbank(prod, chosen_source.output_name)],
                    }
        raise HTTPException(400, 'The provided assembly is not valid.')

    if len(products) == 0:
        raise HTTPException(400, no_products_error_message)

    with id_mode(use_python_internal_id=False):
        return {
            'sources': [p.source.to_pydantic_model(0).model_dump() for p in products],
            'sequences': [format_sequence_genbank(p, output_name) for p in products],
        }
