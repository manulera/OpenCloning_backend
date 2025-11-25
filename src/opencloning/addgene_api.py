from .app_settings import settings
from .http_client import get_http_client
from pydna.dseqrecord import Dseqrecord
from .dna_functions import request_from_addgene
from opencloning_linkml.datamodel import AddgeneIdSource


async def get_plasmid_catalog(name: str) -> str:
    async with get_http_client() as client:
        return client.get(
            'https://api.developers.addgene.org/catalog/plasmid/',
            headers={'Authorization': f'Token {settings.ADDGENE_TOKEN}'},
            params={'name': name},
        )


async def get_plasmid(source: AddgeneIdSource) -> tuple[Dseqrecord, AddgeneIdSource]:
    # TODO replace with a call to the addgene api
    # async with get_http_client() as client:
    #     resp = await client.get(
    #         f'https://api.developers.addgene.org/catalog/plasmid-with-sequences/{plasmid_id}/',
    #         headers={'Authorization': f'Token {settings.ADDGENE_TOKEN}'},
    #     )
    return await request_from_addgene(source)


get_plasmid_catalog('pfa6a gfp kan')
