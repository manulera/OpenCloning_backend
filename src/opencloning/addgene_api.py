import asyncio
from .app_settings import settings
from .http_client import get_http_client


async def get_plasmid_list(name: str) -> str:
    async with get_http_client() as client:
        resp = await client.get(
            'https://api.developers.addgene.org/catalog/plasmid/',
            headers={'Authorization': f'Token {settings.ADDGENE_TOKEN}'},
            params={'name': name},
        )
    # with open('plasmid_list.json', 'w') as f:
    #     json.dump(resp.json(), f, indent=4)

    # print(resp.json()['results'][0]['id'])
    plasmid_id = resp.json()['results'][0]['id']
    async with get_http_client() as client:
        resp = await client.get(
            f'https://api.developers.addgene.org/catalog/plasmid-with-sequences/{plasmid_id}/',
            headers={'Authorization': f'Token {settings.ADDGENE_TOKEN}'},
        )
    # with open('plasmid_details.json', 'w') as f:
    #     json.dump(resp.json(), f, indent=4)


if __name__ == '__main__':
    asyncio.run(get_plasmid_list('pfa6a GFP'))
