from httpx import AsyncClient
import asyncio


# Update seva catalog
async def update_seva_catalog(path: str):
    async with AsyncClient() as client:
        response = await client.get(
            'https://raw.githubusercontent.com/manulera/seva_plasmids_index/refs/heads/master/index.json'
        )
        data = response.json()
        with open(path, 'w') as f:
            for item in data:
                if item['id'] == '':
                    continue
                name = item['Name']
                genbank_link = (
                    item['GenBank_link']
                    .replace('http://www.ncbi.nlm.nih.gov/nuccore/', '')
                    .replace('https://www.ncbi.nlm.nih.gov/nuccore/', '')
                )
                if genbank_link == '':
                    genbank_link = item['GenBank']
                    print(f'Replacing GenBank link for {name} with {genbank_link}')
                f.write(f'{name}\t{genbank_link}\n')


async def update_snapgene_catalog(path: str):
    async with AsyncClient() as client:
        response = await client.get(
            'https://raw.githubusercontent.com/manulera/SnapGene_crawler/refs/heads/master/index.yaml'
        )
        with open(path, 'w') as f:
            f.write(response.text)


asyncio.run(update_seva_catalog('src/opencloning/catalogs/seva.tsv'))
asyncio.run(update_snapgene_catalog('src/opencloning/catalogs/snapgene.yaml'))
