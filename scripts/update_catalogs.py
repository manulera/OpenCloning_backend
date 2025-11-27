from httpx import AsyncClient
import asyncio
import yaml


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
        data = yaml.load(response.text, Loader=yaml.FullLoader)
        for key, plasmid_set in data.items():
            data[key] = [plasmid['subpath'] for plasmid in plasmid_set['plasmids']]

        with open(path, 'w') as f:
            yaml.dump(data, f)


async def update_openDNA_collections_catalog(path: str):
    async with AsyncClient() as client:
        response = await client.get(
            'https://raw.githubusercontent.com/manulera/Open-DNA-Collections/refs/heads/requestable-collection/scripts/index.json'
        )
        data = response.json()
        reshaped_data = {}
        for item in data:
            collection_name = item['collection']
            if collection_name not in reshaped_data:
                reshaped_data[collection_name] = []
            reshaped_data[collection_name].append(
                {
                    'id': item['id'],
                    'name': item['plasmid_name'],
                    'path': item['path'],
                }
            )
        with open(path, 'w') as f:
            yaml.dump(reshaped_data, f)


asyncio.run(update_seva_catalog('src/opencloning/catalogs/seva.tsv'))
asyncio.run(update_snapgene_catalog('src/opencloning/catalogs/snapgene.yaml'))
asyncio.run(update_openDNA_collections_catalog('src/opencloning/catalogs/openDNA_collections.yaml'))
