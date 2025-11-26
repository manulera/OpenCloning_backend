import os
import yaml


def get_seva_catalog():
    seva_catalog_path = os.path.join(os.path.dirname(__file__), 'seva.tsv')
    seva_catalog = dict()
    with open(seva_catalog_path, 'r') as f:
        for line in f:
            name, genbank_link = line.strip().split('\t')
            seva_catalog[name] = genbank_link
    return seva_catalog


def get_snapgene_catalog():
    snapgene_catalog_path = os.path.join(os.path.dirname(__file__), 'snapgene.yaml')
    with open(snapgene_catalog_path, 'r') as f:
        return yaml.load(f, Loader=yaml.FullLoader)


seva_catalog = get_seva_catalog()
snapgene_catalog = get_snapgene_catalog()
