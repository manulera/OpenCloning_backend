genome_region_examples = {
    'full': {
        'summary': 'All parameters provided',
        'value': {
            'id': 1,
            'sequence_accession': 'NC_003424.3',
            'assembly_accession': 'GCF_000002945.2',
            'locus_tag': 'SPOM_SPAPB1A10.09',
            'gene_id': 2543372,
            'start': 1877009,
            'end': 1881726,
            'strand': 1,
        },
    },
    'full_with_genbank_accession': {
        'summary': 'All parameters provided, but sequence accession is GenBank',
        'value': {
            'id': 1,
            'sequence_accession': 'CU329670.1',
            'assembly_accession': 'GCF_000002945.2',
            'locus_tag': 'SPOM_SPAPB1A10.09',
            'gene_id': 2543372,
            'start': 1877009,
            'end': 1881726,
            'strand': 1,
        },
    },
    'id_omitted': {
        'summary': 'Gene ID omitted (filled in response)',
        'value': {
            'id': 1,
            'sequence_accession': 'NC_003424.3',
            'assembly_accession': 'GCF_000002945.2',
            'locus_tag': 'SPOM_SPAPB1A10.09',
            'start': 1877009,
            'end': 1881726,
            'strand': 1,
        },
    },
    'assembly_accession_omitted': {
        'summary': 'Sequence accession only',
        'value': {
            'id': 1,
            'sequence_accession': 'NC_003424.3',
            'start': 1877009,
            'end': 1881726,
            'strand': 1,
        },
    },
    'viral_sequence': {
        'summary': 'Viral sequence not associated with assembly',
        'value': {
            'id': 1,
            'sequence_accession': 'DQ208311.2',
            'start': 20,
            'end': 2050,
            'strand': -1,
        },
    },
}

oligonucleotide_hybridization_examples = {
    'default': {
        'summary': 'Typical example',
        'description': 'blah',
        'value': {
            'source': {
                'id': 1,
                'input': [
                    {'sequence': 2},
                    {'sequence': 3},
                ],
            },
            'primers': [
                {'id': 2, 'name': 'primer1', 'sequence': 'aaGCGGCCGCgtagaactttatgtgcttccttacattggt'},
                {'id': 3, 'name': 'primer2', 'sequence': 'aaGCGGCCGCaccaatgtaaggaagcacataaagttctac'},
            ],
        },
    },
}

benchling_url_examples = {
    'default': {
        'summary': 'Typical example',
        'value': {
            'id': 0,
            'repository_name': 'benchling',
            'repository_id': 'https://benchling.com/siverson/f/lib_B94YxDHhQh-cidar-moclo-library/seq_kryGidaz-c0062_cd.gb',
        },
    },
}

snapgene_plasmid_examples = {
    'default': {
        'summary': 'Typical example',
        'value': {
            'id': 0,
            'repository_name': 'snapgene',
            'repository_id': 'basic_cloning_vectors/pEASY-T1_(linearized)',
        },
    },
}
