"""
See info in README.md
"""

from ..pydantic_models import (
    BaseCloningStrategy as CloningStrategy,
    AssemblySource,
    TextFileSequence,
    PrimerModel,
    SequenceLocationStr,
)
from .._version import __version__
import json
import os
from packaging import version
import copy


def fix_backend_v0_3(input_data: dict) -> CloningStrategy | None:

    data = copy.deepcopy(input_data)
    # Make sure that it is a valid CloningStrategy
    cs = CloningStrategy.model_validate(data)

    # First fix gateway assemblies
    problematic_source_ids = set()

    for source in data['sources']:
        if source['type'] == 'GatewaySource':
            # Take the first assembly value and check that the length of features is 7
            assembly = source['assembly']
            if len(assembly):
                feat2check = (
                    assembly[0]['left_location']
                    if assembly[0]['left_location'] is not None
                    else assembly[0]['right_location']
                )
                if len(SequenceLocationStr(feat2check).to_biopython_location()) != 7:
                    problematic_source_ids.add(source['id'])

        elif 'assembly' in source:
            assembly_source = AssemblySource(
                id=source['id'],
                input=source['input'],
                output=source['output'],
                circular=source['circular'],
                assembly=source['assembly'],
            )
            input_seqs = [
                TextFileSequence.model_validate(s) for s in data['sequences'] if s['id'] in assembly_source.input
            ]
            # Sort input_seqs as in input
            input_seqs.sort(key=lambda x: assembly_source.input.index(x.id))
            if source['type'] == 'PCRSource':
                primer_ids = [assembly_source.assembly[0].sequence, assembly_source.assembly[2].sequence]
                primers = [PrimerModel.model_validate(p) for p in data['primers'] if p['id'] in primer_ids]
                input_seqs = [primers[0], input_seqs[0], primers[1]]

            assembly_plan = assembly_source.get_assembly_plan(input_seqs)
            for join in assembly_plan:
                if len(join[2]) != len(join[3]):
                    problematic_source_ids.add(source['id'])
                    break

    if len(problematic_source_ids) == 0:
        return None

    # Replace problematic sources and their output sequences by templates
    problematic_source_ids.update(sum([cs.all_children_source_ids(s) for s in problematic_source_ids], []))
    for source_id in problematic_source_ids:
        source = next(s for s in data['sources'] if s['id'] == source_id)
        output_seq = next(s for s in data['sequences'] if s['id'] == source['output'])
        remove_keys = ['assembly', 'circular']
        source_keep = {key: value for key, value in source.items() if key not in remove_keys}
        source.clear()
        source.update(source_keep)

        seq_keep = {'id': output_seq['id'], 'type': 'TemplateSequence'}
        output_seq.clear()
        output_seq.update(seq_keep)

    return CloningStrategy.model_validate(data)


def main(file_path: str):
    file_dir = os.path.dirname(file_path)
    file_base = os.path.splitext(os.path.basename(file_path))[0]
    new_file_path = os.path.join(file_dir, f'{file_base}_needs_fixing.json')

    with open(file_path, 'r') as f:
        data = json.load(f)

    if 'backend_version' not in data or data['backend_version'] is None:

        # Fix the data
        cs = fix_backend_v0_3(data)

        if cs is not None:
            cs.backend_version = __version__ if version.parse(__version__) > version.parse('0.3') else '0.3'
            with open(new_file_path, 'w') as f:
                f.write(cs.model_dump_json(indent=2, exclude_none=True))


if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1:
        print('Usage: python assembly_features_spanning_origin.py <file1> <file2> ...')
        sys.exit(1)

    file_paths = sys.argv[1:]

    for file_path in file_paths:
        if file_path.endswith('_needs_fixing.json'):
            print(f'Skipping {file_path}')
            continue
        main(file_path)
