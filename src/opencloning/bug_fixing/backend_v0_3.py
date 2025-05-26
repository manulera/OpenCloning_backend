from ..pydantic_models import (
    BaseCloningStrategy as CloningStrategy,
    AssemblySource,
    TextFileSequence,
    PrimerModel,
)
from .._version import __version__
import json
from opencloning_linkml.migrations import migrate
import os
from packaging import version
import shutil


def fix_backend_v0_3(data: dict) -> dict | None:
    """
    Fix the bug in `SimpleSequenceLocation.from_simple_location`, where origin-spanning
    features were not being read correctly and turned into the entire sequence. This was being used in
    `generate_assemblies` and producing wrong assembly products.

    ## Error fixing
    There were two errors, fixed in https://github.com/manulera/OpenCloning_backend/pull/305

    The first error concerned gateway assemblies. `gateway_overlap` was returning the entire
    overlap, which matched regex like twtGTACAAAaaa (for attB1). That created assemblies in which
    the overlapping part may have mismatches on the w. Now, instead of returning the whole twtGTACAAAaaa
    as overlap, it returns only the common part GTACAAA.

    The second error was a bug in `SimpleSequenceLocation.from_simple_location`, where origin-spanning
    features were not being read correctly and turned into the entire sequence. This was being used in
    `generate_assemblies` and producing wrong assembly products. There is an example in the pull request.
    """

    # Make sure that it is a valid CloningStrategy
    cs = CloningStrategy.model_validate(data)

    # First fix gateway assemblies
    problematic_source_ids = set()

    for source in data['sources']:
        if source['type'] == 'GatewaySource':
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

    return data


def main(file_path: str):
    file_dir = os.path.dirname(file_path)
    file_base = os.path.splitext(os.path.basename(file_path))[0]
    backup_file_path = os.path.join(file_dir, f'{file_base}_backup.json')
    new_file_path = os.path.join(file_dir, f'{file_base}_needs_fixing.json')

    with open(file_path, 'r') as f:
        data = json.load(f)
    if 'backend_version' not in data or data['backend_version'] is None:
        old_version = data['schema_version'] if 'schema_version' in data else None
        migrated_data = migrate(data)
        if migrated_data['schema_version'] != old_version:
            # If the data was migrated, create a backup file
            shutil.copy(file_path, backup_file_path)
            # Overwrite the original file with the migrated data
            with open(file_path, 'w') as f:
                f.write(json.dumps(migrated_data, indent=2))

        # Fix the data
        new_data = fix_backend_v0_3(migrated_data)
        if new_data is not None:
            cs = CloningStrategy.model_validate(new_data)
            cs.backend_version = __version__ if version.parse(__version__) > version.parse('0.3') else '0.3'
            with open(new_file_path, 'w') as f:
                f.write(cs.model_dump_json(indent=2))


if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1:
        print('Usage: python assembly_features_spanning_origin.py <file1> <file2> ...')
        sys.exit(1)

    file_paths = sys.argv[1:]

    for file_path in file_paths:
        if file_path.endswith('_needs_fixing.json') or file_path.endswith('_backup.json'):
            print(f'Skipping {file_path}')
            continue
        main(file_path)
