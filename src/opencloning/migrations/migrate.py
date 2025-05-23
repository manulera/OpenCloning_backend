from .__init__ import migrate
from opencloning_linkml.migrations import migrate as migrate_model

import os

if __name__ == '__main__':
    import json
    import argparse
    import shutil

    parser = argparse.ArgumentParser(description='Migrate JSON files to a target schema version.')
    parser.add_argument('input_files', nargs='+', help='Input JSON files to migrate')
    parser.add_argument('--target-version', default=None, help='Target schema version (optional)')

    args = parser.parse_args()

    for input_file in args.input_files:

        # Create a backup of the original file
        name, extension = os.path.splitext(input_file)
        input_file_old = name + '_old' + extension
        shutil.copy2(input_file, input_file_old)

        with open(input_file, 'r') as f:
            data = json.load(f)

        migrated_model = migrate_model(data)

        migrated_data = migrate(migrated_model, args.target_version)

        # Write migrated data back to the same file
        with open(input_file, 'w') as f:
            json.dump(migrated_data, f, indent=2)

        print(f"Migrated {input_file} to version {migrated_data['schema_version']}")
        print(f"Original file backed up as {input_file_old}")
