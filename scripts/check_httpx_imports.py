#!/usr/bin/env python3
import ast
import sys
from pathlib import Path


def check_file(file_path):
    """Check if httpx is imported in the file and if it's allowed."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    try:
        tree = ast.parse(content)
    except SyntaxError:
        print(f"Error: Could not parse {file_path}")
        return False

    # Check for httpx imports
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for name in node.names:
                if name.name == 'httpx':
                    # Only allow imports in httpClient.py
                    if file_path.name != 'httpClient.py':
                        print(f'Error: httpx import found in {file_path}')
                        print('httpx should only be imported in httpClient.py')
                        return False
    return True


def main():
    files = sys.argv[1:]
    if not files:
        print('No files provided')
        sys.exit(1)

    failed = False
    for file_path in files:
        if not check_file(Path(file_path)):
            failed = True

    if failed:
        sys.exit(1)


if __name__ == '__main__':
    main()
