import re
import io

from bs4 import BeautifulSoup
import httpx

from opencloning.dna_functions import custom_file_parser
from opencloning_linkml.datamodel import SequenceFileFormat


URL = 'https://goldenbraidpro.com/do/domestication/'
SEQ_FILE_PATH = 'ase1_cerevisiae.fasta'


def extract_csrf_token(html: str) -> str | None:
    match = re.search(r'name="csrfmiddlewaretoken"\s+value="([^"]+)"', html)
    return match.group(1) if match else None


def main() -> None:
    # Keep one session so cookies + CSRF flow are preserved.
    with httpx.Client(timeout=30.0) as client:
        first_get = client.get(URL, follow_redirects=True)
        csrf_token = extract_csrf_token(first_get.text)
        if csrf_token is None:
            raise RuntimeError('Could not extract csrfmiddlewaretoken from domestication page')

        print(f"GET status: {first_get.status_code}")
        print('Extracted CSRF from page: True')
        print(f"CSRF token starts with: {csrf_token[:12]}...")
        print(f"Cookies after GET: {dict(client.cookies)}")
        print('-' * 60)

        with open(SEQ_FILE_PATH, 'rb') as f:
            seq_bytes = f.read()
        print(f"Loaded sequence file: {SEQ_FILE_PATH} ({len(seq_bytes)} bytes)")

        response = client.post(
            URL,
            data={
                'csrfmiddlewaretoken': csrf_token,
                'category': 'CDS (B3-B4-B5)',
                'prefix': '',
                'suffix': '',
                'enzymes': ['BsmBI', 'BsaI'],
            },
            files={'seq': (SEQ_FILE_PATH, seq_bytes, 'application/octet-stream')},
            headers={'Referer': URL},
            follow_redirects=True,
        )
        print(f"POST status: {response.status_code}")
        print(f"Content-Type: {response.headers.get('content-type', '')}")

        soup = BeautifulSoup(response.text, 'html.parser')
        record_input = soup.find('input', attrs={'name': 'record'})
        if record_input is None or record_input.get('value') is None:
            raise RuntimeError('Could not find GenBank record input in response HTML')

        dseqrecord = custom_file_parser(io.StringIO(record_input['value']), SequenceFileFormat('genbank'))[0]
        print(f"Parsed Dseqrecord name: {dseqrecord.name}")
        print(f"Parsed Dseqrecord length: {len(dseqrecord)}")
        print(f"Parsed Dseqrecord circular: {dseqrecord.circular}")


if __name__ == '__main__':
    main()
