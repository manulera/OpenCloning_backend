import re
import io

from Bio.Restriction.Restriction import RestrictionBatch
from bs4 import BeautifulSoup
import httpx
from pydna.assembly2 import pcr_assembly, restriction_ligation_assembly
from pydna.parsers import parse as pydna_parse
from pydna.primer import Primer
from pydna.opencloning_models import AddgeneIdSource, CloningStrategy

from opencloning.dna_functions import custom_file_parser
from opencloning_linkml.datamodel import SequenceFileFormat

pUPD2 = pydna_parse('pUPD2.gb')[0]
pUPD2.source = AddgeneIdSource(
    repository_id='68161',
    addgene_sequence_type='depositor-full',
    sequence_file_url='https://media.addgene.org/snapgene-media/v3.38.0/sequences/119013/5a5f4064-22b9-4374-aa25-8d4d8867b261/addgene-plasmid-68161-sequence-119013.gbk',
)


CLONING_TYPE = 'synthesis'
CLONING_TYPE = 'domestication'

URL = f'https://goldenbraidpro.com/do/{CLONING_TYPE}/'
SEQ_FILE_PATH = 'ase1_cerevisiae.fasta'
SEQ_FILE_PATH = 'packages/opencloning/tests/test_files/dummy_EcoRI.fasta'


def extract_csrf_token(html: str) -> str | None:
    match = re.search(r'name="csrfmiddlewaretoken"\s+value="([^"]+)"', html)
    return match.group(1) if match else None


def extract_unique_error_messages(soup: BeautifulSoup) -> list[str]:
    unique_errors: list[str] = []
    seen: set[str] = set()
    for error_list in soup.find_all('ul', class_='errorlist'):
        for li in error_list.find_all('li'):
            message = li.get_text(strip=True)
            if message and message not in seen:
                seen.add(message)
                unique_errors.append(message)
    return unique_errors


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
        response_errors = extract_unique_error_messages(soup)
        if response_errors:
            raise ValueError('; '.join(response_errors))

        if CLONING_TYPE == 'synthesis':
            record_input = soup.find('input', attrs={'name': 'record'})
            if record_input is None or record_input.get('value') is None:
                raise RuntimeError('Could not find GenBank record input in response HTML')

            dseqrecord = custom_file_parser(io.StringIO(record_input['value']), SequenceFileFormat('genbank'))[0]
            print(f"Parsed Dseqrecord name: {dseqrecord.name}")
            print(f"Parsed Dseqrecord length: {len(dseqrecord)}")
            print(f"Parsed Dseqrecord circular: {dseqrecord.circular}")
        elif CLONING_TYPE == 'domestication':
            record_input = soup.find('input', attrs={'name': 'record'})
            if record_input is None or record_input.get('value') is None:
                raise RuntimeError('Could not find GenBank record input in response HTML')

            expected_product = custom_file_parser(io.StringIO(record_input['value']), SequenceFileFormat('genbank'))[0]
            expected_product.seq.circular = True

            protocol_form = soup.find('form', attrs={'action': '/do/domestication/protocol/'})
            if protocol_form is None:
                raise RuntimeError('Could not find domestication protocol form in response HTML')

            protocol_payload: dict[str, str] = {}
            for input_tag in protocol_form.find_all('input'):
                name = input_tag.get('name')
                value = input_tag.get('value')
                if name is not None and value is not None:
                    protocol_payload[name] = value

            protocol_url = 'https://goldenbraidpro.com/do/domestication/protocol/'
            protocol_response = client.post(
                protocol_url,
                data=protocol_payload,
                headers={'Referer': URL},
                follow_redirects=True,
            )
            print(f"Protocol POST status: {protocol_response.status_code}")
            print(f"Protocol Content-Type: {protocol_response.headers.get('content-type', '')}")
            protocol_soup = BeautifulSoup(protocol_response.text, 'html.parser')
            protocol_errors = extract_unique_error_messages(protocol_soup)
            if protocol_errors:
                print(f"Protocol errors ({len(protocol_errors)} unique):")
                for err in protocol_errors:
                    print(f"- {err}")

            with open('domestication_protocol_response.html', 'w', encoding='utf-8') as f:
                f.write(protocol_response.text)
            print('Wrote HTML response to domestication_protocol_response.html')

            pair_pattern = re.compile(
                r'Oligo forward:\s*([ACGTacgt]+)\s*\n\s*Oligo reverse:\s*([ACGTacgt]+)',
                flags=re.MULTILINE,
            )
            oligo_pairs: list[tuple[Primer, Primer]] = [
                (Primer(match.group(1), name=f'fwd_{i}'), Primer(match.group(2), name=f'rev_{i}'))
                for i, match in enumerate(pair_pattern.finditer(protocol_response.text))
            ]

            print(f"Parsed oligo pairs: {len(oligo_pairs)}")
            for i, (fwd, rev) in enumerate(oligo_pairs, start=1):
                print(f"Pair {i}: forward={fwd.seq} reverse={rev.seq}")

            template = pydna_parse(SEQ_FILE_PATH)[0]
            print(f"Parsed template sequence: {template.name} (length={len(template)})")

            pcr_products = []
            for i, (fwd, rev) in enumerate(oligo_pairs, start=1):
                try:
                    product = pcr_assembly(template, fwd, rev, limit=14, mismatches=0)[0]
                except IndexError as exc:
                    raise RuntimeError(f'No PCR product found for oligo pair {i}') from exc
                product.name = f'pcr_product_{i}'
                pcr_products.append(product)
                print(f"PCR product {i}: length={len(product)}")

            enzymes = RestrictionBatch(first=['BsmBI'])
            assemblies = restriction_ligation_assembly(pcr_products + [pUPD2], enzymes, circular_only=True)
            assert len(assemblies) == 1, f"Expected exactly one assembly, got {len(assemblies)}"
            assembly = assemblies[0]
            assert (
                expected_product.seguid() == assembly.seguid()
            ), f"Expected product {expected_product.seguid()} != assembly {assembly.seguid()}"
            assembly.name = 'domesticated_part'
            cs = CloningStrategy.from_dseqrecords([assembly])
            print(f"Cloning strategy: {cs.model_dump_json(indent=2)}")
            with open('cloning_strategy.json', 'w', encoding='utf-8') as f:
                f.write(cs.model_dump_json(indent=2))
            print('Wrote cloning strategy to cloning_strategy.json')


if __name__ == '__main__':
    main()
