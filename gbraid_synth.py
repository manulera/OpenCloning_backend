import re
import io

from Bio.SeqFeature import SeqFeature
from Bio.Restriction.Restriction import RestrictionBatch
from bs4 import BeautifulSoup
import httpx
from pydna.dseqrecord import Dseqrecord
from pydna.assembly2 import pcr_assembly, restriction_ligation_assembly
from pydna.parsers import parse as pydna_parse
from pydna.primer import Primer
from pydna.opencloning_models import AddgeneIdSource, CloningStrategy, Source, SourceInput
from pydna.sequence_regex import dseqrecord_finditer
from pydna.utils import shift_location, SimpleLocation, Location

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
CATEGORY = 'CDS (B3-B4-B5)'
URL = f'https://goldenbraidpro.com/do/{CLONING_TYPE}/'
SEQ_FILE_PATH = 'ase1_cerevisiae.fasta'
NAME = 'ase1_cerevisiae'


def get_location_of_sequence(dseqrecord: Dseqrecord, seq: str) -> Location | None:
    results = list(dseqrecord_finditer(seq, dseqrecord))
    if len(results) != 1:
        return None
    match = results[0]
    loc = SimpleLocation(match.start(), match.end(), 1)
    return shift_location(loc, 0, len(dseqrecord))


def add_seqfeature(dseqrecord: Dseqrecord, expected_product_part_seq: str | None, CATEGORY: str, NAME: str) -> None:
    feat_type = 'CDS' if 'CDS' in CATEGORY.upper() else 'misc_feature'
    if expected_product_part_seq is not None:
        loc = get_location_of_sequence(dseqrecord, re.compile(expected_product_part_seq, re.IGNORECASE))
        if loc is not None:
            dseqrecord.features.append(SeqFeature(type=feat_type, location=loc, qualifiers={'label': [NAME]}))


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
    template = pydna_parse(SEQ_FILE_PATH)[0]
    with httpx.Client(timeout=30.0) as client:
        first_get = client.get(URL, follow_redirects=True)
        csrf_token = extract_csrf_token(first_get.text)
        if csrf_token is None:
            raise RuntimeError('Could not extract csrfmiddlewaretoken from domestication page')

        with open(SEQ_FILE_PATH, 'rb') as f:
            seq_bytes = f.read()

        response = client.post(
            URL,
            data={
                'csrfmiddlewaretoken': csrf_token,
                'category': CATEGORY,
                'prefix': '',
                'suffix': '',
                'enzymes': ['BsmBI', 'BsaI'],
            },
            files={'seq': (SEQ_FILE_PATH, seq_bytes, 'application/octet-stream')},
            headers={'Referer': URL},
            follow_redirects=True,
        )
        soup = BeautifulSoup(response.text, 'html.parser')
        response_errors = extract_unique_error_messages(soup)
        if response_errors:
            raise ValueError('; '.join(response_errors))

        record_input = soup.find('input', attrs={'name': 'record'})
        if record_input is None or record_input.get('value') is None:
            raise RuntimeError('Could not find GenBank record input in response HTML')

        expected_product = custom_file_parser(io.StringIO(record_input['value']), SequenceFileFormat('genbank'))[0]
        expected_product.seq.circular = True
        expected_product_part_seq = None
        if len(expected_product.features) == 1:
            expected_product_part_seq = str(expected_product.features[0].extract(expected_product).seq)

        assembly_inputs = []

        if CLONING_TYPE == 'synthesis':
            main_div = soup.find('div', attrs={'id': 'main'})
            if main_div is None:
                raise RuntimeError('Could not find div#main in synthesis response HTML')

            sequence_pre = main_div.find('pre')
            if sequence_pre is None:
                raise RuntimeError('Could not find synthesis sequence <pre> in div#main')

            # Sequence is line-wrapped in HTML; normalize it into a contiguous string.
            synthesis_sequence = re.sub(r'\s+', '', sequence_pre.get_text())
            if not synthesis_sequence:
                raise RuntimeError('Synthesis sequence <pre> is empty')
            if not re.match(r'^[ACGT]+$', synthesis_sequence):
                raise ValueError(f'Synthesis sequence {synthesis_sequence} contains invalid characters')

            synthesized_sequence = Dseqrecord(synthesis_sequence, circular=False)
            synthesized_sequence.source = Source(input=[SourceInput(sequence=template)])
            add_seqfeature(synthesized_sequence, expected_product_part_seq, CATEGORY, NAME)
            assembly_inputs.append(synthesized_sequence)

        elif CLONING_TYPE == 'domestication':

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
            pair_pattern = re.compile(
                r'Oligo forward:\s*([ACGTacgt]+)\s*\n\s*Oligo reverse:\s*([ACGTacgt]+)',
                flags=re.MULTILINE,
            )
            oligo_pairs: list[tuple[Primer, Primer]] = [
                (Primer(match.group(1), name=f'fwd_{i}'), Primer(match.group(2), name=f'rev_{i}'))
                for i, match in enumerate(pair_pattern.finditer(protocol_response.text))
            ]

            for i, (fwd, rev) in enumerate(oligo_pairs, start=1):
                try:
                    product = pcr_assembly(template, fwd, rev, limit=14, mismatches=0)[0]
                except IndexError as exc:
                    raise RuntimeError(f'No PCR product found for oligo pair {i}') from exc
                product.name = f'pcr_product_{i}'
                assembly_inputs.append(product)

        enzymes = RestrictionBatch(first=['BsmBI'])
        assemblies = restriction_ligation_assembly(assembly_inputs + [pUPD2], enzymes, circular_only=True)
        assert len(assemblies) == 1, f"Expected exactly one assembly, got {len(assemblies)}"
        assembly = assemblies[0]
        assert (
            expected_product.seguid() == assembly.seguid()
        ), f"Expected product {expected_product.seguid()} != assembly {assembly.seguid()}"
        assembly.name = 'domesticated_part'
        if CLONING_TYPE == 'domestication':
            add_seqfeature(assembly, expected_product_part_seq, CATEGORY, NAME)
        cs = CloningStrategy.from_dseqrecords([assembly])

        with open('cloning_strategy.json', 'w', encoding='utf-8') as f:
            f.write(cs.model_dump_json(indent=2))


if __name__ == '__main__':
    main()
