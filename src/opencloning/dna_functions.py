from functools import cmp_to_key
from urllib.error import HTTPError
from Bio.Restriction.Restriction import RestrictionBatch
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from .pydantic_models import TextFileSequence, AddgeneIdSource, SequenceFileFormat, WekWikGeneIdSource, SEVASource
from opencloning_linkml.datamodel import PlannotateAnnotationReport
from pydna.parsers import parse as pydna_parse
import requests
from bs4 import BeautifulSoup
import regex
from Bio.SeqFeature import SimpleLocation, Location
from pydna.utils import shift_location
from pydna.common_sub_strings import common_sub_strings
from Bio.SeqIO import parse as seqio_parse
import io
import warnings
from Bio.SeqIO.InsdcIO import GenBankIterator, GenBankScanner
import re
import httpx
from .ncbi_requests import get_genbank_sequence


def format_sequence_genbank(seq: Dseqrecord, seq_name: str = None) -> TextFileSequence:

    if seq_name is not None:
        seq.name = seq_name
    elif seq.name.lower() == 'exported':
        correct_name(seq)

    return TextFileSequence(
        id=0,
        file_content=seq.format('genbank'),
        sequence_file_format=SequenceFileFormat('genbank'),
        overhang_crick_3prime=seq.seq.ovhg,
        overhang_watson_3prime=seq.seq.watson_ovhg(),
    )


def read_dsrecord_from_json(seq: TextFileSequence) -> Dseqrecord:
    initial_dseqrecord: Dseqrecord = pydna_parse(seq.file_content)[0]
    if seq.overhang_watson_3prime == 0 and seq.overhang_crick_3prime == 0:
        out_dseq_record = initial_dseqrecord
    else:
        out_dseq_record = Dseqrecord(
            Dseq.from_full_sequence_and_overhangs(
                str(initial_dseqrecord.seq), seq.overhang_crick_3prime, seq.overhang_watson_3prime
            ),
            features=initial_dseqrecord.features,
        )
    # We set the id to the integer converted to integer (this is only
    # useful for assemblies)
    out_dseq_record.id = str(seq.id)
    return out_dseq_record


def get_invalid_enzyme_names(enzyme_names_list: list[str | None]) -> list[str]:
    rest_batch = RestrictionBatch()
    invalid_names = list()
    for name in enzyme_names_list:
        # Empty enzyme names are the natural edges of the molecule
        if name is not None:
            try:
                rest_batch.format(name)
            except ValueError:
                invalid_names.append(name)
    return invalid_names


async def get_sequences_from_file_url(
    url: str, format: SequenceFileFormat = SequenceFileFormat('genbank')
) -> list[Dseqrecord]:
    # TODO once pydna parse is fixed it should handle urls that point to non-gb files
    async with httpx.AsyncClient() as client:
        resp = await client.get(url)

    if resp.status_code != 200:
        raise HTTPError(url, 404, 'file requested from url not found', 'file requested from url not found', None)
    if format == SequenceFileFormat('snapgene'):
        return custom_file_parser(io.BytesIO(resp.content), format)
    else:
        return custom_file_parser(io.StringIO(resp.text), format)


def get_sequence_from_snagene_url(url: str) -> Dseqrecord:
    resp = requests.get(url)
    # Check that resp.content is not empty
    if len(resp.content) == 0:
        raise HTTPError(url, 404, 'invalid snapgene id', 'invalid snapgene id', None)
    parsed_seq = next(seqio_parse(io.BytesIO(resp.content), 'snapgene'))
    circularize = 'topology' in parsed_seq.annotations.keys() and parsed_seq.annotations['topology'] == 'circular'
    return Dseqrecord(parsed_seq, circular=circularize)


async def request_from_addgene(source: AddgeneIdSource) -> tuple[Dseqrecord, AddgeneIdSource]:

    url = f'https://www.addgene.org/{source.repository_id}/sequences/'
    async with httpx.AsyncClient() as client:
        resp = await client.get(url)
    if resp.status_code == 404:
        raise HTTPError(url, 404, 'wrong addgene id', 'wrong addgene id', None)
    soup = BeautifulSoup(resp.content, 'html.parser')

    # Get a span.material-name from the soup, see https://github.com/manulera/OpenCloning_backend/issues/182
    plasmid_name = soup.find('span', class_='material-name').text.replace(' ', '_')

    if source.sequence_file_url:
        dseqr = (await get_sequences_from_file_url(source.sequence_file_url))[0]
        dseqr.name = plasmid_name
        return dseqr, source

    sequence_file_url_dict = dict()
    for _type in ['depositor-full', 'depositor-partial', 'addgene-full', 'addgene-partial']:
        sequence_file_url_dict[_type] = []
        if soup.find(id=_type) is not None:
            sequence_file_url_dict[_type] = [
                a.get('href') for a in soup.find(id=_type).findAll(class_='genbank-file-download')
            ]

    # TODO provide addgene sequencing data supporting the sequence
    # We prefer to return addgene full if both available
    products = list()
    sources = list()
    for _type in ['addgene-full', 'depositor-full']:
        if len(sequence_file_url_dict[_type]) > 0:
            for seq_url in sequence_file_url_dict[_type]:
                new_source = source.model_copy()
                new_source.sequence_file_url = seq_url
                new_source.addgene_sequence_type = _type
                sources.append(new_source)
                # There should be only one sequence
                products.append((await get_sequences_from_file_url(seq_url))[0])

    if len(products) == 0:
        # They may have only partial sequences
        raise HTTPError(
            url,
            404,
            f'The requested plasmid does not have full sequences, see https://www.addgene.org/{source.repository_id}/sequences/',
            f'The requested plasmid does not have full sequences, see https://www.addgene.org/{source.repository_id}/sequences/',
            None,
        )

    # Rename the plasmid
    for p in products:
        p.name = plasmid_name
    return products[0], sources[0]


async def request_from_wekwikgene(source: WekWikGeneIdSource) -> tuple[Dseqrecord, WekWikGeneIdSource]:
    url = f'https://wekwikgene.wllsb.edu.cn/plasmids/{source.repository_id}'
    async with httpx.AsyncClient() as client:
        resp = await client.get(url)
    if resp.status_code == 404:
        raise HTTPError(url, 404, 'invalid wekwikgene id', 'invalid wekwikgene id', None)
    soup = BeautifulSoup(resp.content, 'html.parser')
    # Get the sequence file URL from the page
    sequence_file_url = soup.find('a', text=lambda x: x and 'Download Sequence' in x).get('href')
    sequence_name = soup.find('h1', class_='plasmid__info__name').text.replace(' ', '_')
    seq = (await get_sequences_from_file_url(sequence_file_url, 'snapgene'))[0]
    seq.name = sequence_name
    source.sequence_file_url = sequence_file_url
    return seq, source


async def get_seva_plasmid(source: SEVASource) -> tuple[Dseqrecord, SEVASource]:
    if 'ncbi.nlm.nih.gov/nuccore' in source.sequence_file_url:
        genbank_id = source.sequence_file_url.split('/')[-1]
        seq = await get_genbank_sequence(genbank_id)
        seq.name = source.repository_id
    elif source.sequence_file_url.startswith('https://seva-plasmids.com'):
        seq_list = await get_sequences_from_file_url(source.sequence_file_url)
        if len(seq_list) == 0:
            raise ValueError('No sequences found in SEVA file')
        seq = seq_list[0]
    else:
        raise HTTPError(source.sequence_file_url, 404, 'invalid SEVA url', 'invalid SEVA url', None)
    if not seq.circular:
        seq = seq.looped()
    return seq, source


def correct_name(dseq: Dseqrecord):
    # Can set the name from keyword if locus is set to Exported
    if dseq.name.lower() == 'exported' and dseq.locus.lower() == 'exported' and 'keywords' in dseq.annotations:
        dseq.name = dseq.annotations['keywords'][0].replace(' ', '_')


def location_sorter(x, y) -> int:
    """
    Sort by start, then length, then strand.
    """
    if x.parts[0].start != y.parts[0].start:
        return x.parts[0].start - y.parts[0].start
    elif x.parts[-1].end != y.parts[-1].end:
        return x.parts[-1].end - y.parts[-1].end
    return x.strand - y.strand


def get_all_regex_feature_edges(pattern: str, seq: str, is_circular: bool) -> list[tuple[int, int]]:

    subject = 2 * seq if is_circular else seq

    compiled_pattern = regex.compile(pattern, regex.IGNORECASE)
    compiled_pattern_rev = regex.compile('(?r)' + pattern, regex.IGNORECASE)

    matches = list(regex.finditer(compiled_pattern, subject, overlapped=True))
    matches += list(regex.finditer(compiled_pattern_rev, subject, overlapped=True))

    # In circular objects we remove the matches that span the sequence more than once: m.end() - m.start() <= len(seq)
    return list(set([(m.start(), m.end()) for m in matches if (m.end() - m.start() <= len(seq))]))


def find_sequence_regex(pattern: str, seq: str, is_circular: bool) -> list[Location]:

    feature_locations = list()

    # Strand 1
    feature_edges = get_all_regex_feature_edges(pattern, seq, is_circular)
    # We use shift_location to format origin-spanning features in circular DNA
    feature_locations += [shift_location(SimpleLocation(start, end, 1), 0, len(seq)) for start, end in feature_edges]

    # Strand -1
    feature_edges = get_all_regex_feature_edges(pattern, reverse_complement(seq), is_circular)
    feature_locations += [
        shift_location(SimpleLocation(start, end, 1)._flip(len(seq)), 0, len(seq)) for start, end in feature_edges
    ]

    # We return a unique list, cannot use a set because Location is not hashable
    return sorted(
        [x for i, x in enumerate(feature_locations) if x not in feature_locations[:i]], key=cmp_to_key(location_sorter)
    )


# Could be useful at some point
# def seq_overlap_length(dseq: Dseq) -> int:
#     return len(dseq) - abs(dseq.ovhg) - abs(dseq.watson_ovhg())


def oligonucleotide_hybridization_overhangs(
    fwd_oligo_seq: str, rvs_oligo_seq: str, minimal_annealing: int
) -> list[int]:
    """
    Returns possible overhangs between two oligos, and returns an error if mismatches are found.

    see https://github.com/manulera/OpenCloning_backend/issues/302 for notation

    """
    matches = common_sub_strings(fwd_oligo_seq.lower(), reverse_complement(rvs_oligo_seq.lower()), minimal_annealing)

    for pos_fwd, pos_rvs, length in matches:

        if (pos_fwd != 0 and pos_rvs != 0) or (
            pos_fwd + length < len(fwd_oligo_seq) and pos_rvs + length < len(rvs_oligo_seq)
        ):
            raise ValueError('The oligonucleotides can anneal with mismatches')

    # Return possible overhangs
    return [pos_rvs - pos_fwd for pos_fwd, pos_rvs, length in matches]


class MyGenBankScanner(GenBankScanner):
    def _feed_first_line(self, consumer, line):
        # A regex for LOCUS       pKM265       4536 bp    DNA   circular  SYN        21-JUN-2013
        m = re.match(
            r'(?i)LOCUS\s+(?P<name>\S+)\s+(?P<size>\d+ bp)\s+(?P<molecule_type>\S+)(?:\s+(?P<topology>circular|linear))?(?:\s+.+\s+)?(?P<date>\d+-\w+-\d+)?',
            line,
        )
        if m is None:
            raise ValueError('LOCUS line cannot be parsed')
        name, size, molecule_type, topology, date = m.groups()

        consumer.locus(name)
        consumer.size(size[:-3])
        consumer.molecule_type(molecule_type)
        consumer.topology(topology.lower() if topology is not None else None)
        consumer.date(date)


class MyGenBankIterator(GenBankIterator):

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = MyGenBankScanner(debug=0).parse_records(handle)
        return records


def custom_file_parser(
    file_streamer: io.BytesIO | io.StringIO, sequence_file_format: SequenceFileFormat, circularize: bool = False
) -> list[Dseqrecord]:
    """
    Parse a file with SeqIO.parse (specifying the format and using the topology annotation to determine circularity).

    If the format is genbank and the parsing of the LOCUS line fails, fallback to custom regex-based parsing.
    """

    out = list()

    with file_streamer as handle:
        try:
            for parsed_seq in seqio_parse(handle, sequence_file_format):
                circularize = circularize or (
                    'topology' in parsed_seq.annotations.keys() and parsed_seq.annotations['topology'] == 'circular'
                )
                if sequence_file_format == 'genbank' and 'topology' not in parsed_seq.annotations.keys():
                    # If we could not parse the topology from the LOCUS line, raise an error to
                    # fallback to regex-based parsing
                    raise ValueError('LOCUS line does not contain topology')
                out.append(Dseqrecord(parsed_seq, circular=circularize))

        except ValueError as e:
            # If not locus-related error, raise
            if 'LOCUS line does not contain' not in str(e):
                raise e

            # If the error is about the LOCUS line, we try to parse with regex
            warnings.warn(
                'LOCUS line is wrongly formatted, we used a more permissive parser.',
                stacklevel=2,
            )
            # Reset the file handle position to the start since we consumed it in the first attempt
            handle.seek(0)
            out = list()
            for parsed_seq in MyGenBankIterator(handle):
                circularize = circularize or (
                    'topology' in parsed_seq.annotations.keys() and parsed_seq.annotations['topology'] == 'circular'
                )
                out.append(Dseqrecord(parsed_seq, circular=circularize))

    return out


async def get_sequence_from_euroscarf_url(plasmid_id: str) -> Dseqrecord:
    url = f'http://www.euroscarf.de/plasmid_details.php?accno={plasmid_id}'
    async with httpx.AsyncClient() as client:
        try:
            resp = await client.get(url)
        except httpx.ConnectError as e:
            raise HTTPError(url, 504, 'could not connect to euroscarf', 'could not connect to euroscarf', None) from e
    # I don't think this ever happens
    if resp.status_code != 200:
        raise HTTPError(
            url, resp.status_code, 'could not connect to euroscarf', 'could not connect to euroscarf', None
        )
    # Use beautifulsoup to parse the html
    soup = BeautifulSoup(resp.text, 'html.parser')
    # Identify if it's an error (seems to be a php error log without a body tag)
    body_tag = soup.find('body')
    if body_tag is None:
        if 'Call to a member function getName()' in resp.text:
            raise HTTPError(url, 404, 'invalid euroscarf id', 'invalid euroscarf id', None)
        else:
            msg = f'Could not retrieve plasmid details, double-check the euroscarf site: {url}'
            raise HTTPError(url, 503, msg, msg, None)
    # Get the download link
    subpath = soup.find('a', href=lambda x: x and x.startswith('files/dna'))
    if subpath is None:
        msg = f'Could not retrieve plasmid details, double-check the euroscarf site: {url}'
        raise HTTPError(url, 503, msg, msg, None)
    genbank_url = f'http://www.euroscarf.de/{subpath.get("href")}'
    return (await get_sequences_from_file_url(genbank_url))[0]


async def annotate_with_plannotate(
    file_content: str, file_name: str, url: str, timeout: int = 20
) -> tuple[Dseqrecord, PlannotateAnnotationReport, str]:
    async with httpx.AsyncClient() as client:
        try:
            response = await client.post(
                url,
                files={'file': (file_name, file_content, 'text/plain')},
                timeout=timeout,
            )
            if response.status_code != 200:
                detail = response.json().get('detail', 'plannotate server error')
                raise HTTPError(url, response.status_code, detail, detail, None)
            data = response.json()
            dseqr = custom_file_parser(io.StringIO(data['gb_file']), 'genbank')[0]
            report = [PlannotateAnnotationReport.model_validate(r) for r in data['report']]
            return dseqr, report, data['version']
        except httpx.TimeoutException as e:
            raise HTTPError(url, 504, 'plannotate server timeout', 'plannotate server timeout', None) from e
        except httpx.ConnectError as e:
            raise HTTPError(
                url, 500, 'cannot connect to plannotate server', 'cannot connect to plannotate server', None
            ) from e
