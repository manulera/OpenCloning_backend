from fastapi import Body, Query, HTTPException, Response, UploadFile, File
from pydantic import create_model
import io
import warnings
import asyncio
import httpx
from starlette.responses import RedirectResponse
from Bio import BiopythonParserWarning
from typing import Annotated
from urllib.error import HTTPError
from pydna.utils import location_boundaries

from ..get_router import get_router
from ..pydantic_models import (
    TextFileSequence,
    UploadedFileSource,
    RepositoryIdSource,
    AddgeneIdSource,
    WekWikGeneIdSource,
    BenchlingUrlSource,
    SnapGenePlasmidSource,
    EuroscarfSource,
    IGEMSource,
    GenomeCoordinatesSource,
    SequenceFileFormat,
    SEVASource,
    SequenceLocationStr,
)
from ..dna_functions import (
    format_sequence_genbank,
    request_from_addgene,
    request_from_wekwikgene,
    get_sequences_from_file_url,
    get_sequence_from_snagene_url,
    custom_file_parser,
    get_sequence_from_euroscarf_url,
    get_seva_plasmid,
)
from .. import request_examples
from .. import ncbi_requests


router = get_router()


# TODO limit the maximum size of submitted files
@router.post(
    '/read_from_file',
    response_model=create_model(
        'UploadedFileResponse', sources=(list[UploadedFileSource], ...), sequences=(list[TextFileSequence], ...)
    ),
    responses={
        200: {
            'description': 'The sequence was successfully parsed',
            'headers': {
                'x-warning': {
                    'description': 'A warning returned if the file can be read but is not in the expected format or if some sequences were not extracted because they are incompatible with the provided coordinates',
                    'schema': {'type': 'string'},
                },
            },
        },
        422: {
            'description': 'Biopython cannot process this file or provided coordinates are invalid.',
        },
        404: {
            'description': 'The index_in_file is out of range.',
        },
    },
)
async def read_from_file(
    response: Response,
    file: UploadFile = File(...),
    sequence_file_format: SequenceFileFormat | None = Query(
        None,
        description='Format of the sequence file. Unless specified, it will be guessed from the extension',
    ),
    index_in_file: int | None = Query(
        None,
        description='The index of the sequence in the file for multi-sequence files',
    ),
    circularize: bool = Query(
        False,
        description='circularize the sequence read (for GenBank or Snapgene files, it will override the topology indicated in the file)',
    ),
    output_name: str | None = Query(
        None,
        description='Name of the output sequence',
    ),
    start: int | None = Query(None, description='Start position of the sequence to read (0-based)', ge=0),
    end: int | None = Query(
        None,
        description='End position of the sequence to read (0-based)',
        ge=0,
    ),
):
    """Return a json sequence from a sequence file"""

    if sequence_file_format is None:
        extension_dict = {
            'gbk': 'genbank',
            'gb': 'genbank',
            'ape': 'genbank',
            'dna': 'snapgene',
            'fasta': 'fasta',
            'embl': 'embl',
            'fa': 'fasta',
        }
        extension = file.filename.split('.')[-1].lower()
        if extension not in extension_dict:
            raise HTTPException(
                422,
                'We could not guess the format of the file from its extension. Please provide file_format as a query parameter.',
            )

        # We guess the file type from the extension
        sequence_file_format = SequenceFileFormat(extension_dict[extension])

    dseqs = list()
    warning_messages = list()

    file_content = await file.read()
    if sequence_file_format == 'snapgene':
        file_streamer = io.BytesIO(file_content)
    else:
        file_streamer = io.StringIO(file_content.decode())

    try:
        # Capture warnings without converting to errors:
        with warnings.catch_warnings(record=True, category=UserWarning) as warnings_captured:
            dseqs = custom_file_parser(file_streamer, sequence_file_format, circularize)

        # If there were warnings, add them to the response header
        warnings_captured = [w for w in warnings_captured if w.category is not BiopythonParserWarning]

        if warnings_captured:
            warning_messages = [str(w.message) for w in warnings_captured]

    except ValueError as e:
        raise HTTPException(422, f'Biopython cannot process this file: {e}.')

    # This happens when textfiles are empty or contain something else, or when reading a text file as snapgene file,
    # since StringIO does not raise an error when "Unexpected end of packet" is found
    if len(dseqs) == 0:
        raise HTTPException(422, 'Biopython cannot process this file.')

    if index_in_file is not None:
        if index_in_file >= len(dseqs):
            raise HTTPException(404, 'The index_in_file is out of range.')
        dseqs = [dseqs[index_in_file]]

    seq_feature = None
    if start is not None and end is not None:
        extracted_sequences = list()
        for dseq in dseqs:
            try:
                seq_feature = SequenceLocationStr.from_start_and_end(start=start, end=end, seq_len=len(dseq))
                # TODO: We could use extract when this is addressed: https://github.com/biopython/biopython/issues/4989
                location = seq_feature.to_biopython_location()
                i, j = location_boundaries(location)
                extracted_sequence = dseq[i:j]
                # Only add the sequence if the interval is not out of bounds
                if len(extracted_sequence) == len(location):
                    extracted_sequences.append(extracted_sequence)
                else:
                    extracted_sequences.append(None)
            except Exception:
                extracted_sequences.append(None)
        dseqs = extracted_sequences

    # The common part
    parent_source = UploadedFileSource(
        id=0,
        sequence_file_format=sequence_file_format,
        file_name=file.filename,
        circularize=circularize,
        coordinates=seq_feature,
    )

    # If coordinates are provided, we only keep the sequences compatible with those coordinates
    out_sources = list()
    out_sequences = list()
    for i in range(len(dseqs)):
        if dseqs[i] is None:
            continue
        new_source = parent_source.model_copy()
        new_source.index_in_file = index_in_file if index_in_file is not None else i
        out_sources.append(new_source)
        out_sequences.append(format_sequence_genbank(dseqs[i], output_name))

    if len(out_sequences) == 0:
        raise HTTPException(422, 'Provided coordinates are incompatible with sequences in the file.')

    if len(out_sequences) < len(dseqs):
        warning_messages.append(
            'Some sequences were not extracted because they are incompatible with the provided coordinates.'
        )

    if len(warning_messages) > 0:
        response.headers['x-warning'] = '; '.join(warning_messages)

    return {'sequences': out_sequences, 'sources': out_sources}


# TODO: a bit inconsistent that here you don't put {source: {...}} in the request, but
# directly the object.


def repository_id_http_error_handler(exception: HTTPError, source: RepositoryIdSource):

    if exception.code == 500:  # pragma: no cover
        raise HTTPException(
            503, f'{source.repository_name} returned: {exception} - {source.repository_name} might be down'
        )
    elif exception.code == 400 or exception.code == 404:
        raise HTTPException(
            404,
            f'{source.repository_name} returned: {exception} - Likely you inserted a wrong {source.repository_name} id',
        )


# Redirect to the right repository
@router.post(
    '/repository_id',
    response_model=create_model(
        'RepositoryIdResponse',
        sources=(
            list[RepositoryIdSource]
            | list[AddgeneIdSource]
            | list[BenchlingUrlSource]
            | list[EuroscarfSource]
            | list[WekWikGeneIdSource]
            | list[SEVASource],
            ...,
        ),
        sequences=(list[TextFileSequence], ...),
    ),
)
async def get_from_repository_id(
    source: (
        RepositoryIdSource
        | AddgeneIdSource
        | BenchlingUrlSource
        | SnapGenePlasmidSource
        | EuroscarfSource
        | WekWikGeneIdSource
        | SEVASource
    ),
):
    return RedirectResponse(f'/repository_id/{source.repository_name}', status_code=307)


@router.post(
    '/repository_id/genbank',
    response_model=create_model(
        'RepositoryIdResponse', sources=(list[RepositoryIdSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def get_from_repository_id_genbank(source: RepositoryIdSource):
    try:
        # This request already fails if the sequence does not exist
        seq_length = await ncbi_requests.get_sequence_length_from_sequence_accession(source.repository_id)
        if seq_length > 100000:
            raise HTTPException(400, 'sequence is too long (max 100000 bp)')
        seq = await ncbi_requests.get_genbank_sequence(source.repository_id)
    except httpx.ConnectError as exception:
        raise HTTPException(504, f'Unable to connect to NCBI: {exception}')

    return {'sequences': [format_sequence_genbank(seq, source.output_name)], 'sources': [source.model_copy()]}


@router.post(
    '/repository_id/addgene',
    response_model=create_model(
        'AddgeneIdResponse', sources=(list[AddgeneIdSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def get_from_repository_id_addgene(source: AddgeneIdSource):
    try:
        dseq, out_source = await request_from_addgene(source)
    except HTTPError as exception:
        repository_id_http_error_handler(exception, source)
    except httpx.ConnectError:
        raise HTTPException(504, 'unable to connect to Addgene')

    return {'sequences': [format_sequence_genbank(dseq, source.output_name)], 'sources': [out_source]}


@router.post(
    '/repository_id/wekwikgene',
    response_model=create_model(
        'WekWikGeneIdResponse', sources=(list[WekWikGeneIdSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def get_from_repository_id_wekwikgene(source: WekWikGeneIdSource):
    try:
        dseq, out_source = await request_from_wekwikgene(source)
    except HTTPError as exception:
        repository_id_http_error_handler(exception, source)
    except httpx.ConnectError:
        raise HTTPException(504, 'unable to connect to WekWikGene')
    return {'sequences': [format_sequence_genbank(dseq, source.output_name)], 'sources': [out_source]}


@router.post(
    '/repository_id/benchling',
    response_model=create_model(
        'BenchlingUrlResponse', sources=(list[BenchlingUrlSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def get_from_benchling_url(
    source: Annotated[BenchlingUrlSource, Body(openapi_examples=request_examples.benchling_url_examples)]
):
    try:
        dseqs = await get_sequences_from_file_url(source.repository_id)
        return {
            'sequences': [format_sequence_genbank(s, source.output_name) for s in dseqs],
            'sources': [source for s in dseqs],
        }
    except HTTPError as exception:
        repository_id_http_error_handler(exception, source)


@router.post(
    '/repository_id/snapgene',
    response_model=create_model(
        'SnapGenePlasmidResponse', sources=(list[SnapGenePlasmidSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def get_from_repository_id_snapgene(
    source: Annotated[SnapGenePlasmidSource, Body(openapi_examples=request_examples.snapgene_plasmid_examples)]
):
    try:
        plasmid_set, plasmid_name = source.repository_id.split('/')
        url = f'https://www.snapgene.com/local/fetch.php?set={plasmid_set}&plasmid={plasmid_name}'
        dseq = get_sequence_from_snagene_url(url)
        # Unless a name is provided, we use the plasmid name from snapgene
        if source.output_name is None:
            source.output_name = plasmid_name
        return {
            'sequences': [format_sequence_genbank(dseq, source.output_name)],
            'sources': [source],
        }
    except HTTPError as exception:
        repository_id_http_error_handler(exception, source)


@router.post(
    '/repository_id/euroscarf',
    response_model=create_model(
        'EuroscarfResponse', sources=(list[EuroscarfSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def get_from_repository_id_euroscarf(source: EuroscarfSource):
    """
    Return the sequence from a plasmid in Euroscarf. Sometimes plasmid files do not contain correct topology information
    (they indicate linear sequence instead of circular). We force them to be circular.
    """
    try:
        dseq = await get_sequence_from_euroscarf_url(source.repository_id)
        # Sometimes the files do not contain correct topology information, so we loop them
        if not dseq.circular:
            dseq = dseq.looped()
        return {'sequences': [format_sequence_genbank(dseq, source.output_name)], 'sources': [source]}
    except HTTPError as exception:
        repository_id_http_error_handler(exception, source)


@router.post(
    '/repository_id/igem',
    response_model=create_model(
        'IGEMResponse', sources=(list[IGEMSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def get_from_repository_id_igem(source: IGEMSource):
    try:
        dseq = (await get_sequences_from_file_url(source.sequence_file_url))[0]
        return {'sequences': [format_sequence_genbank(dseq, source.output_name)], 'sources': [source]}
    except HTTPError as exception:
        repository_id_http_error_handler(exception, source)


@router.post(
    '/genome_coordinates',
    response_model=create_model(
        'GenomeRegionResponse', sources=(list[GenomeCoordinatesSource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def genome_coordinates(
    source: Annotated[GenomeCoordinatesSource, Body(openapi_examples=request_examples.genome_region_examples)]
):

    # Validate that coordinates make sense
    ncbi_requests.validate_coordinates_pre_request(source.start, source.end, source.strand)

    # Source includes a locus tag in annotated assembly

    async def validate_locus_task():
        if source.locus_tag is not None:

            if source.assembly_accession is None:
                raise HTTPException(422, 'assembly_accession is required if locus_tag is set')

            annotation = await ncbi_requests.get_annotation_from_locus_tag(source.locus_tag, source.assembly_accession)
            gene_range = annotation['genomic_regions'][0]['gene_range']['range'][0]
            gene_strand = 1 if gene_range['orientation'] == 'plus' else -1

            # This field will not be present in all cases, but should be there in reference genomes
            if source.gene_id is not None:
                if 'gene_id' not in annotation:
                    raise HTTPException(400, 'gene_id is set, but not found in the annotation')
                if source.gene_id != int(annotation['gene_id']):
                    raise HTTPException(400, 'gene_id does not match the locus_tag')
            elif 'gene_id' in annotation:
                source.gene_id = int(annotation['gene_id'])

            # The gene should fall within the range (range might be bigger if bases were requested upstream or downstream)
            if (
                int(gene_range['begin']) < source.start
                or int(gene_range['end']) > source.end
                or gene_strand != source.strand
            ):
                raise HTTPException(
                    400,
                    f'wrong coordinates, expected to fall within {source.start}, {source.end} on strand: {source.strand}',
                )

    async def validate_assembly_task():
        if source.assembly_accession is not None:
            # We get the assembly accession (if it exists), and if the user provided one we validate it
            sequence_accessions = await ncbi_requests.get_sequence_accessions_from_assembly_accession(
                source.assembly_accession
            )
            if source.sequence_accession not in sequence_accessions:
                raise HTTPException(
                    400,
                    f'Sequence accession {source.sequence_accession} not contained in assembly accession {source.assembly_accession}, which contains accessions: {", ".join(sequence_accessions)}',
                )

    async def get_sequence_task():
        return await ncbi_requests.get_genbank_sequence(
            source.sequence_accession, source.start, source.end, source.strand
        )

    tasks = [validate_locus_task(), validate_assembly_task(), get_sequence_task()]

    _, _, seq = await asyncio.gather(*tasks)

    # NCBI does not complain for coordinates that fall out of the sequence, so we have to check here
    if len(seq) != source.end - source.start + 1:
        raise HTTPException(400, 'coordinates fall outside the sequence')

    return {'sequences': [format_sequence_genbank(seq, source.output_name)], 'sources': [source.model_copy()]}


@router.post(
    '/repository_id/seva',
    response_model=create_model(
        'SEVASource', sources=(list[SEVASource], ...), sequences=(list[TextFileSequence], ...)
    ),
)
async def get_from_repository_id_seva(source: SEVASource):
    """
    Return the sequence from a plasmid in SEVA.
    """
    try:
        dseq, source = await get_seva_plasmid(source)
        return {'sequences': [format_sequence_genbank(dseq, source.output_name)], 'sources': [source]}
    except HTTPError as exception:
        repository_id_http_error_handler(exception, source)
    except httpx.ConnectError:
        raise HTTPException(504, 'unable to connect to SEVA')
    except Exception as exception:
        raise HTTPException(400, f'Error parsing file: {exception}')
