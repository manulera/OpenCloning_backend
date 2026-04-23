from pydantic import BaseModel


class StubRequest(BaseModel):
    endpoint: str
    method: str
    name: str
    params: dict | None = None
    body: dict | None = None
    headers: dict | None = None
    reset_db: bool = False


stubs = [
    StubRequest(
        name='get_primers',
        endpoint='/primers',
        method='GET',
    ),
    StubRequest(
        name='get_primer',
        endpoint='/primer/7',
        method='GET',
    ),
    StubRequest(
        name='post_primer',
        endpoint='/primer',
        method='POST',
        body={'id': 0, 'name': 'new', 'sequence': 'GGCC'},
        reset_db=True,
    ),
    StubRequest(
        name='patch_primer',
        endpoint='/primer/7',
        method='PATCH',
        body={'name': 'fwd_renamed'},
        reset_db=True,
    ),
    StubRequest(
        name='get_sequences',
        endpoint='/sequences',
        method='GET',
    ),
    StubRequest(
        name='get_sequence',
        endpoint='/sequence/10',
        method='GET',
    ),
    StubRequest(
        name='patch_sequence',
        endpoint='/sequence/10',
        method='PATCH',
        body={'name': 'pcr_product_renamed'},
        reset_db=True,
    ),
    StubRequest(
        name='get_sequence_by_uid',
        endpoint='/sequence/by-uid/example_sequencing-sample',
        method='GET',
    ),
    StubRequest(
        name='get_sequences_by_seguid',
        endpoint='/sequences/by-seguid/ldseguid=oMGruVpBiElY0ffP28XC_BlHXv8',
        method='GET',
    ),
    StubRequest(
        name='get_text_file_sequence',
        endpoint='/sequence/10/text_file_sequence',
        method='GET',
    ),
    StubRequest(
        name='get_cloning_strategy',
        endpoint='/sequence/10/cloning_strategy',
        method='GET',
    ),
    StubRequest(
        name='get_sequence_primers',
        endpoint='/sequence/10/primers',
        method='GET',
    ),
]
