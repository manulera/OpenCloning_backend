from fastapi import Form, File, UploadFile, HTTPException
from typing import Annotated, Literal
from tempfile import TemporaryDirectory
import os
from fastapi.responses import FileResponse
from .pombe_get_primers import main as pombe_primers
from .pombe_clone import main as pombe_clone
from .pombe_summary import main as pombe_summary
from .pombe_gather import main as pombe_gather
import shutil
import traceback
from ...get_router import get_router
from fastapi import Request

router = get_router()


@router.get('/batch_cloning/pombe')
async def get_batch_cloning_page(request: Request):
    return FileResponse(os.path.join(os.path.dirname(__file__), 'index.html'))


@router.post('/batch_cloning/pombe')
async def post_batch_cloning(
    gene_list: str = Form(...),
    plasmid_file: UploadFile | None = File(None),
    addgene_id: str | None = Form(None),
    plasmid_option: Annotated[Literal['addgene', 'file'], Form(...)] = None,
    checking_primer_forward: str = Form(..., pattern=r'^[ACGTacgt]+$', min_length=1),
    checking_primer_reverse: str = Form(..., pattern=r'^[ACGTacgt]+$', min_length=1),
):

    plasmid = plasmid_file if plasmid_option == 'file' else addgene_id
    if plasmid is None:
        raise HTTPException(status_code=400, detail='No plasmid provided')

    genes = [gene.strip() for gene in gene_list.split() if gene.strip()]

    if not genes:
        raise HTTPException(status_code=400, detail='No valid genes provided')

    with TemporaryDirectory() as temp_dir:
        if plasmid_option == 'file':
            # Write the plasmid to the temp dir
            with open(os.path.join(temp_dir, plasmid_file.filename), 'wb') as f:
                shutil.copyfileobj(plasmid_file.file, f)

        # Write the checking primers
        with open(os.path.join(temp_dir, 'checking_primers.fa'), 'w') as f:
            f.write(f'>common_insert_fwd\n{checking_primer_forward}\n>common_insert_rvs\n{checking_primer_reverse}')

        for gene in genes:
            try:
                await pombe_primers(gene, temp_dir)
            except Exception:
                raise HTTPException(status_code=404, detail=f'Primers for {gene} not found')
            try:
                if plasmid_option == 'file':
                    with open(os.path.join(temp_dir, plasmid_file.filename), 'rb') as f:
                        await pombe_clone(
                            gene, 'GCF_000002945.2', temp_dir, UploadFile(file=f, filename=plasmid_file.filename)
                        )
                else:
                    await pombe_clone(gene, 'GCF_000002945.2', temp_dir, addgene_id)
            except Exception:
                # Show the stack trace in console
                print(f'Error occurred while cloning {gene}:')
                traceback.print_exc()
                raise HTTPException(status_code=400, detail=f'Clone for {gene} failed')
        try:
            pombe_summary(temp_dir)
            pombe_gather(temp_dir)
        except Exception:
            raise HTTPException(status_code=400, detail='Summary failed')

        # zip the temp dir and return it
        zip_filename = f'{temp_dir}_archive'
        shutil.make_archive(zip_filename, 'zip', temp_dir)
        zip_file = f'{zip_filename}.zip'
        return FileResponse(zip_file, filename='batch_cloning_output.zip')
