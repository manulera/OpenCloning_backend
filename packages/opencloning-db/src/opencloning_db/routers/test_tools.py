"""Test-only backend utilities."""

from __future__ import annotations

import os
import subprocess

from fastapi import APIRouter, Header, HTTPException, Response, status

router = APIRouter(prefix='/__test', tags=['test-tools'])


@router.post('/reset-db', status_code=status.HTTP_204_NO_CONTENT)
def reset_db(x_test_reset_token: str | None = Header(default=None)) -> Response:
    """Reset DB by delegating to ``opencloning-cli db test reset``."""

    if not os.getenv('OPENCLONING_TESTING'):
        raise HTTPException(status_code=403, detail='Forbidden')

    result = subprocess.run(
        ['opencloning-cli', 'db', 'test', 'reset'],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or 'Reset command failed'
        raise HTTPException(status_code=500, detail=message)

    return Response(status_code=status.HTTP_204_NO_CONTENT)
