"""
This module contains the settings for the app that can be set via environment variables.
"""

import os
from pydantic import BaseModel


def parse_bool(value: str) -> bool:
    return value in {'1', 'TRUE', 'true', 'True'}


# API settings ===============================================
SERVE_FRONTEND = parse_bool(os.environ['SERVE_FRONTEND']) if 'SERVE_FRONTEND' in os.environ else False
BATCH_CLONING = parse_bool(os.environ['BATCH_CLONING']) if 'BATCH_CLONING' in os.environ else True
RECORD_STUBS = parse_bool(os.environ['RECORD_STUBS']) if 'RECORD_STUBS' in os.environ else False
ALLOWED_ORIGINS = ['http://localhost:3000', 'http://localhost:5173']
if os.environ.get('ALLOWED_ORIGINS') is not None:
    # Remove trailing slash from each origin if ends with one
    ALLOWED_ORIGINS = [origin.rstrip('/') for origin in os.environ['ALLOWED_ORIGINS'].split(',')]


# External services settings =================================
NCBI_API_KEY = os.environ.get('NCBI_API_KEY')
PLANNOTATE_URL = os.environ['PLANNOTATE_URL'] if 'PLANNOTATE_URL' in os.environ else None
PLANNOTATE_TIMEOUT = int(os.environ['PLANNOTATE_TIMEOUT']) if 'PLANNOTATE_TIMEOUT' in os.environ else 20
# Handle trailing slash:
if PLANNOTATE_URL is not None and not PLANNOTATE_URL.endswith('/'):
    PLANNOTATE_URL += '/'

PROXY_URL = os.environ.get('PROXY_URL')
PROXY_CERT_FILE = os.environ.get('PROXY_CERT_FILE')


class Settings(BaseModel):
    SERVE_FRONTEND: bool
    BATCH_CLONING: bool
    RECORD_STUBS: bool
    NCBI_API_KEY: str | None
    ALLOWED_ORIGINS: list[str]
    PLANNOTATE_URL: str | None
    PLANNOTATE_TIMEOUT: int
    PROXY_URL: str | None
    # Must be a full path to the proxy certificate file
    PROXY_CERT_FILE: str | None


settings = Settings(
    SERVE_FRONTEND=SERVE_FRONTEND,
    BATCH_CLONING=BATCH_CLONING,
    RECORD_STUBS=RECORD_STUBS,
    NCBI_API_KEY=NCBI_API_KEY,
    ALLOWED_ORIGINS=ALLOWED_ORIGINS,
    PLANNOTATE_URL=PLANNOTATE_URL,
    PLANNOTATE_TIMEOUT=PLANNOTATE_TIMEOUT,
    PROXY_URL=PROXY_URL,
    PROXY_CERT_FILE=PROXY_CERT_FILE,
)
