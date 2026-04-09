from fastapi import APIRouter
from .api_config_utils import RecordStubRoute
from .app_settings import settings


def get_router():
    if settings.RECORD_STUBS:
        return APIRouter(route_class=RecordStubRoute)
    else:
        return APIRouter()
