# Response is imported into other modules from here
from httpx import Response, AsyncClient, ConnectError, TimeoutException  # noqa: F401
import ssl
import certifi
from .app_settings import settings

ctx = None
proxy = None
if settings.PROXY_URL:
    ctx = ssl.create_default_context(cafile=certifi.where())
    # See https://www.python-httpx.org/advanced/ssl/
    if settings.PROXY_CERT_FILE:
        ctx.load_verify_locations(cafile=settings.PROXY_CERT_FILE)
    proxy = settings.PROXY_URL


def get_http_client():
    return AsyncClient(proxy=proxy, verify=ctx)
