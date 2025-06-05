# Response is imported into other modules from here
from httpx import (  # noqa: F401
    Response,
    AsyncClient,
    ConnectError,
    TimeoutException,
    AsyncHTTPTransport,
    Request,
    RequestError,
)
import ssl
import certifi
from .app_settings import settings
import re

white_listed_urls = {
    r'^https://www.addgene.org/',
    r'^https://media.addgene.org/',
    r'^https://wekwikgene.wllsb.edu.cn',
    r'^https://seva-plasmids.com/',
    r'^https://api.ncbi.nlm.nih.gov/datasets/v2alpha/',
    r'^https://eutils.ncbi.nlm.nih.gov/entrez/eutils/',
    r'^https://www.snapgene.com/local/fetch.php',
    r'^https://benchling.com/',
    r'^https://raw.githubusercontent.com/manulera/annotated-igem-distribution',
    r'^http://www.euroscarf.de/',
}

if settings.PLANNOTATE_URL:
    white_listed_urls.add(settings.PLANNOTATE_URL)


class WhiteListTransport(AsyncHTTPTransport):
    async def handle_async_request(self, request: Request) -> Response:
        if any(re.match(url, str(request.url)) for url in white_listed_urls):
            return await super().handle_async_request(request)

        raise RequestError(f'Request to {request.url} is not whitelisted')


transport = WhiteListTransport()

ctx = None
proxy = None
if settings.PROXY_URL:
    ctx = ssl.create_default_context(cafile=certifi.where())
    # See https://www.python-httpx.org/advanced/ssl/
    if settings.PROXY_CERT_FILE:
        ctx.load_verify_locations(cafile=settings.PROXY_CERT_FILE)
    proxy = settings.PROXY_URL


def get_http_client():
    return AsyncClient(proxy=proxy, verify=ctx, transport=transport)
