import glob
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

from .get_router import get_router
from .api_config_utils import custom_http_exception_handler as _custom_http_exception_handler

from .endpoints.primer_design import router as primer_design_router
from .endpoints.external_import import router as import_router
from .endpoints.other import router as other_router
from .endpoints.annotation import router as annotation_router
from .endpoints.assembly import router as assembly_router
from .endpoints.no_assembly import router as no_assembly_router
from .endpoints.no_input import router as no_input_router
from .app_settings import settings

# =====================================================

# Instance of the API object
app = FastAPI()

router = get_router()
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=['*'],
    allow_headers=['*'],
    expose_headers=['x-warning'],
)


@app.exception_handler(500)
async def custom_http_exception_handler(request: Request, exc: Exception):
    return await _custom_http_exception_handler(request, exc, app, settings.ALLOWED_ORIGINS)


if not settings.SERVE_FRONTEND:

    @router.get('/')
    async def greeting(request: Request):
        html_content = """
            <html>
                <head>
                    <title>Welcome to OpenCloning API</title>
                </head>
                <body>
                    <h1>Welcome to OpenCloning API</h1>
                    <p>You can access the endpoints documentation <a href="./docs">here</a></p>
                </body>
            </html>
            """
        return HTMLResponse(content=html_content, status_code=200)

else:
    app.mount('/assets', StaticFiles(directory='frontend/assets'), name='assets')
    app.mount('/examples', StaticFiles(directory='frontend/examples'), name='examples')

    @router.get('/')
    async def get_frontend_index(request: Request):
        return FileResponse('frontend/index.html')

    frontend_files = (
        glob.glob('frontend/*.json')
        + glob.glob('frontend/*.ico')
        + glob.glob('frontend/*.png')
        + glob.glob('frontend/*.txt')
    )
    frontend_files = [f.split('/')[-1] for f in frontend_files]

    @router.get('/{name:path}')
    async def get_other_frontend_files(name: str):
        """Catch-all for frontend files"""
        if name in frontend_files:
            return FileResponse(f'frontend/{name}')
        raise HTTPException(404)


app.include_router(import_router, tags=['External Import'])
app.include_router(assembly_router, tags=['Assembly'])
app.include_router(no_assembly_router, tags=['No Assembly'])
app.include_router(other_router, tags=['Other'])
app.include_router(no_input_router, tags=['No Input'])
app.include_router(primer_design_router, tags=['Primer Design'])
app.include_router(annotation_router, tags=['Annotation'])

if settings.BATCH_CLONING:
    from .batch_cloning import router as batch_cloning_router
    from .batch_cloning.ziqiang_et_al2024 import router as ziqiang_et_al2024_router
    from .batch_cloning.pombe import router as pombe_router

    app.include_router(batch_cloning_router, tags=['Batch Cloning'])
    app.include_router(ziqiang_et_al2024_router, tags=['Batch Cloning'])
    app.include_router(pombe_router, tags=['Batch Cloning'])


# This router must be added last because when SERVE_FRONTEND is True,
# it contains a catch-all route. The catch-all route '/{name:path}' matches any URL path, so if this
# section were placed earlier, it would intercept all requests before they could reach their intended
# API endpoints. For example, requests to '/docs' or '/version' would incorrectly return 404 errors
# instead of reaching their proper handlers.
app.include_router(router, tags=['General'])
