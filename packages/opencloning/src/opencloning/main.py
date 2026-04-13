import glob
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

from .get_router import get_router

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
_app = FastAPI()
app = CORSMiddleware(
    _app,
    allow_origins=settings.ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=['*'],
    allow_headers=['*'],
    expose_headers=['x-warning'],
)

router = get_router()


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
    from .app_settings import frontend_config, FrontendConfig

    _app.mount('/assets', StaticFiles(directory='frontend/assets'), name='assets')
    _app.mount('/examples', StaticFiles(directory='frontend/examples'), name='examples')

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

    @router.get('/config.json', response_model=FrontendConfig)
    async def get_config_json():
        """Return frontend config file built from env vars"""
        return frontend_config

    @router.get('/{name:path}')
    async def get_other_frontend_files(name: str):
        """Catch-all for frontend files"""
        if name in frontend_files:
            return FileResponse(f'frontend/{name}')
        raise HTTPException(404)


_app.include_router(import_router, tags=['External Import'])
_app.include_router(assembly_router, tags=['Assembly'])
_app.include_router(no_assembly_router, tags=['No Assembly'])
_app.include_router(other_router, tags=['Other'])
_app.include_router(no_input_router, tags=['No Input'])
_app.include_router(primer_design_router, tags=['Primer Design'])
_app.include_router(annotation_router, tags=['Annotation'])

if settings.BATCH_CLONING:
    from .batch_cloning import router as batch_cloning_router
    from .batch_cloning.ziqiang_et_al2024 import router as ziqiang_et_al2024_router
    from .batch_cloning.pombe import router as pombe_router

    _app.include_router(batch_cloning_router, tags=['Batch Cloning'])
    _app.include_router(ziqiang_et_al2024_router, tags=['Batch Cloning'])
    _app.include_router(pombe_router, tags=['Batch Cloning'])


# This router must be added last because when SERVE_FRONTEND is True,
# it contains a catch-all route. The catch-all route '/{name:path}' matches any URL path, so if this
# section were placed earlier, it would intercept all requests before they could reach their intended
# API endpoints. For example, requests to '/docs' or '/version' would incorrectly return 404 errors
# instead of reaching their proper handlers.
_app.include_router(router, tags=['General'])
