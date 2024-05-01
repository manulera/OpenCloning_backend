import json
import datetime
from fastapi.exceptions import RequestValidationError, HTTPException
from typing import Callable
from fastapi.routing import APIRoute
from fastapi import Request, Response
from fastapi.responses import JSONResponse
import os
import glob


class RecordStubRoute(APIRoute):
    """Subclass of APIRoute that stores the request and response of a route in the folder `stubs`"""

    def get_route_handler(self) -> Callable:
        original_route_handler = super().get_route_handler()

        async def custom_route_handler(request: Request) -> Response:
            if request.method != 'POST':
                return await original_route_handler(request)
            formatted_request = {
                'path': request.url.path,
                'method': request.method,
                'body': await request.json(),
                'headers': dict(request.headers),
            }
            try:
                response: JSONResponse = await original_route_handler(request)
            except RequestValidationError as exc:
                errors = exc.errors()
                for e in errors:
                    if 'ctx' in e:
                        # This is a ValueError Object that cannot be serialized
                        del e['ctx']
                response = JSONResponse(content=errors, status_code=422)
            except HTTPException as exc:
                print('>>exception', exc)
                detail = {'detail': exc.detail}
                response = JSONResponse(content=detail, status_code=exc.status_code)

            if type(response) is JSONResponse:
                formatted_response = {
                    'statusCode': response.status_code,
                    'body': json.loads(response.body),
                    'headers': dict(response.headers),
                }

                formatted_time = datetime.datetime.now().strftime('%Y_%m_%d-%H_%M_%S')
                stub_folder = f'stubs{formatted_request["path"]}/{formatted_time}'
                existing_folders = glob.glob(f'{stub_folder}_*')
                stub_folder = f'{stub_folder}_{len(existing_folders)}'
                if not os.path.exists(stub_folder):
                    os.makedirs(stub_folder)
                with open(f'{stub_folder}/request.json', 'w') as f:
                    json.dump(formatted_request, f, indent=4)
                with open(f'{stub_folder}/response.json', 'w') as f:
                    json.dump(formatted_response, f, indent=4)
                with open(f'{stub_folder}/response_body.json', 'w') as f:
                    json.dump(formatted_response['body'], f, indent=4)

            print(9 * ' ', '> stub written to', stub_folder)
            return response

        return custom_route_handler
