import unittest
import json
from fastapi import FastAPI, Request
from fastapi.testclient import TestClient

import opencloning.main as _main
from opencloning.api_config_utils import custom_http_exception_handler


client = TestClient(_main.app)


class GreetingTest(unittest.TestCase):
    def test_greeting(self):
        response = client.get('/')
        self.assertEqual(response.status_code, 200)
        self.assertIn('<a href="./docs">here</a>', response.text)


class InternalServerErrorTest(unittest.IsolatedAsyncioTestCase):

    async def test_internal_server_error(self):
        request = Request(scope={'type': 'http', 'headers': [(b'origin', b'http://localhost:3000')]})
        print(_main.settings.ALLOWED_ORIGINS)
        response = await _main.app.exception_handlers[500](request, None)
        self.assertEqual(response.headers['access-control-allow-credentials'], 'true')
        self.assertEqual(response.headers['access-control-allow-origin'], 'http://localhost:3000')


class CustomHttpExceptionHandlerTest(unittest.IsolatedAsyncioTestCase):
    async def test_internal_server_error_with_cors_headers(self):
        dummy_app = FastAPI()

        # Simulate a request from a specific origin not in the allowed list
        request = Request(scope={'type': 'http', 'headers': [(b'origin', b'http://example.com')]})
        response = await custom_http_exception_handler(request, None, dummy_app, ['http://localhost:3000', 'http://localhost:5173'])
        self.assertEqual(response.status_code, 500)
        self.assertEqual(json.loads(response.body), {'message': 'internal server error'})
        self.assertEqual(response.headers['access-control-allow-credentials'], 'true')
        self.assertNotIn('access-control-allow-origin', response.headers)

        # All origins allowed: per CORS spec the specific origin is echoed back
        # (Access-Control-Allow-Origin cannot be '*' when credentials are allowed)
        request = Request(scope={'type': 'http', 'headers': [(b'origin', b'http://example.com')]})
        response = await custom_http_exception_handler(request, None, dummy_app, ['*'])
        self.assertEqual(response.status_code, 500)
        self.assertEqual(json.loads(response.body), {'message': 'internal server error'})
        self.assertEqual(response.headers['access-control-allow-credentials'], 'true')
        self.assertEqual(response.headers['access-control-allow-origin'], 'http://example.com')

        # Origin in the allowed list: specific origin is returned
        request = Request(scope={'type': 'http', 'headers': [(b'origin', b'http://localhost:3000')]})
        response = await custom_http_exception_handler(request, None, dummy_app, ['http://localhost:3000', 'http://localhost:5173'])
        self.assertEqual(response.status_code, 500)
        self.assertEqual(json.loads(response.body), {'message': 'internal server error'})
        self.assertEqual(response.headers['access-control-allow-credentials'], 'true')
        self.assertEqual(response.headers['access-control-allow-origin'], 'http://localhost:3000')

        # All origins allowed with a cookie header: specific origin is still returned
        request = Request(scope={'type': 'http', 'headers': [(b'origin', b'http://localhost:3000'), (b'cookie', b'session=abc123')]})
        response = await custom_http_exception_handler(request, None, dummy_app, ['*'])
        self.assertEqual(response.status_code, 500)
        self.assertEqual(json.loads(response.body), {'message': 'internal server error'})
        self.assertEqual(response.headers['access-control-allow-credentials'], 'true')
        self.assertEqual(response.headers['access-control-allow-origin'], 'http://localhost:3000')
