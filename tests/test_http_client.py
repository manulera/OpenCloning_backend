import unittest
from importlib import reload
from opencloning import httpClient
from opencloning import app_settings
from pytest import MonkeyPatch
from unittest.mock import patch, AsyncMock
import os

mock_cert_file = os.path.join(os.path.dirname(__file__), 'test_files/dummy_certificate.crt')


class TestHttpClientProxy(unittest.IsolatedAsyncioTestCase):

    def tearDown(self):
        MonkeyPatch().delenv('PROXY_URL', raising=False)
        MonkeyPatch().delenv('PROXY_CERT_FILE', raising=False)

    async def test_proxy(self):
        MonkeyPatch().setenv('PROXY_URL', 'proxy_url')
        MonkeyPatch().setenv('PROXY_CERT_FILE', mock_cert_file)
        reload(app_settings)
        reload(httpClient)
        mock_client = AsyncMock()
        with patch('httpx.AsyncClient', return_value=mock_client) as mock_client_class:
            async with httpClient.get_http_client():
                pass

            # Verify the client was initialized with correct proxy settings
            mock_client_class.assert_called_once()
            called_kwargs = mock_client_class.call_args.kwargs
            self.assertIn('proxy', called_kwargs)
            self.assertIn('verify', called_kwargs)
            self.assertEqual(called_kwargs['proxy'], 'proxy_url')

        # To ensure that the certificate is actually being loaded,
        # provide invalid path and test the error
        MonkeyPatch().setenv('PROXY_CERT_FILE', 'invalid_path')
        reload(app_settings)
        with self.assertRaises(FileNotFoundError):
            reload(httpClient)
