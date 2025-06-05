import unittest
import pytest
from importlib import reload
import os

from opencloning import app_settings


class TestAppSettings(unittest.TestCase):
    def setUp(self):
        # Store original environment variables
        self.original_env = {
            'SERVE_FRONTEND': os.getenv('SERVE_FRONTEND'),
            'BATCH_CLONING': os.getenv('BATCH_CLONING'),
            'RECORD_STUBS': os.getenv('RECORD_STUBS'),
            'NCBI_API_KEY': os.getenv('NCBI_API_KEY'),
            'ALLOWED_ORIGINS': os.getenv('ALLOWED_ORIGINS'),
            'PLANNOTATE_URL': os.getenv('PLANNOTATE_URL'),
            'PLANNOTATE_TIMEOUT': os.getenv('PLANNOTATE_TIMEOUT'),
            'PROXY_URL': os.getenv('PROXY_URL'),
            'PROXY_CERT_FILE': os.getenv('PROXY_CERT_FILE'),
        }

    def tearDown(self):
        # Restore original environment variables
        monkeypatch = pytest.MonkeyPatch()
        for key, value in self.original_env.items():
            if value is None:
                monkeypatch.delenv(key, raising=False)
            else:
                monkeypatch.setenv(key, value)
        reload(os)
        reload(app_settings)

    def test_default_values(self):
        self.assertEqual(app_settings.settings.SERVE_FRONTEND, False)
        self.assertEqual(app_settings.settings.BATCH_CLONING, True)
        self.assertEqual(app_settings.settings.RECORD_STUBS, False)
        # self.assertEqual(app_settings.settings.NCBI_API_KEY, None) > This is different in the CI
        self.assertEqual(app_settings.settings.ALLOWED_ORIGINS, ['http://localhost:3000', 'http://localhost:5173'])
        self.assertEqual(app_settings.settings.PLANNOTATE_URL, None)
        self.assertEqual(app_settings.settings.PLANNOTATE_TIMEOUT, 20)
        self.assertEqual(app_settings.settings.PROXY_URL, None)
        self.assertEqual(app_settings.settings.PROXY_CERT_FILE, None)

    def test_settings_from_env(self):
        monkeypatch = pytest.MonkeyPatch()
        monkeypatch.setenv('SERVE_FRONTEND', '1')
        monkeypatch.setenv('BATCH_CLONING', '0')
        monkeypatch.setenv('RECORD_STUBS', '1')
        monkeypatch.setenv('NCBI_API_KEY', 'test')
        monkeypatch.setenv('ALLOWED_ORIGINS', 'hello,bye')
        monkeypatch.setenv('PLANNOTATE_URL', 'http://dummy/url')
        monkeypatch.setenv('PLANNOTATE_TIMEOUT', '30')
        monkeypatch.setenv('PROXY_URL', 'http://dummy/url')
        monkeypatch.setenv('PROXY_CERT_FILE', 'dummy/cert.pem')

        reload(app_settings)

        self.assertEqual(app_settings.settings.SERVE_FRONTEND, True)
        self.assertEqual(app_settings.settings.BATCH_CLONING, False)
        self.assertEqual(app_settings.settings.RECORD_STUBS, True)
        self.assertEqual(app_settings.settings.NCBI_API_KEY, 'test')
        self.assertEqual(app_settings.settings.ALLOWED_ORIGINS, ['hello', 'bye'])
        self.assertEqual(app_settings.settings.PLANNOTATE_URL, 'http://dummy/url/')  # Trailing slash added
        self.assertEqual(app_settings.settings.PLANNOTATE_TIMEOUT, 30)
        self.assertEqual(app_settings.settings.PROXY_URL, 'http://dummy/url')
        self.assertEqual(app_settings.settings.PROXY_CERT_FILE, 'dummy/cert.pem')

        # Test boolean inputs
        monkeypatch.setenv('SERVE_FRONTEND', 'True')
        reload(app_settings)
        self.assertEqual(app_settings.settings.SERVE_FRONTEND, True)

        monkeypatch.setenv('SERVE_FRONTEND', 'TRUE')
        reload(app_settings)
        self.assertEqual(app_settings.settings.SERVE_FRONTEND, True)

        monkeypatch.setenv('SERVE_FRONTEND', 'true')
        reload(app_settings)
        self.assertEqual(app_settings.settings.SERVE_FRONTEND, True)
