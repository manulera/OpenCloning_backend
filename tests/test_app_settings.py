import unittest
import pytest
from importlib import reload
import os

from opencloning import app_settings


class TestAppSettings(unittest.TestCase):
    def setUp(self):
        self.monkeypatch = pytest.MonkeyPatch()
        self.monkeypatch.delenv('ADDGENE_USERNAME', raising=False)
        self.monkeypatch.delenv('ADDGENE_PASSWORD', raising=False)
        reload(app_settings)

    def tearDown(self):
        self.monkeypatch.undo()
        reload(os)
        reload(app_settings)

    def test_default_values(self):
        self.assertEqual(app_settings.settings.SERVE_FRONTEND, False)
        self.assertEqual(app_settings.settings.BATCH_CLONING, True)
        self.assertEqual(app_settings.settings.RECORD_STUBS, False)
        # self.assertEqual(app_settings.settings.NCBI_API_KEY, None) > This is different in the CI
        self.assertEqual(
            app_settings.settings.ALLOWED_ORIGINS,
            ['http://localhost:3000', 'http://localhost:5173', 'http://localhost:3002'],
        )
        self.assertEqual(app_settings.settings.PLANNOTATE_URL, None)
        self.assertEqual(app_settings.settings.PLANNOTATE_TIMEOUT, 20)
        self.assertEqual(app_settings.settings.PROXY_URL, None)
        self.assertEqual(app_settings.settings.PROXY_CERT_FILE, None)
        self.assertEqual(app_settings.settings.ADDGENE_USERNAME, None)
        self.assertEqual(app_settings.settings.ADDGENE_PASSWORD, None)

    def test_settings_from_env(self):
        self.monkeypatch = pytest.MonkeyPatch()
        self.monkeypatch.setenv('SERVE_FRONTEND', '1')
        self.monkeypatch.setenv('BATCH_CLONING', '0')
        self.monkeypatch.setenv('RECORD_STUBS', '1')
        self.monkeypatch.setenv('NCBI_API_KEY', 'test')
        self.monkeypatch.setenv('ALLOWED_ORIGINS', 'hello,bye')
        self.monkeypatch.setenv('PLANNOTATE_URL', 'http://dummy/url')
        self.monkeypatch.setenv('PLANNOTATE_TIMEOUT', '30')
        self.monkeypatch.setenv('PROXY_URL', 'http://dummy/url')
        self.monkeypatch.setenv('PROXY_CERT_FILE', 'dummy/cert.pem')
        self.monkeypatch.setenv('ADDGENE_USERNAME', 'dummy-user')
        self.monkeypatch.setenv('ADDGENE_PASSWORD', 'dummy-password')

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
        self.assertEqual(app_settings.settings.ADDGENE_USERNAME, 'dummy-user')
        self.assertEqual(app_settings.settings.ADDGENE_PASSWORD, 'dummy-password')

        # Test boolean inputs
        self.monkeypatch.setenv('SERVE_FRONTEND', 'True')
        reload(app_settings)
        self.assertEqual(app_settings.settings.SERVE_FRONTEND, True)

        self.monkeypatch.setenv('SERVE_FRONTEND', 'TRUE')
        reload(app_settings)
        self.assertEqual(app_settings.settings.SERVE_FRONTEND, True)

        self.monkeypatch.setenv('SERVE_FRONTEND', 'true')
        reload(app_settings)
        self.assertEqual(app_settings.settings.SERVE_FRONTEND, True)
