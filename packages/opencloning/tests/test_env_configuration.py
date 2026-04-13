import unittest
from pytest import MonkeyPatch
from importlib import reload
from fastapi.testclient import TestClient
from opencloning.utils import TemporaryFolderOverride
import shutil
import glob
import os
import tempfile
import opencloning.main as main
import opencloning.app_settings as app_settings


test_files = os.path.join(os.path.dirname(__file__), 'test_files')


class TestTemporaryFolderOverride(unittest.TestCase):

    def test_temporary_folder_override(self):
        # File that did not exist before
        with tempfile.TemporaryDirectory() as tempdir:
            override_folder = os.path.join(tempdir, 'dummy_test_file_folder_override')
            with TemporaryFolderOverride(override_folder):
                assert os.path.exists(override_folder)
            assert not os.path.exists(override_folder)

        # File that existed before
        with tempfile.TemporaryDirectory() as tempdir:
            override_folder = os.path.join(tempdir, 'dummy_test_file_folder_override')
            os.mkdir(override_folder)
            # Add a file to the folder
            with open(os.path.join(override_folder, 'dummy_file.txt'), 'w') as f:
                f.write('Dummy file')
            with TemporaryFolderOverride(override_folder):
                assert os.path.exists(override_folder)
                assert not os.path.exists(os.path.join(override_folder, 'dummy_file.txt'))

            assert os.path.exists(override_folder)
            assert os.path.exists(os.path.join(override_folder, 'dummy_file.txt'))


class TestServeFrontend(unittest.TestCase):
    # DO this before each test
    def setUp(self):

        self.folder_override = TemporaryFolderOverride('frontend')
        self.folder_override.__enter__()

        for f in glob.glob(f'{test_files}/dummy_frontend/*'):
            if os.path.isdir(f):
                shutil.copytree(f, f'./frontend/{f.split("/")[-1]}', dirs_exist_ok=True)
            else:
                shutil.copy2(f, './frontend')

        # Has to be imported here to get the right environment variable
        MonkeyPatch().setenv('SERVE_FRONTEND', '1')
        reload(app_settings)
        reload(main)
        client = TestClient(main._app)
        self.client = client

    # DO this after each test
    def tearDown(self):
        self.folder_override.__exit__(None, None, None)
        MonkeyPatch().setenv('SERVE_FRONTEND', '0')
        reload(app_settings)
        reload(main)

    def test_serve_frontend(self):
        # The index is served at the root
        response = self.client.get('/')
        self.assertEqual(response.status_code, 200)
        self.assertIn('<noscript>You need to enable JavaScript to run this app.</noscript>', response.text)

        # The rest of files can be accessed:
        for file in ['config.json', 'favicon.ico', 'robots.txt', 'logo192.png', 'logo512.png', 'manifest.json']:
            response = self.client.get(f'/{file}')
            self.assertEqual(response.status_code, 200)

        # Can still make a normal request (e.g. enzymes)
        response = self.client.get('/restriction_enzyme_list')
        self.assertEqual(response.status_code, 200)
        self.assertIn('EcoRI', response.text)

        # If requesting a file that does not exist, it should return a 404
        response = self.client.get('/dummy_file.json')
        self.assertEqual(response.status_code, 404)

    def test_frontend_config_endpoint_uses_env_vars(self):
        monkeypatch = MonkeyPatch()
        monkeypatch.setenv('BACKEND_URL', 'https://api.example.org/')
        monkeypatch.setenv('DATABASE', 'production')
        monkeypatch.setenv('SHOW_APP_BAR', '0')
        monkeypatch.setenv('NO_EXTERNAL_REQUESTS', '1')
        monkeypatch.setenv('ENABLE_ASSEMBLER', 'False')
        monkeypatch.setenv('ENABLE_PLANNOTATE', 'false')
        monkeypatch.setenv('STATIC_CONTENT_PATH', '/srv/static')
        reload(app_settings)
        reload(main)
        client = TestClient(main._app)

        response = client.get('/config.json')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            response.json(),
            {
                'backendUrl': 'https://api.example.org/',
                'database': 'production',
                'showAppBar': False,
                'noExternalRequests': True,
                'enableAssembler': False,
                'enablePlannotate': False,
                'staticContentPath': '/srv/static',
            },
        )
        monkeypatch.undo()
