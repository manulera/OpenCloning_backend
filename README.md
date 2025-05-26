[![Python tests](https://github.com/manulera/OpenCloning_backend/actions/workflows/ci.yml/badge.svg)](https://github.com/manulera/OpenCloning_backend/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/manulera/OpenCloning_backend/graph/badge.svg?token=CFIB2H6WMO)](https://codecov.io/gh/manulera/OpenCloning_backend)

# OpenCloning Backend API

This API is part of a bigger application, before going further, please go to the [main project readme](https://github.com/manulera/OpenCloning), where you can find an introduction.

This python API is built with [FastAPI](https://fastapi.tiangolo.com/) and is for *in silico* cloning.

## Summary

Read [main project readme](https://github.com/manulera/OpenCloning) first.

This API provides a series of entry points. The API documentation can be accessed [here](https://api.opencloning.org/docs). You can use the documentation page to try some request directly on the browser. Otherwise, the API is open for you to make requests from a python script or command line at: [https://api.opencloning.org/](https://api.opencloning.org/).

## Scripting

The API functions can also be used to write python scripts to automate cloning. See the [scripting examples](examples/scripting) for more information.

## Migrating between model versions and fixing model bugs

* The data model changes, so the json files you created may not be compatible with the newest version of the library, which uses the latest data mode. You can easily fix this using `python -m opencloning_linkml.migrations.migrate file.json
` see [full documentation](https://github.com/OpenCloning/OpenCloning_LinkML?tab=readme-ov-file#migration-from-previous-versions-of-the-schema).
* Before version 0.3, there was a bug for assembly fields that included locations spanning the origin. See the details and how to fix it in the documentation of [this file](./src/opencloning/bug_fixing/README.md).

## Getting started

If you want to quickly set up a local instance of the frontend and backend of the application, check [getting started in 5 minutes](https://github.com/manulera/OpenCloning#timer_clock-getting-started-in-5-minutes) in the main repository.

### Running locally

You can install this as a python package:

```bash
# Create a virtual environment
python -m venv .venv
# Activate the virtual environment
source .venv/bin/activate
# Install the package from pypi
pip install opencloning
# Run the API (uvicorn should be installed in the virtual environment)
uvicorn opencloning.main:app
```

### Running locally if you want to contribute

For the management of the dependencies `poetry` is used, if you don't have it, visit https://python-poetry.org/.

In the project directory:

```bash
# This should install the dependencies and create a virtual environment
poetry install

# Install the pre-commit hooks
pre-commit install

# Activate the virtual environment
poetry shell

```

The virtual environment is installed in the project folder. This is convenient if you are using an IDE for development. For settings of vscode see the folder `.vscode`.

Now you should be able to run the api by running:

```bash
# The --reload argument will reload the API if you make changes to the code
uvicorn opencloning.main:app --reload --reload-exclude='.venv'
```

Then you should be able to open the API docs at [http://127.0.0.1:8000/docs](http://127.0.0.1:8000/docs) to know that your API is working.

### Running locally with docker :whale:

If you want to serve the full site (backend and frontend) with docker, check [getting started in 5 minutes](https://github.com/manulera/OpenCloning#timer_clock-getting-started-in-5-minutes) in the main repository.

If you want to serve only the backend from a docker container, an image is available at [manulera/opencloningbackend](https://hub.docker.com/r/manulera/opencloningbackend). The image is built from the Dockerfile in the root of this repository and exposes the port 3000. To run it:

```bash
docker build -t manulera/opencloningbackend .
docker run -d --name backendcontainer -p 8000:8000 manulera/opencloningbackend

```

If you don't want to download the repository and build the image, you can fetch the latest image from dockerhub (same image that is used in [https://api.opencloning.org/](https://api.opencloning.org/))

```bash
docker pull manulera/opencloningbackend
docker run -d --name backendcontainer -p 8000:8000 manulera/opencloningbackend
```

The api will be running at `http://localhost:8000`, so you should be able to access the docs at [http://localhost:8000/docs](http://localhost:8000/docs).

### Connecting to the frontend

If you want to receive requests from the [frontend](https://github.com/manulera/OpenCloning_frontend), or from another web application you may have to include the url of the frontend application in the CORS exceptions. By default, if you run the dev server with `uvicorn opencloning.main:app --reload --reload-exclude='.venv'`, the backend will accept requests coming from `http://localhost:3000`, which is the default address of the frontend dev server (ran with `yarn start`).

If you want to change the allowed origins, you can do so via env variables (comma-separated). e.g.:

```
ALLOWED_ORIGINS=http://localhost:3000,http://localhost:3001 uvicorn opencloning.main:app --reload --reload-exclude='.venv'
```

Similarly, the frontend should be configured to send requests to the backend address, [see here](https://github.com/manulera/OpenCloning_frontend#connecting-to-the-backend).

#### Serving the frontend from the backend

You may prefer to handle everything from a single server. You can do so by:
* Build the [frontend](https://github.com/manulera/OpenCloning_frontend) with `yarn build`.
* Copy the folder `build` from the frontend to the root directory of the backend, and rename it to `frontend`.
* Set the environment variable `SERVE_FRONTEND=1` when running the backend. By default this will remove all allowed origins, but you can still set them with `ALLOWED_ORIGINS`.
* Set the value of `backendUrl` in `frontend/config.js` to `/`.
* Now, when you go to the root of the backend (e.g. `http://localhost:8000`), you should receive the frontend instead of the greeting page of the API.

You can see how this is done in this [docker image](https://github.com/manulera/OpenCloning/blob/master/Dockerfile) and [docker-compose file](https://github.com/manulera/OpenCloning/blob/master/docker-compose.yml).

## Contributing :hammer_and_wrench:

Check [contribution guidelines in the main repository](https://github.com/manulera/OpenCloning/blob/master/CONTRIBUTING.md) for general guidelines.

For more specific tasks:
* Creating a new type of source: follow the [new source issue template](.github/ISSUE_TEMPLATE/new-source.md). You can create an issue like that [here](https://github.com/manulera/OpenCloning_backend/issues/new?assignees=&labels=new-source&projects=&template=new-source.md&title=New+source%3A+%3Cname-of-source%3E).

## Running the tests locally

```
pytest -v -ks
```

## Notes

### Ping a particular library version from github:

```
poetry add git+https://github.com/BjornFJohansson/pydna#4fd760d075f77cceeb27969e017e04b42f6d0aa3
```

When installing the last version, sometimes poetry may not be able to access the latest version
and complain with `Could not find a matching version`, in that case you can do:

```
poetry cache clear pypi --all
```

### Generating API stubs

For the frontend, it may be useful to produce stubs (I use them for writing the tests). See how this is implemented
by looking at the `RecordStubRoute` class in `api_config_utils.py`. To run the dev server and record stubs:

```bash
RECORD_STUBS=1 uvicorn opencloning.main:app --reload --reload-exclude='.venv'
```

This will record the stubs (requests and responses) in the `stubs` folder.
