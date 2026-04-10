# Backend for OpenCloning
# https://github.com/manulera/OpenCloning_backend

# BUILDER IMAGE
FROM manulera/opencloningbackend-base:python_3.12-alpine3.21 AS builder

RUN adduser -s /bin/bash -D backend
USER backend
WORKDIR /home/backend

ENV PIP_DISABLE_PIP_VERSION_CHECK=on
ENV UV_COMPILE_BYTECODE=1

ENV VIRTUAL_ENV="/home/backend/venv"
ENV UV_PROJECT_ENVIRONMENT=$VIRTUAL_ENV
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN pip install --no-cache-dir uv

COPY pyproject.toml uv.lock ./
COPY packages/opencloning/pyproject.toml packages/opencloning/README.md packages/opencloning/
COPY packages/opencloning/src packages/opencloning/src

ARG PACKAGE_VERSION="0.1.0"
ENV SETUPTOOLS_SCM_PRETEND_VERSION="${PACKAGE_VERSION}"

## Docker-test-comment
# The above comment is used to build another Dockerfile to run the tests in the container during CI.

RUN uv sync --frozen --package opencloning --no-default-groups --no-editable

# FINAL IMAGE
FROM python:3.12-alpine3.21

# You need bash to run mafft and runtime libraries for MARS
RUN apk update --no-cache && apk add --no-cache bash libstdc++ libgomp libgcc

# create a user to run the app
RUN adduser -s /bin/bash -D backend
USER backend
WORKDIR /home/backend

ENV VIRTUAL_ENV="/home/backend/venv"
COPY --from=builder $VIRTUAL_ENV $VIRTUAL_ENV
COPY --from=builder /usr/local/bin/mars /usr/local/bin/mars
COPY --from=builder /usr/local/bin/mafft /usr/local/bin/mafft

ENV PATH="/usr/local/bin/mafft/bin:$VIRTUAL_ENV/bin:$PATH"
# For example, ROOT_PATH="/syc"
ENV ROOT_PATH=""
ENV USE_HTTPS=false

COPY ./docker_entrypoint.sh ./docker_entrypoint.sh

CMD ["sh", "./docker_entrypoint.sh"]
