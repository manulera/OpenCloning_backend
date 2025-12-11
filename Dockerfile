# Backend for OpenCloning
# https://github.com/manulera/OpenCloning_backend

# BUILDER IMAGE
FROM python:3.12-alpine3.21 AS builder

RUN apk update --no-cache && apk add --no-cache build-base bash cmake git

# Build mafft from source
RUN git clone https://gitlab.com/manulera/mafft.git
RUN cd mafft && git checkout 3cf6cfa8b81756b27eed28ed368474b9ee50e2da
RUN cd mafft/core && \
    sed -i 's/^PREFIX = .*/PREFIX = \/usr\/local\/bin\/mafft/' Makefile && \
    make && \
    make install

# Build MARS from source
RUN wget https://github.com/manulera/MARS/archive/refs/tags/v0.2.6.tar.gz && \
tar -xzf v0.2.6.tar.gz && \
cd MARS-0.2.6 && \
./pre-install.sh && \
make -f Makefile && \
mv mars /usr/local/bin/mars


RUN adduser -s /bin/bash -D backend
USER backend
WORKDIR /home/backend

ENV PIP_DISABLE_PIP_VERSION_CHECK=on

# Poetry
# https://python-poetry.org/docs/configuration/#using-environment-variables
# make poetry install to this location
ENV POETRY_HOME="/home/backend/.bin"
# do not ask any interactive question
ENV POETRY_NO_INTERACTION=1
# never create virtual environment automatically, only use env prepared by us
ENV POETRY_VIRTUALENVS_CREATE=false
# this is where our dependencies and virtual environment will live
ENV VIRTUAL_ENV="/home/backend/venv"
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$POETRY_HOME/bin:$VIRTUAL_ENV/bin:$PATH"

RUN pip install --no-cache-dir poetry

COPY ./src ./src
COPY ./poetry.lock .
COPY ./pyproject.toml .
COPY ./README.md .

# Set version in pyproject.toml before installing
ARG PACKAGE_VERSION="0.1.0"
RUN sed -i "s/^version = .*/version = \"${PACKAGE_VERSION}\"/" pyproject.toml

## Docker-test-comment
# The above comment is used to build another Dockerfile to run the tests in the container during CI.

RUN poetry install --only main

# FINAL IMAGE
FROM python:3.12-alpine3.21

# You need bash to run mafft and runtime libraries for MARS
RUN apk update --no-cache && apk add --no-cache bash libstdc++ libgomp libgcc

# directly output things to stdout/stderr, without buffering
ENV PYTHONUNBUFFERED=1

# create a user to run the app
RUN adduser -s /bin/bash -D backend
USER backend
WORKDIR /home/backend

ENV VIRTUAL_ENV="/home/backend/venv"
COPY --from=builder $VIRTUAL_ENV $VIRTUAL_ENV
COPY --from=builder /usr/local/bin/mars /usr/local/bin/mars
COPY --from=builder /usr/local/bin/mafft /usr/local/bin/mafft

COPY ./src ./src
ENV PATH="/usr/local/bin/mafft/bin:$VIRTUAL_ENV/bin:$PATH"
# For example, ROOT_PATH="/syc"
ENV ROOT_PATH=""
ENV USE_HTTPS=false

COPY ./docker_entrypoint.sh ./docker_entrypoint.sh

CMD ["sh", "./docker_entrypoint.sh"]
