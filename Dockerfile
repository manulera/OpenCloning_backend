# Backend for OpenCloning
# https://github.com/manulera/OpenCloning_backend

# BUILDER IMAGE
FROM python:3.12-slim-bookworm AS builder
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y && apt-get install -y gcc git g++ wget build-essential cmake unzip

# Build mafft from source
# mafft creates an entire directory, that's what we copy later and that's why it needs
# to be added to the path
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.525-without-extensions-src.tgz
RUN tar -xzf mafft-7.525-without-extensions-src.tgz
RUN cd mafft-7.525-without-extensions/core && \
    sed -i 's/^PREFIX = .*/PREFIX = \/usr\/local\/bin\/mafft/' Makefile && \
    make && \
    make install

# Build MARS from source
RUN wget https://github.com/manulera/MARS/archive/refs/tags/v0.2.tar.gz && \
tar -xzf v0.2.tar.gz && \
cd MARS-0.2 && \
./pre-install.sh && \
make -f Makefile && \
mv mars /usr/local/bin/mars


RUN useradd -ms /bin/bash backend
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

RUN poetry install --only main

# FINAL IMAGE
FROM python:3.12-slim-bookworm

RUN apt-get update && apt-get upgrade -y && apt-get install --no-install-recommends -y mafft libgomp1 && rm -rf /var/lib/apt/lists/*

# Remove perl, it's not needed and introduced vulnerabilities
RUN apt-get purge --allow-remove-essential -y perl perl-base perl-modules-5.36 && apt-get autoremove -y

# directly output things to stdout/stderr, without buffering
ENV PYTHONUNBUFFERED=1

# create a user to run the app
RUN useradd -ms /bin/bash backend
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
