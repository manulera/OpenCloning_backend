# Backend for OpenCloning
# https://github.com/manulera/OpenCloning_backend

# BUILDER IMAGE
FROM python:3.12-slim-bookworm AS builder
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y && apt-get install -y gcc git g++ wget make

# Build mafft from source
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.525-without-extensions-src.tgz
RUN tar -xzf mafft-7.525-without-extensions-src.tgz
RUN cd mafft-7.525-without-extensions/core && \
    sed -i 's/^PREFIX = .*/PREFIX = \/mafft/' Makefile && \
    make && \
    make install

# FINAL IMAGE
FROM python:3.12-slim-bookworm

# directly output things to stdout/stderr, without buffering
ENV PYTHONUNBUFFERED=1

# create a user to run the app
RUN useradd -ms /bin/bash backend
USER backend
WORKDIR /home/backend

COPY --from=builder /mafft /mafft
# Add mafft to PATH
ENV PATH="/mafft/bin:$PATH"
