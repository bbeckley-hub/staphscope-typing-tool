# syntax=docker/dockerfile:1.6

FROM mambaorg/micromamba:1.5.8-jammy

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="StaphScope - Advanced Staphylococcus aureus Typing & Lineage Analysis Platform"

USER root

# --- System dependencies ------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        jq \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/staphscope

# --- Copy environment spec FIRST so code edits don't invalidate the env layer -
COPY environment.yml /opt/staphscope/environment.yml

# --- Create the env with a BuildKit cache mount ------------------------------
# No `micromamba clean` here — /opt/conda/pkgs IS the cache and must be kept
# for retry-safe rebuilds. The cache mount survives --no-cache.
RUN --mount=type=cache,target=/opt/conda/pkgs \
    micromamba create -y -n staphscope -f /opt/staphscope/environment.yml

# --- Copy the project ---------------------------------------------------------
COPY . /opt/staphscope/

# --- Editable install of the project itself ----------------------------------
RUN micromamba run -n staphscope pip install --no-deps -e /opt/staphscope

# --- Use the env as the default shell for subsequent RUN commands ------------
SHELL ["micromamba", "run", "-n", "staphscope", "/bin/bash", "-c"]

# --- One-time database setup --------------------------------------------------
# abricate: all bundled resistance/virulence/plasmid databases
RUN abricate --setupdb

# AMRFinderPlus: latest reference DB
RUN cd /opt/staphscope/staphscope/modules/amr_module && \
    python amrfinder_standalone.py --update-db

# S. aureus MLST scheme — baked in so non-root Docker users get a working MLST
RUN staphscope --pull-mlst-db

# --- Runtime environment ------------------------------------------------------
# HOME=/tmp is world-writable so `-u $(id -u):$(id -g)` can write to $HOME
ENV HOME=/tmp
ENV MAMBA_ROOT_PREFIX=/opt/conda

# --- Inline entrypoint --------------------------------------------------------
# Sets up the full conda env activation (PATH, LD_LIBRARY_PATH, CONDA_PREFIX)
# that bioconda binaries like amrfinder, diamond, blastn, and prodigal need
# at runtime. Also prepends `staphscope` if the first argument starts with a
# dash, so `docker run image -i x.fna -o out` works directly.
RUN printf '%s\n' \
    '#!/bin/bash' \
    'set -e' \
    'export CONDA_PREFIX=/opt/conda/envs/staphscope' \
    'export CONDA_DEFAULT_ENV=staphscope' \
    'export PATH=${CONDA_PREFIX}/bin:${PATH}' \
    'export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}' \
    'if [ $# -gt 0 ] && [ "${1#-}" != "$1" ]; then' \
    '    set -- staphscope "$@"' \
    'fi' \
    'exec "$@"' \
    > /usr/local/bin/entrypoint.sh \
    && chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["staphscope", "-h"]
