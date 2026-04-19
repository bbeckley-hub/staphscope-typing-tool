FROM continuumio/miniconda3:latest

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="StaphScope - Advanced Staphylococcus aureus Typing & Lineage Analysis Platform"

# Install system dependencies (procps for resource monitoring, jq for JSON parsing, git for version info)
RUN apt-get update && apt-get install -y \
    procps \
    jq \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/staphscope

# Copy entire project
COPY . /opt/staphscope/

# Create the Conda environment from environment.yml
RUN conda env create -f environment.yml && \
    conda clean -afy

# Make the environment the default for RUN commands
SHELL ["conda", "run", "-n", "staphscope", "/bin/bash", "-c"]

# Run abricate database setup (one-time)
RUN abricate --setupdb

# Run AMR database setup (downloads latest AMRfinderPlus database)
RUN cd /opt/staphscope/staphscope/modules/amr_module && \
    python amrfinder_standalone.py --update-db

# Set entrypoint to staphscope command
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "staphscope", "staphscope"]
CMD ["-h"]
