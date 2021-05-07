FROM nfcore/base:1.13.3
LABEL authors="Phil Ewels" \
      description="Docker image containing all software requirements for the nf-core/methylseq pipeline"

# Install libtbb system dependency for bowtie2
RUN apt-get update \
      && apt-get install -y libtbb-dev \
      && apt-get clean -y \
      && rm -rf /var/lib/apt/lists/*

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-methylseq-1.6.1/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-methylseq-1.6.1 > nf-core-methylseq-1.6.1.yml
