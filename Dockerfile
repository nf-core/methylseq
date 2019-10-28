FROM nfcore/base:dev
LABEL authors="Phil Ewels" \
      description="Docker image containing all requirements for nf-core/methylseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-methylseq-1.4dev > nf-core-methylseq-1.4dev.yml
ENV PATH /opt/conda/envs/nf-core-methylseq-1.4dev/bin:$PATH
