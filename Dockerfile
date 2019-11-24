FROM nfcore/base:1.7
LABEL authors="Phil Ewels" \
      description="Docker image containing all software requirements for the nf-core/methylseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-methylseq-1.5dev > nf-core-methylseq-1.5dev.yml
ENV PATH /opt/conda/envs/nf-core-methylseq-1.5dev/bin:$PATH
