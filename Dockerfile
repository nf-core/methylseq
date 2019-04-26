FROM nfcore/base
LABEL authors="phil.ewels@scilifelab.se" \
    description="Docker image containing all requirements for the nf-core/methylseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-methylseq-1.4dev/bin:$PATH
