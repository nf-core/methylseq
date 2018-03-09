FROM continuumio/miniconda
LABEL authors="phil.ewels@scilifelab.se" \
    description="Docker image containing all requirements for the nf-core/MethylSeq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml
ENV PATH /opt/conda/envs/nfcore-methylseq/bin:$PATH
