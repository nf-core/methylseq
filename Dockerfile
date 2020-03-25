FROM nfcore/base:1.7
LABEL authors="Phil Ewels" \
      description="Docker image containing all software requirements for the nf-core/methylseq pipeline"
RUN     apt-get update -y && \
        apt-get install -y --no-install-recommends apt-utils && \
        apt-get install zlib1g-dev -y && \
        apt-get install libbz2-dev -y && \
        apt-get install liblzma-dev -y && \
        apt-get install libncurses5-dev -y && \
        apt-get install curl -y 

RUN cd / && \
        curl -OL $(curl -s https://api.github.com/repos/zwdzwd/biscuit/releases/latest | grep browser_download_url | grep linux_amd64 | cut -d '"' -f 4) && \
        chmod 755 biscuit*linux_amd64 && \
        mv biscuit*linux_amd64 biscuit

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-methylseq-1.4.1 > nf-core-methylseq-1.4.1.yml
ENV PATH /opt/conda/envs/nf-core-methylseq-1.4.1/bin:$PATH
