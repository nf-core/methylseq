From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Phil Ewels <phil.ewels@scilifelab.se>
    DESCRIPTION Docker image containing all requirements for the nf-core/methylseq pipeline
    VERSION 1.1dev

%environment
    PATH=/opt/conda/envs/nfcore-methylseq-1.1/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
