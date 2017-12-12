# openjdk:8 moved from debian jessie to stretch after u131, which breaks everything (bowtie)
FROM openjdk:8u121

LABEL authors="phil.ewels@scilifelab.se,denis.moreno@scilifelab.se" \
    description="Docker image containing all requirements for NGI-MethylSeq pipeline"

# Install container-wide requrements gcc, pip, zlib, libssl, make, libncurses, fortran77, g++, R
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        g++ \
        gcc \
        gfortran \
        libbz2-dev \
        libcurl4-openssl-dev \
        libgsl0-dev \
        liblzma-dev \
        libncurses5-dev \
        libpcre3-dev \
        libreadline-dev \
        libssl-dev \
        libtbb-dev \
        make \
        python-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pip
RUN curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py && \
    python /opt/get-pip.py && \
    rm /opt/get-pip.py

RUN curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -o /opt/fastqc_v0.11.5.zip && \
    unzip /opt/fastqc_v0.11.5.zip -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm /opt/fastqc_v0.11.5.zip

# Install cutadapt
RUN pip install cutadapt

# Install TrimGalore
RUN mkdir /opt/TrimGalore && \
    curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.4.zip -o /opt/TrimGalore/0.4.4.zip && \
    unzip /opt/TrimGalore/0.4.4.zip -d /opt/TrimGalore && \
    ln -s /opt/TrimGalore/TrimGalore-0.4.4/trim_galore /usr/local/bin/trim_galore && \
    rm /opt/TrimGalore/0.4.4.zip

# Install SAMTools
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 -o /opt/samtools-1.5.tar.bz2 && \
    tar xvjf /opt/samtools-1.5.tar.bz2 -C /opt/ && \
    cd /opt/samtools-1.5 && \
    make && \
    make install && \
    rm /opt/samtools-1.5.tar.bz2

# Install Bowtie2
RUN mkdir /opt/bowtie2 && \
    curl -fsSL https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.2/bowtie2-2.3.2-source.zip -o /opt/bowtie2/bowtie2.zip && \
    unzip /opt/bowtie2/bowtie2.zip -d /opt/bowtie2 && \
    cd /opt/bowtie2/bowtie2-2.3.2/ && \
    make && \
    ln -s /opt/bowtie2/bowtie2-2.3.2/bowtie2 /usr/local/bin/bowtie2 && \
    ln -s /opt/bowtie2/bowtie2-2.3.2/bowtie2-build /usr/local/bin/bowtie2-build && \
    rm /opt/bowtie2/bowtie2.zip

# Install Bismark
RUN mkdir /opt/Bismark && \
    curl -fsSL https://github.com/FelixKrueger/Bismark/archive/0.18.2.zip -o /opt/Bismark/bismark.zip && \
    unzip /opt/Bismark/bismark.zip -d /opt/Bismark && \
    rm /opt/Bismark/bismark.zip
ENV PATH="/opt/Bismark/Bismark-0.18.2:${PATH}"

# Install Qualimap
RUN mkdir /opt/Qualimap && \
    curl -fsSL https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip -o /opt/Qualimap/qualimap.zip && \
    unzip /opt/Qualimap/qualimap.zip -d /opt/Qualimap && \
    ln -s /opt/Qualimap/qualimap_v2.2.1/qualimap /usr/local/bin/qualimap && \
    rm /opt/Qualimap/qualimap.zip

# Install MultiQC
RUN pip install git+https://github.com/ewels/MultiQC.git

# Create root directories for UPPMAX and c3se hebbe
RUN mkdir /pica /lupus /crex1 /crex2 /proj /scratch /sw \
          /c3se /local /apps
