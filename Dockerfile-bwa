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

# Install PicardTools
RUN curl -fsSL https://github.com/broadinstitute/picard/releases/download/2.0.1/picard-tools-2.0.1.zip -o /opt/picard-tools-2.0.1.zip && \
    unzip /opt/picard-tools-2.0.1.zip -d /opt/ && \
    rm /opt/picard-tools-2.0.1.zip
ENV PICARD_HOME /opt/picard-tools-2.0.1

# Install Qualimap
RUN mkdir /opt/Qualimap && \
    curl -fsSL https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip -o /opt/Qualimap/qualimap.zip && \
    unzip /opt/Qualimap/qualimap.zip -d /opt/Qualimap && \
    ln -s /opt/Qualimap/qualimap_v2.2.1/qualimap /usr/local/bin/qualimap && \
    rm /opt/Qualimap/qualimap.zip

# Install BWA
RUN mkdir /opt/bwa && \
    curl -fsSL https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2 -o /opt/bwa/bwa.tar.bz2 && \
    tar xvjf /opt/bwa/bwa.tar.bz2 -C /opt/bwa/ && \
    cd /opt/bwa/bwa-0.7.15/ && \
    make && \
    ln -s /opt/bwa/bwa-0.7.15/bwa /usr/local/bin/bwa && \
    rm /opt/bwa/bwa.tar.bz2

# Install bwa-meth
RUN pip install toolshed && \
    pip install git+git://github.com/brentp/bwa-meth.git

# Install MethylDackel
RUN mkdir /opt/MethylDackel && \
    curl -fsSL https://github.com/dpryan79/MethylDackel/archive/0.2.1.zip -o /opt/MethylDackel/MethylDackel.zip && \
    unzip /opt/MethylDackel/MethylDackel.zip -d /opt/MethylDackel && \
    cd /opt/MethylDackel/MethylDackel-0.2.1 && \
    make && \
    make install prefix=/usr/local/bin && \
    rm /opt/MethylDackel/MethylDackel.zip

# Install MultiQC
RUN pip install multiqc

# Create root directories for UPPMAX and c3se hebbe
RUN mkdir /pica /lupus /crex1 /crex2 /proj /scratch /sw \
          /c3se /local /apps
