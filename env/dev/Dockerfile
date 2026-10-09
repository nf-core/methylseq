FROM mambaorg/micromamba:0.25.1
COPY --chown=$MAMBA_USER:$MAMBA_USER conda_environment.yml /tmp/conda_environment.yml
RUN micromamba install -y -n base -f /tmp/conda_environment.yml && \
    micromamba clean --all --yes
USER root
RUN apt-get update && apt-get -y install sudo
RUN sudo apt-get update -qq && \
    sudo apt-get install -y zip procps ghostscript

WORKDIR /work
# COPY --chown=$MAMBA_USER:$MAMBA_USER .aws/ /home/mambauser/.aws/