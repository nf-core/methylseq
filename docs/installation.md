# nf-core/methylseq installation

To start using the nf-core/methylseq pipeline, follow the steps below:

1. [Install Nextflow](#1-install-nextflow)
2. [Install the pipeline](#2-install-the-pipeline)
    * [Automatic](#21-automatic)
    * [Offline](#22-offline)
    * [Development](#23-development)
3. [Pipeline configuration](#3-pipeline configuration)
    * [Software deps: Docker and Singularity](#31-software-deps-docker-and-singularity)
    * [Software deps: Bioconda](#32-software-deps-bioconda)
    * [Configuration profiles](#33-configuration-profiles)
4. [Reference genomes](#4-reference-genomes)
5. [Appendices](#appendices)
    * [Running on UPPMAX](#running-on-uppmax)

## 1) Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v8+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

See [nextflow.io](https://www.nextflow.io/) for further instructions on how to install and configure Nextflow.

## 2) Install the pipeline

#### 2.1) Automatic
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `nf-core/methylseq` is specified as the pipeline name.

#### 2.2) Offline
The above method requires an internet connection so that Nextflow can download the pipeline files. If you're running on a system that has no internet connection, you'll need to download and transfer the pipeline files manually:

```bash
#TODO: Find correct URL and command
curl https://github.com/nf-core/methylseq.git
nextflow run /path/to/pipelines/nf-core-methylseq
```

#### 2.3) Development

If you would like to make changes to the pipeline, it's best to make a fork on GitHub and then clone the files. Once cloned you can run the pipeline directly as above.


## 3) Pipeline configuration
By default, the pipeline runs with the `standard` configuration profile. This uses a number of sensible defaults for process requirements and is suitable for running on a simple (if powerful!) basic server. You can see this configuration in [`conf/base.config`](../conf/base.config).

Be warned of two important points about this default configuration:

1. The default profile uses the `local` executor
    * All jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.
2. Nextflow will expect all software to be installed and available on the `PATH`

#### 3.1) Software deps: Docker and Singularity
Running the pipeline with the option `-with-singularity` or `-with-docker` tells Nextflow to enable either [Singularity](http://singularity.lbl.gov/) or Docker for this run. An image containing all of the software requirements will be automatically fetched and used (https://hub.docker.com/r/nf-core/methylseq).

If running offline with Singularity, you'll need to download and transfer the Singularity image first:

```bash
# Switch 1.4 with the version of the pipeline that you're using
singularity pull --name nfcore-methylseq-1.4.simg shub://nfcore/methylseq:1.4
```

Once transferred, use `-with-singularity` but specify the path to the image file:

```
nextflow run /path/to/nf-core-methylseq -with-singularity /path/to/nfcore-methylseq-1.4.simg
```

#### 3.2) Software deps: bioconda

If you're unable to use either Docker or Singularity but you have conda installed, you can use the bioconda environment that comes with the pipeline. Running this command will create a new conda environment with all of the required software installed:

```bash
conda env create -f environment.yml
conda clean -a # Recommended, not essential
source activate nfcore-methylseq
```

The `environment.yml` file is the one that is packaged with the pipeline (you may need to download this file).


#### 3.3) Configuration profiles

Nextflow can be configured to run on a wide range of different computational infrastructures. In addition to the above pipeline-specific parameters it is likely that you will need to define system-specific options. For more information, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/).

Whilst most parameters can be specified on the command line, it is usually sensible to create a configuration file for your environment.

If you are the only person to be running this pipeline, you can create your config file as `~/.nextflow/config` and it will be applied every time you run Nextflow. Alternatively, save the file anywhere and reference it when running the pipeline with `-c path/to/config`.

If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `-profile <name>` when running the pipeline.

The pipeline comes with several such config profiles - see the installation appendices and usage documentation for more information.

## 4) Reference Genomes
The nf-core/methylseq pipeline needs a reference genome for read alignment. Support for many common genomes is built in if running on UPPMAX or AWS, by using [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

If you don't want to use the illumina iGenomes you can supply either a Bismark / BWA reference or a FASTA file. This can be done on the command line (see the [usage docs](usage.md#)).
Alternatively, you can add the paths to your NextFlow config under a relevant id and just specify this id with `--genome ID` when you run the pipeline:

```groovy
params {
  genomes {
    'YOUR-ID' {
      bismark  = '<PATH TO BISMARK REF>/BismarkIndex'
      bwa_meth  = '<PATH TO BWAMETH REF>/genome.fa'
      fasta  = '<PATH TO FASTA FILE>/genome.fa' // used if above is not specified
    }
    'OTHER-GENOME' {
      // [..]
    }
  }
  // Optional - default genome. Ignored if --genome 'OTHER-GENOME' specified on command line
  genome = 'YOUR-ID'
}
```


## Appendices

#### Running on UPPMAX
To run the pipeline on the [Swedish UPPMAX](https://www.uppmax.uu.se/) clusters (`rackham`, `irma`, `bianca` etc), use the command line flag `-profile uppmax`. This tells Nextflow to submit jobs using the SLURM job executor with Singularity for software dependencies.

Note that you will need to specify your UPPMAX project ID when running a pipeline. To do this, use the command line flag `--project <project_ID>`. The pipeline will exit with an error message if you try to run it pipeline with the default UPPMAX config profile without a project.

**Optional Extra:** To avoid having to specify your project every time you run Nextflow, you can add it to your personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```groovy
params.project = 'project_ID' // eg. b2017123
```
