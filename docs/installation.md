# nf-core/methylseq Installation

To start using the nf-core/methylseq pipeline, there are three steps described below:

1. [Install Nextflow](#1-install-nextflow)
2. [Install the pipeline](#2-install-the-pipeline)
3. Configure the pipeline
    * [Swedish UPPMAX System](#31-configuration-uppmax)
    * [Other Clusters](#32-configuration-other-clusters)
    * [Docker](#33-configuration-docker)
    * [Amazon AWS](#34-configuration-amazon-ec2)

## 1) Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v7+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

See [nextflow.io](https://www.nextflow.io/) and [NGI-NextflowDocs](https://github.com/SciLifeLab/NGI-NextflowDocs) for further instructions on how to install and configure Nextflow.

## 2) Install the Pipeline
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `nf-core/methylseq` is specified as the pipeline name.

If you prefer, you can download the files yourself from GitHub and run them directly:

```bash
git clone https://github.com/nf-core/methylseq.git
nextflow run nf-core/methylseq/bismark.nf
nextflow run nf-core/methylseq/bwa-meth.nf # Alternative bwa-meth pipeline
```

## 3) Pipeline Configuration
By default, the pipeline runs with the `standard` configuration profile. This uses a number of sensible defaults for process requirements and is suitable for running on a simple (if powerful!) basic server. You can see this configuration in [`conf/base.config`](../conf/base.config).

Be warned of two important points about this default configuration:

1. The default profile uses the `local` executor
    * All jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.
2. Nextflow will expect all software to be installed and available on the `PATH`

The default resource requests are all measured against the following default limits:

* `max_memory` = 128.GB
* `max_cpus` = 16
* `max_time` = 240.h

To adjust these limits, specify them on the command line with two hyphens, eg. `--max_memory '64.GB'`.

Note that these limits are the maximum to be used _per task_. Nextflow will automatically attempt to parallelise as many jobs as possible given the available resources.

## 3.2) Configuration: Docker and Singularity
Running the pipeline with the option `-with-docker` tells Nextflow to enable docker for this run. The above base profile already specifies the container to be used (hub.docker.com/r/nf-core/methylseq).

If you don't have Docker available and do have [Singularity](http://singularity.lbl.gov/), instead use `-with-singularity` on the command line.

## 3.3) Configuration: Amazon Cloud
The pipeline comes with a configuration profile to be used with Amazon AWS. To use it, run the pipeline with `-profile aws`. This tells Nextflow to use the `ignite` job executor and Docker. Note that there are very many different ways to run Nextflow on AWS. This is beyond the scope of this documentation, please see the main Nextflow documentation.

## 3.4) Configuration: UPPMAX
To run the pipeline on the [Swedish UPPMAX](https://www.uppmax.uu.se/) clusters (`rackham`, `irma`, `bianca` etc), use the command line flag `-profile uppmax`. This tells Nextflow to submit jobs using the SLURM job executor with Singularity for software dependencies.

Note that you will need to specify your UPPMAX project ID when running a pipeline. To do this, use the command line flag `--project <project_ID>`. The pipeline will exit with an error message if you try to run it pipeline with the default UPPMAX config profile without a project.

**Optional Extra:** To avoid having to specify your project every time you run Nextflow, you can add it to your personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```groovy
params.project = 'project_ID' // eg. b2017123
```

## 3.5) Configuration: Other clusters
It is entirely possible to run this pipeline on other clusters, though you will need to set up your own config file so that the script knows where to find your reference files and how your cluster works.

> If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `-profile <name>` when running the pipeline.

If you are the only person to be running this pipeline, you can create your config file as `~/.nextflow/config` and it will be applied every time you run Nextflow. Alternatively, save the file anywhere and reference it when running the pipeline with `-c path/to/config`.

To specify your cluster environment, add the following line to your config file

```groovy
process {
  executor = 'YOUR_SYSTEM_TYPE'
}
```

This can be also be done on the command line:

```
-e.process={executor='SLURM'}
```

Many different cluster types are supported by Nextflow. For more information, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html).

Note that you may need to specify cluster options, such as a project or queue. To do so, use the `clusterOptions` config option.

```groovy
process {
  executor = 'SLURM'
  clusterOptions = '-A myproject'
}
```

### Reference Genomes
The nf-core/methylseq pipeline needs a reference genome for read alignment. Support for many common genomes is built in if running on UPPMAX or AWS, by using [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

If you don't want to use the illumina iGenomes you can supply either a Bismark reference or a FASTA file. If a Bismark reference is specified, the pipeline won't have to generate it and will be finished quite a bit faster. If a FASTA file is supplied then the Bismark reference will be built when the pipeline starts. Use the command line option `--saveReference` to keep the generated references so that they can be added to your config and used again in the future. Use  `--bismark_index` or `--fasta` to specify the paths to the reference.

Alternatively, you can add the paths to your NextFlow config under a relevant id and just specify this id with `--genome ID` when you run the pipeline:

```groovy
params {
  genomes {
    'YOUR-ID' {
      bismark  = '<PATH TO BISMARK REF>/BismarkIndex'
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
