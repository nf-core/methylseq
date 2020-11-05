# nf-core/methylseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/methylseq/usage](https://nf-co.re/methylseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
  * [Bismark and bwa-meth workflow](#bismark-and-bwa-meth-workflow)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

### Bismark and bwa-meth workflow

The nf-core/methylseq package is actually two pipelines in one. The default workflow uses [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as alignment tool: unless specified otherwise, nf-core/methylseq will run this pipeline.

Since bismark v0.21.0 it is also possible to use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) as alignment tool. To run this workflow, invoke the pipeline with the command line flag `--aligner bismark_hisat`. HISAT2 also supports splice-aware alignment if analysis of RNA is desired (e.g. [SLAMseq](https://science.sciencemag.org/content/360/6390/800) experiments), a file containing a list of known splicesites can be provided with `--known_splices`.

The second workflow uses [BWA-Meth](https://github.com/brentp/bwa-meth) and [MethylDackel](https://github.com/dpryan79/methyldackel) instead of Bismark. To run this workflow, run the pipeline with the command line flag `--aligner bwameth`.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/methylseq --input '*_R{1,2}.fastq.gz' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/methylseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/methylseq releases page](https://github.com/nf-core/methylseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

### `--unmapped`

Use the `--unmapped` flag to set the `--unmapped` flag with Bismark align and save the unmapped reads to FastQ files.

### `--save_align_intermeds`

### `--min_depth`

Specify to specify a minimum read coverage for MethylDackel to report a methylation call.

### `--meth_cutoff`

Use this to specify a minimum read coverage to report a methylation call during Bismark's `bismark_methylation_extractor` step.

### `--ignore_flags`

Specify to run MethylDackel with the `--ignore_flags` flag to ignore SAM flags.

### `--methyl_kit`

Specify to run MethylDackel with the `--methyl_kit` flag to produce files suitable for use with the methylKit R package.

### `--known_splices`

Specify to run Bismark with the `--known-splicesite-infile` flag to run splice-aware alignment using HISAT2. A `.gtf` file has to be provided from which a list of known splicesites is created by the pipeline. (only works with `--aligner bismark_hisat`)

### `--local_alignment`

Specify to run Bismark with the `--local` flag to allow soft-clipping of reads. This should only be used with care in certain single-cell applications or PBAT libraries, which may produce chimeric read pairs. (See [Wu et al.](https://doi.org/10.1093/bioinformatics/btz125) (doesn't work with `--aligner bwameth`)

### `--bismark_align_cpu_per_multicore`

The pipeline makes use of the `--multicore` option for Bismark align. When using this option,
Bismark uses a large number of CPUs for every `--multicore` specified. The pipeline
calculates the number of `--multicore` based on the resources available to the task.
It divides the available CPUs by 3, or by 5 if any of `--single_cell`, `--zymo` or `--non_directional`
are specified. This is based on usage for a typical mouse genome.

You may find when running the pipeline that Bismark is not using this many CPUs. To fine tune the
usage and speed, you can specify an integer with `--bismark_align_cpu_per_multicore` and the pipeline
will divide the available CPUs by this value instead.

See the [bismark documentation](https://github.com/FelixKrueger/Bismark/tree/master/Docs#alignment)
for more information.

### `--bismark_align_mem_per_multicore`

Exactly as above, but for memory. By default, the pipeline divides the available memory by `13.GB`,
or `18.GB` if any of `--single_cell`, `--zymo` or `--non_directional` are specified.

Note that the final `--multicore` value is based on the lowest limiting factor of both CPUs and memory.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/methylseq`](https://hub.docker.com/r/nfcore/methylseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/methylseq`](https://hub.docker.com/r/nfcore/methylseq/)
* `podman`
  * A generic configuration profile to be used with [Podman](https://podman.io/)
  * Pulls software from Docker Hub: [`nfcore/methylseq`](https://hub.docker.com/r/nfcore/methylseq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity or Podman.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `bismark_align` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: bismark_align {
    memory = 32.GB
  }
}
```

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--multiqc_config`

If you would like to supply a custom config file to MultiQC, you can specify a path with `--multiqc_config`. This is used _instead of_ the [config file](../conf/multiqc_config.yaml) that comes with the pipeline.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--project`

UPPMAX profile only: Cluster project for SLURM job submissions.

### `--clusterOptions`

UPPMAX profile only: Submit arbitrary SLURM options.
