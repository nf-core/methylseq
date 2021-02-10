

# nf-core/methylseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/methylseq/usage](https://nf-co.re/methylseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
  * [Bismark, bwa-meth and biscuit workflow](#bismark-and-bwa-meth-workflow)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
<<<<<<< HEAD
* [Main arguments](#main-arguments)
  * [`-profile`](#-profile)
  * [`--reads`](#--reads)
  * [`--single_end`](#--single_end)
* [Reference genomes](#reference-genomes)
  * [`--genome` (using iGenomes)](#--genome-using-igenomes)
  * [`--fasta`](#--fasta)
  * [`--igenomes_ignore`](#--igenomes_ignore)
  * [Supplying reference indices](#supplying-reference-indices)
  * [`--save_reference`](#--save_reference)
* [Additional parameters](#additional-parameters)
  * [Adapter Trimming](#adapter-trimming)
    * [`--rrbs`](#--rrbs)
    * [`--pbat`](#--pbat)
    * [`--skip_trimming`](#--skip_trimming)
  * [`--skip_deduplication`](#--skip_deduplication)
  * [`--non_directional`](#--non_directional)
  * [`--comprehensive`](#--comprehensive)
  * [`--cytosine_report`](#--cytosine_report)
  * [`--relax_mismatches` and `--num_mismatches`](##--relax_mismatches-and---num_mismatches)
  * [`--unmapped`](#--unmapped)
  * [`--save_trimmed`](#--save_trimmed)
  * [`--save_align_intermeds`](#--save_align_intermeds)
  * [`--save_pileup_file`](#--save_pileup_file)
  * [`--epiread`](#--epiread)
  * [`--assets_dir`](#--assets_dir)
  * [`--min_depth`](#--min_depth)
  * [`--meth_cutoff`](#--meth_cutoff)
  * [`--ignore_flags`](#--ignore_flags)
  * [`--methyl_kit`](#--methyl_kit)
  * [`--known_splices`](#--known_splices)
  * [`--slamseq`](#--slamseq)
  * [`--local_alignment`](#--local_alignment)
  * [`--bismark_align_cpu_per_multicore`](#--bismark_align_cpu_per_multicore)
  * [`--bismark_align_mem_per_multicore`](#--bismark_align_mem_per_multicore)
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
  * [`--awscli`](#--awscli)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--multiqc_config`](#--multiqc_config)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

### Bismark, bwa-meth and biscuit workflow

The nf-core/methylseq package is actually three pipelines in one. The default workflow uses [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as alignment tool: unless specified otherwise, nf-core/methylseq will run this pipeline.
=======

## Introduction

The nf-core/methylseq package is actually two pipelines in one. The default workflow uses [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as alignment tool: unless specified otherwise, nf-core/methylseq will run this pipeline.
>>>>>>> 9218f1199bca434af49b54963eea91cfed572597

Since bismark v0.21.0 it is also possible to use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) as alignment tool. To run this workflow, invoke the pipeline with the command line flag `--aligner bismark_hisat`. HISAT2 also supports splice-aware alignment if analysis of RNA is desired (e.g. [SLAMseq](https://science.sciencemag.org/content/360/6390/800) experiments), a file containing a list of known splicesites can be provided with `--known_splices`.

The second workflow uses [BWA-Meth](https://github.com/brentp/bwa-meth) and [MethylDackel](https://github.com/dpryan79/methyldackel) instead of Bismark. To run this workflow, run the pipeline with the command line flag `--aligner bwameth`.

The third workflow uses [biscuit]([https://github.com/huishenlab/biscuit](https://github.com/huishenlab/biscuit)) . This workflow uses biscuit as an aligner, and biscuit-QC for quality control.  To run this workflow, run the pipeline with the command line flag `--aligner biscuit`

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

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: bismark_align {
    memory = 32.GB
  }
}
```

<<<<<<< HEAD
### Supplying reference indices

If you don't want to use the Illumina iGenomes references, you can supply your own reference genome.

The minimum requirement is just a FASTA file - the pipeline will automatically generate the relevant reference index from this. You can use the command line option `--save_reference` to keep the generated references so that they can be added to your config and used again in the future. The bwa-meth and biscuit workflows always need a FASTA file, for methylation calling. The FASTA is also required for the generation of Picard metrics.

### `--fasta`

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
# Single multifasta for genome
--fasta /path/to/genome.fa

# Bismark index directory
--bismark_index /path/to/ref/BismarkIndex/

# bwa-meth index filename base
# where for example the index files are called:
# /path/to/ref/genome.fa.bwameth.c2t.bwt
--bwa_meth_index /path/to/ref/genome.fa

# biscuit index filename base
# where for example the index files are called:
# /path/to/ref/genome.fa.bis.amb
--bwa_meth_index /path/to/ref/genome.fa

# Genome Fasta index file
--fasta_index /path/to/genome.fa.fai
```

### `--igenomes_ignore`

Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.

### `--save_reference`

Supply this parameter to save any generated reference genome files to your results folder. These can then be used for future pipeline runs, reducing processing times.

## Additional parameters

### Adapter Trimming

Bisulfite libraries often require additional base pairs to be removed from the ends of the reads before alignment. You can specify these custom trimming parameters as follows:

* `--clip_r1 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).
* `--clip_r2 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).
* `--three_prime_clip_r1 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been performed.
* `--three_prime_clip_r2 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.

The pipeline also accepts a number of presets for common bisulfite library preparation methods:

| Parameter       | 5' R1 Trim | 5' R2 Trim | 3' R1 Trim | 3' R2 Trim |
|-----------------|------------|------------|------------|------------|
| `--pbat`        | 6          | 9          | 6          | 9          |
| `--single_cell` | 6          | 6          | 6          | 6          |
| `--epignome`    | 8          | 8          | 8          | 8          |
| `--accel`       | 10         | 15         | 10         | 10         |
| `--zymo`        | 10         | 15         | 10         | 10         |
| `--cegx`        | 6          | 6          | 2          | 2          |

### `--rrbs`

Specifying `--rrbs` will pass on the `--rrbs` parameter to TrimGalore! See the [TrimGalore! documentation](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md#rrbs-specific-options-mspi-digested-material) to read more about the effects of this option.

This parameter also makes the pipeline skip the deduplication step.

### `--skip_trimming`

Specifying `--skip_trimming` will skip the adapter trimming step. Use this if your input FastQ files have already been trimmed outside of the workflow.

### `--skip_deduplication`

By default, the pipeline includes a deduplication step after alignment. Use `--skip_deduplication` on the command line to skip this step. This is automatically set if using `--rrbs` for the workflow.

### `--pbat`

Using the `--pbat` parameter will affect the trimming (see above) and also set the `--pbat` flag when aligning with Bismark and biscuit.
For bismark, it tells the aligner to align complementary strands (the opposite of `--directional`).
For biscuit, it tells the aligner to switch between reads in paired-end (or align to synthesized strand on single-end).

### `--non_directional`

By default, Bismark and biscuit assume that libraries are directional and do not align against complementary strands. If your library prep was not directional, use `--non_directional` to align against all four possible strands.

Note that the `--single_cell` and `--zymo` parameters both set the `--non_directional` workflow flag automatically.

### `--comprehensive`

By default, the pipeline only produces data for cytosine methylation states in CpG context. Specifying `--comprehensive` makes the pipeline give results for all cytosine contexts. Note that for large genomes (e.g. Human), these can be massive files. This is only recommended for small genomes (especially those that don't exhibit strong CpG context methylation specificity).

If specified, this flag instructs the Bismark methylation extractor to use the `--comprehensive` and `--merge_non_CpG` flags. This produces coverage files with information from about all strands and cytosine contexts merged into two files - one for CpG context and one for non-CpG context.

If using the bwa-meth workflow, the flag makes MethylDackel report CHG and CHH contexts as well.

if using the biscuit aligner, the flag generate the bedgraph file extracting all possible types from the pileup file (including c, cg, ch, hcg, gch).

### `--cytosine_report`

By default, Bismark does not produce stranded calls. With this option the output considers all Cs on both forward and reverse strands and reports their position, strand, trinucleotide context and methylation state.

### `--relax_mismatches` and `--num_mismatches`

By default, Bismark is pretty strict about which alignments it accepts as valid. If you have good reason to believe that your reads will contain more mismatches than normal, these flags can be used to relax the stringency that Bismark uses when accepting alignments. This can greatly improve the number of aligned reads you get back, but may negatively impact the quality of your data.

`--num_mismatches` is `0.2` by default in Bismark, or `0.6` if `--relax_mismatches` is specified. `0.6` will allow a penalty of `bp * -0.6` - for 100bp reads, this is `-60`. Mismatches cost `-6`, gap opening `-5` and gap extension `-2`. So, `-60` would allow 10 mismatches or ~ 8 x 1-2bp indels.

### `--unmapped`

Use the `--unmapped` flag to set the `--unmapped` flag with Bismark align and save the unmapped reads to FastQ files.

### `--save_trimmed`

By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete.

### `--save_align_intermeds`

By default intermediate BAM files will not be saved. The final BAM files created after the deduplication step are always. Set to true to also copy out BAM files from the initial Bismark alignment step. If `--skip_deduplication` or `--rrbs` is specified then BAMs from the initial alignment will always be saved.

### `--save_pileup_file`

When running with biscuit aligner, the methylation extraction is based on vcf file. By default these vcf files will not be saved. Set to true to also copy out the vcf-file and the index-vcf file.

### `--epiread`

[Epiread]([https://github.com/zhou-lab/biscuit/wiki/Convert-to-epiread-format](https://github.com/zhou-lab/biscuit/wiki/Convert-to-epiread-format)) format is a compact way of storing CpG retention pattern on the same read. This option will tell the biscuit workflow to generate epiread file. 

### `--assets_dir`

Path to a directory containing needed file for biscuit-QC step. The needed files for hg38,hg19 and mm10 can be found in [here](https://www.cse.huji.ac.il/~ekushele/assets.html). 
**This parameter is mandatory when running the pipeline using biscuit workflow**

### `--min_depth`

Specify to specify a minimum read coverage for MethylDackel or biscuit vcf2bed to report a methylation call.
=======
See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.
>>>>>>> 9218f1199bca434af49b54963eea91cfed572597

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
