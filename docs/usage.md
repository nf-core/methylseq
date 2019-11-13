# nf-core/methylseq: Usage

## Table of contents

<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:false -->
* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
  * [Bismark and bwa-meth workflow](#bismark-and-bwa-meth-workflow)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
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
  * [`--skip_preseq`](#--skip_preseq)
  * [`--non_directional`](#--non_directional)
  * [`--comprehensive`](#--comprehensive)
  * [`--relax_mismatches` and `--num_mismatches`](#--relax_mismatches-and---num_mismatches)
  * [`--unmapped`](#--unmapped)
  * [`--save_trimmed`](#--save_trimmed)
  * [`--save_align_intermeds`](#--save_align_intermeds)
  * [`--min_depth`](#--min_depth)
  * [`--meth_cutoff`](#--meth_cutoff)
  * [`--ignore_flags`](#--ignore_flags)
  * [`--methyl_kit`](#--methyl_kit)
  * [`--known_splices`](#--known_splices)
  * [`--slamseq`](#--slamseq)
  * [`--local_alignment`](#--local_alignment)
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
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
  * [`--multiqc_config`](#--multiqc_config)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--project`](#--project)
  * [`--clusterOptions`](#--clusteroptions)
<!-- TOC END -->

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

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
nextflow run nf-core/methylseq --reads '*_R{1,2}.fastq.gz' -profile docker
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

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/methylseq`](http://hub.docker.com/r/nfcore/methylseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/methylseq`](http://hub.docker.com/r/nfcore/methylseq/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--single_end`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--single_end --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

## Reference genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta   = '<path to the genome fasta file>' // Used if no index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### Supplying reference indices

If you don't want to use the Illumina iGenomes references, you can supply your own reference genome.

The minimum requirement is just a FASTA file - the pipeline will automatically generate the relevant reference index from this. You can use the command line option `--save_reference` to keep the generated references so that they can be added to your config and used again in the future. The bwa-meth workflow always needs a FASTA file, for methylation calling.

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
  * Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been
* `--three_prime_clip_r2 <NUMBER>`
  * Instructs Trim Galore to re move bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.

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

### `--skip_preseq`

Use `--skip_preseq` to skip the QC step where the pipeline runs Preseq. This often fails with low complexity libraries.

### `--pbat`

Using the `--pbat` parameter will affect the trimming (see above) and also set the `--pbat` flag when aligning with Bismark. It tells Bismark to align complementary strands (the opposite of `--directional`).

### `--non_directional`

By default, Bismark assumes that libraries are directional and does not align against complementary strands. If your library prep was not directional, use `--non_directional` to align against all four possible strands.

Note that the `--single_cell` and `--zymo` parameters both set the `--non_directional` workflow flag automatically.

### `--comprehensive`

By default, the pipeline only produces data for cytosine methylation states in CpG context. Specifying `--comprehensive` makes the pipeline give results for all cytosine contexts. Note that for large genomes (e.g. Human), these can be massive files. This is only recommended for small genomes (especially those that don't exhibit strong CpG context methylation specificity).

If specified, this flag instructs the Bismark methylation extractor to use the `--comprehensive` and `--merge_non_CpG` flags. This produces coverage files with information from about all strands and cytosine contexts merged into two files - one for CpG context and one for non-CpG context.

If using the bwa-meth workflow, the flag makes MethylDackel report CHG and CHH contexts as well.

### `--relax_mismatches` and `--num_mismatches`

By default, Bismark is pretty strict about which alignments it accepts as valid. If you have good reason to believe that your reads will contain more mismatches than normal, these flags can be used to relax the stringency that Bismark uses when accepting alignments. This can greatly improve the number of aligned reads you get back, but may negatively impact the quality of your data.

`--num_mismatches` is `0.2` by default in Bismark, or `0.6` if `--relax_mismatches` is specified. `0.6` will allow a penalty of `bp * -0.6` - for 100bp reads, this is `-60`. Mismatches cost `-6`, gap opening `-5` and gap extension `-2`. So, `-60` would allow 10 mismatches or ~ 8 x 1-2bp indels.

### `--unmapped`

Use the `--unmapped` flag to set the `--unmapped` flag with Bismark align and save the unmapped reads to FastQ files.

### `--save_trimmed`

By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete.

### `--save_align_intermeds`

By default intermediate BAM files will not be saved. The final BAM files created after the deduplication step are always. Set to true to also copy out BAM files from the initial Bismark alignment step. If `--skip_deduplication` or `--rrbs` is specified then BAMs from the initial alignment will always be saved.

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

### `--slamseq`

Specify to run Bismark with the `--slam` flag to run bismark in [SLAM-seq mode](https://github.com/FelixKrueger/Bismark/blob/master/CHANGELOG.md#slam-seq-mode) (only works with `--aligner bismark_hisat`)

### `--local_alignment`

Specify to run Bismark with the `--local` flag to allow soft-clipping of reads. This should only be used with care in certain single-cell applications or PBAT libraries, which may produce chimeric read pairs. (See [Wu et al.](https://doi.org/10.1093/bioinformatics/btz125) (doesn't work with `--aligner bwameth`)

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack/).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

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
