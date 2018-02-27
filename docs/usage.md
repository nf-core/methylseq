# nf-core/MethylSeq Usage

## Bismark and bwa-meth workflow

The nf-core/MethylSeq package is actually two pipelines in one. The default (and most polished) workflow uses [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) for processing and can be found in `bismark.nf`. Unless specified otherwise, nf-core/MethylSeq will run this pipeline.

The second included workflow uses [BWA-Meth](https://github.com/brentp/bwa-meth) and [MethylDackyl](https://github.com/dpryan79/methyldackel) instead of Bismark and can be found in `bwa-meth.nf`. To run this workflow, you must download the repository files and explicitly call `nextflow run /path/to/files/bwa-meth.nf`. Note that this workflow has not been tested to the same extent and may have bugs.

All of the documentation refers to the Bismark workflow at this stage, though most of it also applies to the bwa-meth workflow.

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/MethylSeq --genome GRCh37 --reads '*_R{1,2}.fastq.gz' -profile docker
```

This will launch the pipeline with the `docker` configuration profile (Swedish UPPMAX users use `-profile uppmax`). See below for more information about profiles.

Note that the pipeline will create files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```


### `-profile`
Use this parameter to choose a configuration profile. `-profile docker` is likely useful for users outside of Sweden. Available profiles are:

* `standard`
    * The default profile - this is used if `-profile` is not specified at all
    * Uses sensible defaults for requirements, runs using the `local` executor (native system calls) and expects all software to be installed and available on the `PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles below.
* `uppmax`
    * Designed to be used on the Swedish [UPPMAX](http://uppmax.uu.se/) clusters such as `milou`, `rackham`, `bianca` and `irma`
    * Launches jobs using the SLURM executor.
    * Uses [Singularity](http://singularity.lbl.gov/) to provide all software requirements
    * Comes with locations for illumina iGenome reference paths built in
    * Use with `--project` to provide your UPPMAX project ID.
* `uppmax_devel`
    * Uses the milou [devel partition](http://www.uppmax.uu.se/support/user-guides/slurm-user-guide/#tocjump_030509106905141747_8) for testing the pipeline quickly.
    * Not suitable for proper analysis runs
* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Runs using the `local` executor and pulls software from dockerhub: [`scilifelab/ngi-rnaseq`](http://hub.docker.com/r/scilifelab/ngi-rnaseq/)
* `aws`
    * A starter configuration for running the pipeline on Amazon Web Services.
    * Specifies docker configuration and uses the `spark` job executor
    * Requires additional configuration to run - see the documentation dedicated to this topic.
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile.

## Input Data

### `--reads`
Location of the input FastQ files:

```bash
 --reads 'path/to/data/sample_*_{1,2}.fastq'
```

**NB: Must be enclosed in quotes!**

Note that the `{1,2}` parentheses are required to specify paired end data. The file path should be in quotation marks to prevent shell glob expansion.

If left unspecified, the pipeline will assume that the data is in a directory called `data` in the working directory (`data/*{1,2}.fastq.gz`).

### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example: `--singleEnd --reads '*.fastq'`

It is not possible to run a mixture of single-end and paired-end files in one run.

## Reference Genomes

### `--genome`
The reference genome to use for the analysis, needs to be one of the genome specified in the config file. This is `False` by default and needs to be specified (unless index files are supplied, see below).

See [`conf/uppmax.config`](conf/uppmax.config) for a list of the supported reference genomes and their keys. Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

If you're not running on UPPMAX (the default profile), you can create your own config file with paths to your reference genomes. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to add this.

The syntax for this reference configuration is as follows:

```groovy
params {
  genomes {
    'GRCh37' {
      bismark = '<path to the bismark index folder>'
      fasta   = '<path to the genome fasta file>' // Used if no bismark index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### `--bismark_index`, `--fasta`
If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--bismark_index '[path to Bismark index]'
--fasta '[path to Fasta reference]'
```

### `--saveReference`
Supply this parameter to save any generated reference genome files to your results folder. These can then be used for future pipeline runs, reducing processing times.

## Adapter Trimming
The pipeline accepts a number of parameters to change how the trimming is done, according to your data type.
The following settings are set for these command line flags:

| Parameter       | 5' R1 Trim | 5' R2 Trim | 3' R1 Trim | 3' R2 Trim |
|-----------------|------------|------------|------------|------------|
| `--pbat`        | 6          | 9          | 6          | 9          |
| `--single_cell` | 6          | 6          | 6          | 6          |
| `--epignome`    | 8          | 8          | 8          | 8          |
| `--accel`       | 10         | 15         | 10         | 10         |
| `--zymo`        | 10         | 15         | 10         | 10         |
| `--cegx`        | 6          | 6          | 2          | 2          |

You can specify custom trimming parameters as follows:

* `--clip_r1 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).
* `--clip_r2 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).
* `--three_prime_clip_r1 <NUMBER>`
  * Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been
* `--three_prime_clip_r2 <NUMBER>`
  * Instructs Trim Galore to re move bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.

### `--rrbs`
Specifying `--rrbs` will pass on the `--rrbs` parameter to TrimGalore! See the [TrimGalore! documentation](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md#rrbs-specific-options-mspi-digested-material) to read more about the effects of this option.

### `--notrim`
Specifying `--notrim` will skip the adapter trimming step. Use this if your input FastQ files have already been trimmed outside of the workflow.

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete.


## Bismark Parameters

### `--pbat`
Using the `--pbat` parameter will affect the trimming (see above) and also set the `--pbat` flag when aligning with Bismark. It tells Bismark to align complementary strands (the opposite of `--directional`).

### `--non_directional`
By default, Bismark assumes that libraries are directional and does not align against complementary strands. If your library prep was not direcional, use `--non_directional` to align against all four possible strands.

Note that the `--single_cell` and `--zymo` parameters both set the `--non_directional` workflow flag automatically.

### `--comprehensive`
This flag instructs the Bismark methylation extractor to use the `--comprehensive` and `--merge_non_CpG` flags. This produces coverage files with information from about all strands and cytosine contexts merged into two files - one for CpG context and one for non-CpG context. Note that for large genomes (eg. Human), these can be massive files. This is only recommended for small genomes (especially those that don't exhibit strong CpG context methylation specificity).

### `--relaxMismatches` and `--numMismatches`

By default, Bismark is pretty strict about which alignments it accepts as valid. If you have good reason to believe that your reads will contain more mismatches than normal, these flags can be used to relax the stringency that Bismark uses when accepting alignments. This can greatly improve the number of aligned reads you get back, but may negatively impact the quality of your data.

`--numMismatches` is `0.2` by default in Bismark, or `0.6` if `--relaxMismatches` is specified. `0.6` will allow a penalty of `bp * -0.6` - for 100bp reads, this is `-60`. Mismatches cost `-6`, gap opening `-5` and gap extension `-2`. So, `-60` would allow 10 mismatches or ~ 8 x 1-2bp indels.

### `--unmapped`
Use the `--unmapped` flag to set the `--unmapped` flag with Bismark align and save the unmapped reads to FastQ files.

### `--nodedup`
By default, the pipeline includes a deduplication step after alignment. Use `--nodedup` on the command line to skip this step. This is automatically set if using `--rrbs` for the workflow.

### `--saveAlignedIntermediates`
By default intermediate BAM files will not be saved. The final BAM files created after the Bismark deduplication step are always saved. Set to true to also copy out BAM files from the initial Bismark alignment step. If `--nodedup` or `--rrbs` is specified then BAMs from the initial alignment will always be saved.

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits on UPPMAX with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

## Other command line parameters
### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes.

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults.

### `--clusterOptions`
Submit arbitrary SLURM options (UPPMAX profile only). For instance, you could use `--clusterOptions '-p devcore'`
to run on the development node (though won't work with default process time requests).

### `--multiqc_config`
If you would like to supply a custom config file to MultiQC, you can specify a path with `--multiqc_config`. This is used instead of the config file specific to the pipeline.

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
