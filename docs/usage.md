# nf-core/methylseq Usage

## Bismark and bwa-meth workflow

The nf-core/methylseq package is actually two pipelines in one. The default workflow uses [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/): unless specified otherwise, nf-core/methylseq will run this pipeline.

The second workflow uses [BWA-Meth](https://github.com/brentp/bwa-meth) and [MethylDackyl](https://github.com/dpryan79/methyldackel) instead of Bismark. To run this workflow, run the pipeline with the command line flag `--aligner bwameth`.

## Running the pipeline
The main command for running the pipeline is as follows:

```bash
nextflow run nf-core/methylseq [parameters]
```

Note that the pipeline will create files in your working directory:

```bash
work/           # Directory containing the nextflow working files
results/        # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
.nextflow/      # Nextflow cache and history information
```


## Inputs and outputs

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

### `--outdir`
The output directory where the results will be saved.

## Software
There are different ways to provide the required software dependencies for the pipeline. The recommended method is to use either Docker or Singularity. The command line flags below are only needed if you are using the `standard` config profile (see below). It is not required with the the other profiles.

If using Bioconda or manual installation, no command line option is required.

### `-with-singularity`
Flag to enable use of singularity. The image will automatically be pulled from the internet. If running offline, follow the option with the path to the image file.

### `-with-docker`
Flag to enable docker. The image will automatically be pulled from Dockerhub.

## Reference Genomes

### `--genome`
To make it easier to run an analysis with a common reference genome, the pipeline comes bundled with configuration paths for the AWS-iGenomes resource. If using these references, you can specify just the genome ID using `--genome`.

The default AWS-iGenome paths are for files hosted on AWS s3. If running regularly, it is recommended that you download your own copy of the reference genomes and set the path root using `--igenomes_base` (or `params.igenomes_base` in a config file).

If using the `uppmax` config profile (see below), the iGenomes base is already set to a local copy of the files held on the UPPMAX clusters.

See [`conf/igenomes.config`](conf/igenomes.config) for a list of all of the supported reference genomes and their keys. Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

### Supplying reference indices

If you don't want to use the illumina iGenomes references, you can supply your own reference genome.

The minimum requirement is just a FASTA file - the pipeline will automatically generate the relevant reference index from this. You can use the command line option `--saveReference` to keep the generated references so that they can be added to your config and used again in the future. The bwa-meth workflow always needs a FASTA file, for methylation calling.

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

### `--saveReference`
Supply this parameter to save any generated reference genome files to your results folder. These can then be used for future pipeline runs, reducing processing times.

## Adapter Trimming
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

### `--notrim`
Specifying `--notrim` will skip the adapter trimming step. Use this if your input FastQ files have already been trimmed outside of the workflow.

### `--nodedup`
By default, the pipeline includes a deduplication step after alignment. Use `--nodedup` on the command line to skip this step. This is automatically set if using `--rrbs` for the workflow.

## Preset configurations

### `-profile`
Use this parameter to choose a configuration profile. See the [installation documentation](installation.md#33-configuration-profiles) for more information about profiles.

Profiles available with the pipeline are:

* `standard`
    * The default profile, used if `-profile` is not specified.
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
* `aws`
    * A starter configuration for running the pipeline on Amazon Web Services.
    * Specifies docker configuration and uses the `spark` job executor
    * Requires additional configuration to run - see the documentation dedicated to this topic.
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile.

## Additional parameters

### `--pbat`
Using the `--pbat` parameter will affect the trimming (see above) and also set the `--pbat` flag when aligning with Bismark. It tells Bismark to align complementary strands (the opposite of `--directional`).

### `--non_directional`
By default, Bismark assumes that libraries are directional and does not align against complementary strands. If your library prep was not direcional, use `--non_directional` to align against all four possible strands.

Note that the `--single_cell` and `--zymo` parameters both set the `--non_directional` workflow flag automatically.

### `--comprehensive`
By default, the pipeline only produces data for cytosine methylation states in CpG context. Specifying `--comprehensive` makes the pipeline give results for all cytosine contexts. Note that for large genomes (eg. Human), these can be massive files. This is only recommended for small genomes (especially those that don't exhibit strong CpG context methylation specificity).

If specified, this flag instructs the Bismark methylation extractor to use the `--comprehensive` and `--merge_non_CpG` flags. This produces coverage files with information from about all strands and cytosine contexts merged into two files - one for CpG context and one for non-CpG context.

If using the bwa-meth workflow, the flag makes MethylDackel report CHG and CHH contexts as well.

### `--relaxMismatches` and `--numMismatches`

By default, Bismark is pretty strict about which alignments it accepts as valid. If you have good reason to believe that your reads will contain more mismatches than normal, these flags can be used to relax the stringency that Bismark uses when accepting alignments. This can greatly improve the number of aligned reads you get back, but may negatively impact the quality of your data.

`--numMismatches` is `0.2` by default in Bismark, or `0.6` if `--relaxMismatches` is specified. `0.6` will allow a penalty of `bp * -0.6` - for 100bp reads, this is `-60`. Mismatches cost `-6`, gap opening `-5` and gap extension `-2`. So, `-60` would allow 10 mismatches or ~ 8 x 1-2bp indels.

### `--unmapped`
Use the `--unmapped` flag to set the `--unmapped` flag with Bismark align and save the unmapped reads to FastQ files.

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete.

### `--saveAlignedIntermediates`
By default intermediate BAM files will not be saved. The final BAM files created after the deduplication step are always. Set to true to also copy out BAM files from the initial Bismark alignment step. If `--nodedup` or `--rrbs` is specified then BAMs from the initial alignment will always be saved.

### `--mindepth`
Specify to specify a minimum read coverage for MethylDackel to report a methylation call.

### `--ignoreFlags`
Specify to run MethylDackel with the `--ignoreFlags` flag to ignore SAM flags.


## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits on UPPMAX with an error code of `143` or `137` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Maximum resource requests
All resource requests are checked against the following default limits (`standard` config profile shown):

```bash
--max_memory '128.GB'
--max_cpus '16'
--max_time '240.h'
```

If a task requests more than this amount, it will be reduced to this threshold.

To adjust these limits, specify them on the command line, eg. `--max_memory '64.GB'`.

Note that these limits are the maximum to be used _per task_. Nextflow will automatically attempt to parallelise as many jobs as possible given the available resources.

## Other command line parameters

### `-name` _(single dash)_
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic. This name is also used in MultiQC reports if specified.

### `-resume` _(single dash)_
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `-c` _(single dash)_
Specify the path to a specific config file (this is a core NextFlow command). You can use this in combination with configuration profiles to override defaults.

### `--multiqc_config`
If you would like to supply a custom config file to MultiQC, you can specify a path with `--multiqc_config`. This is used _instead of_ the [config file](../conf/multiqc_config.yaml) that comes with the pipeline.

### `--project`
UPPMAX profile only: Cluster project for SLURM job submissions.

### `--clusterOptions`
UPPMAX profile only: Submit arbitrary SLURM options.
