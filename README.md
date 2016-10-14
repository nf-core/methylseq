# NGI-MethylSeq
Methylation (BS-Seq) Best Practice pipeline for bisulfite sequencing data analysis at
the [National Genomics Infastructure](https://ngisweden.scilifelab.se/)
at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

# UNDER DEVELOPMENT
This pipeline is currently under active development and has very little in the way of testing.
There's a pretty good chance that it won't run. Any help debugging is very welcome! Please
fork, make changes and submit a pull request.

## Installation
### NextFlow installation
To use this pipeline, you need to have a working version of NextFlow installed. You can find more
information about this pipeline tool at [nextflow.io](http://www.nextflow.io/). The typical installation
of NextFlow looks like this:

```
curl -fsSL get.nextflow.io | bash
mv ./nextflow ~/bin
```
Note that if you're running on a Swedish UPPMAX cluster you can load NextFlow as an environment module:
```
module load nextflow
```

### NextFlow configuration
Next, you need to set up a config file so that NextFlow knows how to run and where to find reference
indexes. You can find an example configuration file for UPPMAX (milou) with this repository:
[`example_uppmax_config`](https://github.com/SciLifeLab/NGI-MethylSeq/blob/master/example_uppmax_config).

Copy this file to `~/.nextflow/config` and edit the line `'-A YOUR_PROJECT_ID'` to contain your
UPPMAX project identifier.

It is entirely possible to run this pipeline on other clusters - just note that you may need to customise
the `process` environment (eg. if you're using a cluster system other than SLURM) and the paths to reference
files.

### Pipeline installation
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub when run if
`SciLifeLab/NGI-MethylSeq` is specified as the pipeline name.

If you prefer, you can download the files yourself from GitHub and run them directly:
```
git clone https://github.com/SciLifeLab/NGI-MethylSeq.git
nextflow NGI-MethylSeq/main.nf
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```
nextflow SciLifeLab/NGI-MethylSeq --reads '*_R{1,2}.fastq.gz' --genome GRCm38
```

### `--reads`
Location of the input FastQ files:
```
 --reads 'path/to/data/sample_*_{1,2}.fastq'
```

Note that the `{1,2}` parentheses are required to specify paired end data. Running `--reads '*.fastq'` will treat
all files as single end. Also, note that the file path should be in quotation marks to prevent shell glob expansion.

If left unspecified, the pipeline will assume that the data is in a directory called `data` in the working directory.

### `--genome`
The reference genome to use of the analysis, needs to be one of the genome specified in the config file.
The human `GRCh37` genome is set as default.
```
--genome 'GRCm38'
```
The `example_uppmax_config` file currently has the location of references for most of the
[Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)
held on UPPMAX.

### Trimming Parameters
The pipeline accepts a number of parameters to change how the trimming is done, according to your data type.
The following settings are set for these command line flags:

| Parameter       | 5' R1 Trim | 5' R2 Trim | 3' R1 Trim | 3' R2 Trim |
|-----------------|------------|------------|------------|------------|
| `--pbat`        | 6          | 6          | 0          | 0          |
| `--single_cell` | 9          | 9          | 0          | 0          |
| `--epignome`    | 6          | 6          | 6          | 6          |
| `--accel`       | 10         | 15         | 10         | 10         |
| `--cegx`        | 6          | 6          | 2          | 2          |

You can specify custom trimming parameters as follows:

* `--clip_r1 <NUMBER>`
* `--clip_r2 <NUMBER>`
* `--three_prime_clip_r1 <NUMBER>`
* `--three_prime_clip_r2 <NUMBER>`

Finally, specifying `--rrbs` will pass on the `--rrbs` parameter to TrimGalore!

## Bismark Parameters
Using the `--pbat` parameter will affect the trimming (see above) and also set the `--pbat` flag when
aligning with Bismark.

Using the `--single_cell` parameter will set the `--non_directional` flag when aligning with Bismark.
This can also be set with `--non_directional` (doesn't affect trimming).

Use the `--unmapped` flag to set the `--unmapped` flag with Bismark align and save the unmapped reads.

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes.

## Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
Written by Phil Ewels (@ewels) and Rickard Hammarén (@Hammarn).

<p align="center"><a href="stand_alone/http://www.scilifelab.se/" target="_blank"><img src="Docs/images/SciLifeLab_logo.png" title="SciLifeLab"></a></p>
