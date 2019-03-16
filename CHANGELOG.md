# nf-core/methylseq

## v1.4dev

#### New features
* Added support for centralized configuration profiles [nf-core/configs](https://github.com/nf-core/configs)
* Add `meth_cutoff = 5` default for Bismark
  * Now, 5 reads (previously 1) are required to call a methylation site with `bismark_methylation_extractor`.
  * New `--meth_cutoff` option makes this configurable.
  * Use `--meth_cutoff 1` on the command line or `params.meth_cutoff = 1` in a config file to revert to previous behaviour.

#### Bug fixes
* Fixed a bug that caused conda dependencies to be resolved very slowly

## [v1.3](https://github.com/nf-core/methylseq/releases/tag/1.3) - 2019-02-01

#### New features
* Added [preseq](http://smithlabresearch.org/software/preseq/) analysis to calculate sample complexity.
  * This new step can help decide sufficient sequencing depth has been reached.

#### Bug fixes
* Fixed new bug that meant pipeline only worked with one sample at a time [#66](https://github.com/nf-core/methylseq/issues/66)
  * Introduced in previous release. TrimGalore onwards would only process one sample.


## [v1.2](https://github.com/nf-core/methylseq/releases/tag/1.2) - 2019-01-02

#### New features
* Trim 9bp from both ends of both reads for PBAT mode.
* Save `where_are_my_files.txt` to results directory to inform the user about missing intermediate files [#42](https://github.com/nf-core/methylseq/issues/42)

#### Software updates
* Fastqc `0.11.7` > `0.11.8`
* Bowtie2 `2.3.4.2` > `2.3.4.3`
* Bismark `0.19.1` > `0.20.0`
* Qualimap `2.2.2a` > `2.2.2b`
* Picard `2.18.11` > `2.18.21`
* MultiQC `1.6` > `1.7`

#### Bug fixes
* Fixed error when running the pipeline with `--unmapped`
  * Previously, could result in error `Error ~ No such variable: bismark_unmapped`
* Fixed error where single-sample reports could mix up log files [#48](https://github.com/nf-core/methylseq/issues/48)
* Fixed bug in MultiQC process that skipped results from some tools
* Supply available memory as argument to Picard MarkDuplicates


## [v1.1](https://github.com/nf-core/methylseq/releases/tag/1.1) - 2018-08-09

* Tests simplified - now work by simply using the `test` config profile
    * eg: `nextflow run nf-core/methylseq -profile test,docker`
    * Removed previous `run_test.sh` script and data
* New `Singularity` build script for direct compatibility with [singularity-hub](https://singularity-hub.org/)
* Minor improvements to the docs
* A number of boilerplate nf-core code updates
* Updated `process$name` nextflow syntax to avoid warnings in new versions of nextflow
* Updated software tools
    * `trim-galore` `v0.4.5` update to `0.5.0`
    * `samtools` `v1.8` update to `1.9`
    * `bowtie2` `v2.3.4.1` update to `2.3.4.2`
    * `multiqc` `v1.5` update to `1.6`
    * `picard` `v2.18.2` update to `2.18.11`
    * `bwameth` `v0.2.0` update to `0.2.2`



## [v1.0](https://github.com/nf-core/methylseq/releases/tag/1.0) - 2018-04-17
Version 1.0 marks the first release of this pipeline under the nf-core flag. It also marks a significant step up in the maturity of the workflow, with everything now in a single script and both aligner workflows fully supported.

* Renamed and moved [SciLifeLab/NGI-MethylSeq](https://github.com/SciLifeLab/NGI-MethylSeq/) to [nf-core/methylseq](https://github.com/nf-core/methylseq/)
* Merged bwa-meth and bismark pipeline scripts, now chosen with `--aligner` flag
* Refactored multi-core parameters for Bismark alignment and methylation extraction
* Rewrote most of the documentation
* Changed the Docker container to use Bioconda installations

---

Previous to these releases, this pipeline was called [SciLifeLab/NGI-MethylSeq](https://github.com/SciLifeLab/NGI-MethylSeq):

## v0.4dev
* Fixed MultiQC channel bug
* Integrated config for QBiC Tuebingen
* Numerous small container bugfixes
* Refactored how the config is loaded
* Fix for resource limit function, improved resource request defaults
* Fix for iGenomes base path in configs

## [v0.3.1](https://github.com/SciLifeLab/NGI-MethylSeq/releases/tag/0.3.1) - 2017-09-05
* Include base profile name and documentation about Singularity.
* Testing automated docker hub image tagging for releases.


## [v0.3](https://github.com/SciLifeLab/NGI-MethylSeq/releases/tag/0.3) - 2017-09-01
* Fix `--rrbs` mode ([#24](https://github.com/SciLifeLab/NGI-MethylSeq/issues/24))
* Fixed fairly major bug where only a single sample would run past alignment
* Merged test scripts and rewrote to use command line flags / new travis script.
* Refactored software version collection code to be more resilient and cleaner / easier to maintain.
* Dropped support for environment modules and added support for use of Singularity on UPPMAX

## [v0.2](https://github.com/SciLifeLab/NGI-MethylSeq/releases/tag/0.2) - 2017-07-17
First (semi-) stable release of the new NGI-MethylSeq pipeline, as we head towards deployment in production.
