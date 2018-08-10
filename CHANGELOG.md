# nf-core/methylseq

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
