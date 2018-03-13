# NGI-MethylSeq

## v0.4dev
* Renaming and moving [SciLifeLab/NGI-MethylSeq](https://github.com/SciLifeLab/NGI-MethylSeq/) to [nf-core/methylseq](https://github.com/nf-core/methylseq/)
* Refactoring multi-core parameters for Bismark alignment and methylation extraction
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
