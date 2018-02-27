# ![nf-core/MethylSeq](docs/images/MethylSeq_logo.png)

[![Build Status](https://travis-ci.org/nf-core/MethylSeq.svg?branch=master)](https://travis-ci.org/nf-core/MethylSeq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)

### Introduction

`nf-core/MethylSeq` is a bioinformatics best-practice analysis pipeline used for Methylation (BS-Seq) data analysis.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

### Choice of workflows

There are two separate workflows contained in this repository - one using [Bismark](https://github.com/FelixKrueger/Bismark) and one using [bwa-meth](https://github.com/brentp/bwa-meth) / [MethylDackel](https://github.com/dpryan79/methyldackel). The Bismark pipeline is being actively developed and maintained, the bwa-meth workflow is not _(currently)_. The Nextflow manifest specifies the Bismark pipeline as the default workflow, so the bwa-meth script will be ignored unless explicitly run.

### Documentation
The nf-core/MethylSeq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation and configuration](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)


### Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden. Written by Phil Ewels (@ewels) and Rickard Hammarén (@Hammarn).

---

[![SciLifeLab](https://raw.githubusercontent.com/nf-core/MethylSeq/master/docs/images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](https://raw.githubusercontent.com/nf-core/MethylSeq/master/docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
