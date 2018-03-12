# ![nf-core/methylseq](docs/images/methylseq_logo.png)

[![Build Status](https://travis-ci.org/nf-core/methylseq.svg?branch=master)](https://travis-ci.org/nf-core/methylseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)

### Introduction

**nf-core/methylseq** is a bioinformatics best-practice analysis pipeline used for Methylation (BS-Seq) data analysis.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

### Choice of workflows

There are two separate workflows contained in this repository - one using [Bismark](https://github.com/FelixKrueger/Bismark) and one using [bwa-meth](https://github.com/brentp/bwa-meth) / [MethylDackel](https://github.com/dpryan79/methyldackel). The Bismark pipeline is being actively developed and maintained, the bwa-meth workflow is not _(currently)_. The Nextflow manifest specifies the Bismark pipeline as the default workflow, so the bwa-meth script will be ignored unless explicitly run.

### Documentation
The nf-core/methylseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation and configuration](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)


### Credits
These scripts were originally written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.

* Main authors:
  * Phil Ewels ([@ewels](https://github.com/ewels/))
  * Rickard Hammarén ([@Hammarn](https://github.com/Hammarn/))
* Contributors:
  * Alexander Peltzer ([@apeltzer](https://github.com/apeltzer/))


### Participating Institutes
**nf-core/methylseq** is used by a number of core sequencing and bioinformatics facilities. Some of these are listed below. If you use this pipeline too, please let us know in an issue and we will add you to the list.

<table>
  <tr>
    <td width="200">
      <img src="docs/images/SciLifeLab_logo.png" width=200">
      <img src="docs/images/NGI_logo.png" width=200">
    </td>
    <td>
      SciLifeLab National Genomics Infrastructure (NGI), Sweden <br>
      https://ngisweden.scilifelab.se/
    </td>
  </tr>
  <tr>
    <td width="200"><img src="https://raw.githubusercontent.com/SciLifeLab/NGI-RNAseq/master/docs/images/QBiC_logo.png" width="120"></td>
    <td>
      Quantitative Biology Center (QBiC), Germany <br>
      https://portal.qbic.uni-tuebingen.de/portal/
    </td>
  </tr>
</table>
