# NGI-MethylSeq Output

NGI-MethylSeq is the new RNA-seq Best Practice pipeline used by the [National Genomics Infrastructure](https://ngisweden.scilifelab.se/) at [SciLifeLab](https://www.scilifelab.se/platforms/ngi/) in Stockholm, Sweden.

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

Note that NGI-MethylSeq contains two workflows - one for Bismark, one for bwa-meth. This document describes the output from the default Bismark workflow.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [TrimGalore](#trimgalore) - adapter trimming
* [Bismark](#bismark)
  * [Alignment](#alignment) - aligning reads to reference genome
  * [Deduplication](#deduplication) - deduplicating reads
  * [Methylation Extraction](#methylation-extraction) - calling cytosine methylation steps
  * [Report](#report) - single-sample analysis report
  * [Summary](#summary) - multi-sample analysis summary report
* [Qualimap](#qualimap) - QC package for genome alignments
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## TrimGalore
The NGI-MethylSeq BP 2.0 pipeline uses [TrimGalore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for removal of adapter contamination and trimming of low quality regions. TrimGalore uses [Cutadapt](https://github.com/marcelm/cutadapt) for adapter trimming and runs FastQC after it finishes.

MultiQC reports the percentage of bases removed by TrimGalore in the _General Statistics_ table, along with a line plot showing where reads were trimmed.

**Output directory: `results/trim_galore`**

Contains FastQ files with quality and adapter trimmed reads for each sample, along with a log file describing the trimming.

* `sample_val_1.fq.gz`, `sample_val_2.fq.gz`
  * Trimmed FastQ data, reads 1 and 2.
  * NB: Only saved if `--saveTrimmed` has been specified.
* `logs/sample_val_1.fq.gz_trimming_report.txt`
  * Trimming report (describes which parameters that were used)
* `FastQC/sample_val_1_fastqc.zip`
  * FastQC report for trimmed reads

Single-end data will have slightly different file names and only one FastQ file per sample.

## Bismark
The default NGI-MethylSeq workflow uses [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) to process raw sequencing reads into cytosine methylation calls.

### Alignment
Bismark converts all Cytosines to Thymine _in-silico_ and then aligns against a three-letter reference genome. It produces a BAM file of genomic alignments.

**Output directory: `results/bismark_alignments/`**

* `sample.bam`
  * Aligned reads in BAM format.
  * Only saved if `--saveAlignedIntermediates`, `--nodedup` or `--rrbs` is specified when running the pipeline.
* `logs/sample_PE_report.txt`
  * Log file giving summary statistics about alignment.
* `unmapped/unmapped_reads_1.fq.gz`, `unmapped/unmapped_reads_2.fq.gz`
  * Unmapped reads in FastQ format.
  * Only saved if `--unmapped` specified when running the pipeline.

### Deduplication
This step removes alignments with identical mapping position to avoid technical duplication in the results. Note that it is skipped if `--saveAlignedIntermediates`, `--nodedup` or `--rrbs` is specified when running the pipeline.

**Output directory: `results/bismark_deduplicated/`**

* `deduplicated.bam`
  * BAM file with only unique alignments.
* `logs/deduplication_report.txt`
  * Log file giving summary statistics about deduplication.

### Methylation Extraction
The bismark methylation extractor tool takes a BAM file aligned with Bismark and generates files containing methylation information about cytosines. It produces a few different output formats, described below.

Note that the output may vary a little depending on whether you specify `--comprehensive` or `--non_directional` when running the pipeline.

Filename abbreviations stand for the following reference alignment strands:

* `OT` - original top strand
* `OB` - original bottom strand
* `CTOT` - complementary to original top strand
* `CTOB` - complementary to original bottom strand

(`CTOT` and `CTOB` are not aligned unless `--non_directional` specified).

**Output directory: `results/bismark_methylation_calls/`**

* `methylation_calls/XXX_context_sample.txt.gz`
  * Individual methylation calls, sorted into files according to cytosine context.
* `methylation_coverage/sample.bismark.cov.gz`
  * Coverage text file summarising cytosine methylation values.
* `bedGraph/sample.bedGraph.gz`
  * Methylation statuses in [bedGraph](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) format, with 0-based genomic start and 1- based end coordinates.
* `m-bias/sample.M-bias.txt`
  * QC data showing methylation bias across read lengths. See the [bismark documentation](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html#m-bias-plot) for more information.
* `logs/sample_splitting_report.txt`
  * Log file giving summary statistics about methylation extraction.

### Report
Bismark generates a HTML report describing results for each sample.

**Output directory: `results/bismark_reports`**

### Summary
Bismark generates summary text and HTML reports giving an overview of results for all samples in a project.

**Output directory: `results/bismark_summary`**

## Qualimap

[Qualimap BamQC](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#bam-qc) is a general-use quality-control tool that generates a number of statistics about aligned BAM files. It's not specific to bisulfite data, but it produces several useful stats - for example, insert size and coverage statistics.

**Output directory: `results/qualimap`**

* `sample/qualimapReport.html`
  * Qualimap HTML report
* `sample/genome_results.txt`, `sample/raw_data_qualimapReport/*.txt`
  * Text-based statistics that can be loaded into downstream programs


## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
