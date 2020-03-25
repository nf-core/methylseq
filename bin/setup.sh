#!/usr/bin/env bash

## edit this
##### required #####
## samtools fai-indexed reference
export BISCUIT_REFERENCE="$1"

assets_dir="$2"
#this_dir=$(dirname ${BASH_SOURCE[0]})
##### optional #####
## use <unset> if the file is nonexistent, the corresponding
## QC section will be skipped
## CpGs
export BISCUIT_CPGBED="$assets_dir/cpg.bed.gz"
## CpG islands
export BISCUIT_CGIBED="$assets_dir/cgi.bed.gz"
## repeat masker bed file
export BISCUIT_RMSK="$assets_dir/rmsk.bed.gz"
## merged exon bed file
export BISCUIT_EXON="$assets_dir/exon.bed.gz"
## genes
export BISCUIT_GENE="$assets_dir/genes.bed.gz"
## locations for the top 100bp bins in GC content
export BISCUIT_TOPGC_BED="$assets_dir/windows100bp.gc_content.top10p.bed.gz"
## locations for the bottom 100bp bins in GC content
export BISCUIT_BOTGC_BED="$assets_dir/windows100bp.gc_content.bot10p.bed.gz"

### QC operations to perform ###
export BISCUIT_QC_BASECOV=true
export BISCUIT_QC_DUPLICATE=true
export BISCUIT_QC_CPGCOV=true
export BISCUIT_QC_CPGDIST=true
export BISCUIT_QC_CGICOV=true
export BISCUIT_QC_UNIFORMITY=true
export BISCUIT_QC_CPGUNIF=true
export BISCUIT_QC_BSCONV=true
export BISCUIT_QC_CGICOV=true
export BISCUIT_QC_MAPPING=true
export BISCUIT_QC_BETAS=true
