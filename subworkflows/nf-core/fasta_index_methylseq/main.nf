nextflow.enable.types = true

include { UNTAR as UNTAR_BISMARK    } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_BWAMETH    } from '../../../modules/nf-core/untar/main'
include { GUNZIP                    } from '../../../modules/nf-core/gunzip/main'
include { BISMARK_GENOMEPREPARATION as BISMARK_GENOMEPREPARATION_BOWTIE } from '../../../modules/nf-core/bismark/genomepreparation/main'
include { BISMARK_GENOMEPREPARATION as BISMARK_GENOMEPREPARATION_HISAT } from '../../../modules/nf-core/bismark/genomepreparation/main'
include { BWAMETH_INDEX             } from '../../../modules/nf-core/bwameth/index/main'
include { BWA_INDEX                 } from '../../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx/main'

def isGzipped(file: Path) -> Boolean {
    return file.name.endsWith('.gz')
}

workflow FASTA_INDEX_METHYLSEQ {

    take:
    fasta: Path
    fasta_index: Path?
    bismark_index: Path?
    bwameth_index: Path?
    bwamem_index: Path?
    aligner: String             // bismark, bismark_hisat, bwameth or bwamem
    collecthsmetrics: Boolean   // whether to run picard collecthsmetrics
    methurator: Boolean         // whether to run methurator
    use_mem2: Boolean           // generate mem2 index if no index provided, and bwameth is selected
    genomeprep_args: String     // args for bismark genome preparation

    main:

    val_fasta_index   = null
    val_bismark_index = null
    val_bwameth_index = null
    val_bwamem_index  = null

    // Check if fasta file is gzipped and decompress if needed
    val_fasta = isGzipped(fasta)
        ? GUNZIP( fasta, '' )
        : channel.value(fasta)

    // Aligner: bismark or bismark_hisat
    if( aligner =~ /bismark/ ){
        /*
         * Generate bismark index if not supplied
         */
        if (bismark_index) {
            val_bismark_index = isGzipped(bismark_index)
                ? UNTAR_BISMARK( bismark_index, '' )
                : channel.value(bismark_index)
        } else if( aligner == "bismark_hisat") {
            val_bismark_index = BISMARK_GENOMEPREPARATION_HISAT( val_fasta, genomeprep_args )
        } else {
            val_bismark_index = BISMARK_GENOMEPREPARATION_BOWTIE( val_fasta, genomeprep_args )
        }
    }

    // Aligner: bwameth
    else if ( aligner == 'bwameth' ){
        /*
         * Generate bwameth index if not supplied
         */
        if (bwameth_index) {
            val_bwameth_index = isGzipped(bwameth_index)
                ? UNTAR_BWAMETH( bwameth_index, '' )
                : channel.value(bwameth_index)
        } else {
            val_bwameth_index = BWAMETH_INDEX( val_fasta, use_mem2 )
        }
    }

    else if ( aligner == 'bwamem' ){
        /*
         * Generate BWA index from FASTA file
         */
        if (bwamem_index) {
            val_bwamem_index = isGzipped(bwamem_index)
                ? UNTAR_BISMARK( bwamem_index, '' )
                : channel.value(bwamem_index)
        } else {
            log.info "BWA index not provided. Generating BWA index from FASTA file."
            val_bwamem_index = BWA_INDEX( val_fasta, '' )
        }
    }

    /*
    * Generate fasta index if not supplied for bwameth workflow or picard collecthsmetrics tool or methurator tool
    */
    if (aligner == 'bwameth' || aligner == 'bwamem' || collecthsmetrics || methurator) {
        if (fasta_index) {
            val_fasta_index = channel.value(fasta_index)
        } else {
            log.info "Fasta index not provided. Generating fasta index from FASTA file."
            val_fasta_index = SAMTOOLS_FAIDX( val_fasta, null, false, '' ).map { r -> r.fai }
        }
    }

    emit:
    fasta         : Value<Path>  = val_fasta
    fasta_index   : Value<Path>? = val_fasta_index
    bismark_index : Value<Path>? = val_bismark_index
    bwameth_index : Value<Path>? = val_bwameth_index
    bwamem_index  : Value<Path>? = val_bwamem_index
}
