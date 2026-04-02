nextflow.preview.types = true

include { UNTAR                     } from '../../../modules/nf-core/untar/main'
include { GUNZIP                    } from '../../../modules/nf-core/gunzip/main'
include { BISMARK_GENOMEPREPARATION } from '../../../modules/nf-core/bismark/genomepreparation/main'
include { BWAMETH_INDEX             } from '../../../modules/nf-core/bwameth/index/main'
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx/main'

include { MethylseqParams } from '../../../workflows/methylseq/main'

def isGzipped(file: Path) -> Boolean {
    return file.name.endsWith('.gz')
}

workflow FASTA_INDEX_BISMARK_BWAMETH {

    take:
    fasta: Path?
    fasta_index: Path?
    bismark_index: Path?
    bwameth_index: Path?
    use_mem2: Boolean           // generate mem2 index if no index provided, and bwameth is selected
    params: MethylseqParams

    main:

    val_fasta         = null
    val_fasta_index   = null
    val_bismark_index = null
    val_bwameth_index = null

    // Check if fasta file is gzipped and decompress if needed
    if( fasta ) {
        val_fasta = isGzipped(fasta)
            ? GUNZIP( fasta )
            : channel.value(fasta)
    }

    // Aligner: bismark or bismark_hisat
    if( params.aligner =~ /bismark/ ){
        /*
         * Generate bismark index if not supplied
         */
        if (bismark_index) {
            // Handle channel-based bismark index
            val_bismark_index = isGzipped(bismark_index)
                ? UNTAR( bismark_index )
                : channel.value(bismark_index)
        } else {
            val_bismark_index = BISMARK_GENOMEPREPARATION( val_fasta )
        }
    }

    // Aligner: bwameth
    else if ( params.aligner == 'bwameth' ){
        /*
         * Generate bwameth index if not supplied
         */
        if (bwameth_index) {
            // Handle channel-based bwameth index
            val_bwameth_index = isGzipped(bwameth_index)
                ? UNTAR( bwameth_index )
                : channel.value(bwameth_index)
        } else {
            val_bwameth_index = BWAMETH_INDEX( val_fasta, use_mem2 )
        }
    }

    /*
    * Generate fasta index if not supplied for bwameth workflow or picard collecthsmetrics tool
    */
    if (params.aligner == 'bwameth' || params.collecthsmetrics) {
        // already exising fasta index
        if (fasta_index) {
            val_fasta_index = channel.value(fasta_index)
        } else {
            val_faidx_inputs = val_fasta.map { fa -> record(fasta: fa, get_sizes: false) }
            val_fasta_index = SAMTOOLS_FAIDX( val_faidx_inputs ).map { r -> r.fai }
        }
    }

    emit:
    fasta         : Value<Path>? = val_fasta
    fasta_index   : Value<Path>? = val_fasta_index
    bismark_index : Value<Path>? = val_bismark_index
    bwameth_index : Value<Path>? = val_bwameth_index
}
