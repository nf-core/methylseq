#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
             B S - S E Q   M E T H Y L A T I O N   B E S T - P R A C T I C E
========================================================================================
 New Methylation (BS-Seq) Best Practice Analysis Pipeline. Started June 2016.
 #### Homepage / Documentation
 https://github.com/SciLifeLab/NGI-MethylSeq
 #### Authors
 Phil Ewels <phil.ewels@scilifelab.se>
----------------------------------------------------------------------------------------
*/


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1

// Configurable variables
params.genome = false
params.bismark_index = params.genome ? params.genomes[ params.genome ].bismark ?: false : false
params.reads = "data/*_{1,2}.fastq.gz"
params.outdir = './results'
params.nodedup = false
params.relaxMismatches = false
params.numMismatches = 0.2
// -0.2 will allow a penalty of bp * -0.2
// This is ~ 20 mismatches or ~ 3 indels of 1-2 bp for a 100bp read (or a combination thereof)

// Validate inputs
if( params.bismark_index ){
    bismark_index = file(params.bismark_index)
    if( !bismark_index.exists() ) exit 1, "Bismark index not found: ${params.bismark_index}"
} else {
    exit 1, "No reference genome specified! Please use --genome or --bismark_index"
}

params.pbat = false
params.single_cell = false
params.epignome = false
params.accel = false
params.cegx = false
if(params.pbat){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 0
    params.three_prime_clip_r2 = 0
} else if(params.single_cell){
    params.clip_r1 = 9
    params.clip_r2 = 9
    params.three_prime_clip_r1 = 0
    params.three_prime_clip_r2 = 0
} else if(params.epignome){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 6
    params.three_prime_clip_r2 = 6
} else if(params.accel){
    params.clip_r1 = 10
    params.clip_r2 = 15
    params.three_prime_clip_r1 = 10
    params.three_prime_clip_r2 = 10
} else if(params.cegx){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 2
    params.three_prime_clip_r2 = 2
} else {
    params.clip_r1 = 0
    params.clip_r2 = 0
    params.three_prime_clip_r1 = 0
    params.three_prime_clip_r2 = 0
}

def single

log.info "===================================="
log.info " NGI-MethylSeq : Bisulfite-Seq Best Practice v${version}"
log.info "===================================="
log.info "Reads          : ${params.reads}"
log.info "Genome         : ${params.genome}"
log.info "Bismark Index  : ${params.bismark_index}"
log.info "Current home   : $HOME"
log.info "Current user   : $USER"
log.info "Current path   : $PWD"
log.info "Script dir     : $baseDir"
log.info "Working dir    : $workDir"
log.info "Output dir     : ${params.outdir}"
log.info "===================================="
log.info "Deduplication  : ${params.nodedup ? 'No' : 'Yes'}"
if(params.rrbs){        log.info "RRBS Mode      : On" }
if(params.pbat){        log.info "Trim Profile   : PBAT" }
if(params.single_cell){ log.info "Trim Profile   : Single Cell" }
if(params.epignome){    log.info "Trim Profile   : Epignome" }
if(params.accel){       log.info "Trim Profile   : Accel" }
if(params.cegx){        log.info "Trim Profile   : CEGX" }
log.info "Output dir     : ${params.outdir}"
log.info "Trim R1        : ${params.clip_r1}"
log.info "Trim R2        : ${params.clip_r2}"
log.info "Trim 3' R1     : ${params.three_prime_clip_r1}"
log.info "Trim 3' R2     : ${params.three_prime_clip_r2}"
log.info "Config Profile : ${workflow.profile}"
log.info "===================================="

// Validate inputs
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: -1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_fastqc; read_files_trimming }

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    set val(name), file(reads) from read_files_fastqc
    
    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    
    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy'
    
    input:
    set val(name), file(reads) from read_files_trimming
    
    output:
    file '*fq.gz' into trimmed_reads
    file '*trimming_report.txt' into trimgalore_results
    
    script:
    single = reads instanceof Path
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    rrbs = params.rrbs ? "--rrbs" : ''
    if (single) {
        """
        trim_galore --gzip $rrbs $c_r1 $tpc_r1 $reads
        """
    } else {
        """
        trim_galore --paired --gzip $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    }
}

/*
 * STEP 3 - align with Bismark
 */
process bismark_align {
    tag "$trimmed_reads"
    publishDir "${params.outdir}/bismark/aligned", mode: 'copy'
    
    input:
    file index from bismark_index
    file trimmed_reads
    
    output:
    file '*.bam' into bam, bam_2
    file '*report.txt' into bismark_align_log_1, bismark_align_log_2, bismark_align_log_3
    file '*.{png,gz}' into bismark_align_results
    file '*.{fq, fastq}' into bismark_unmapped
    
    script:
    pbat = params.pbat ? "--pbat" : ''
    non_directional = params.single_cell || params.non_directional ? "--non_directional" : ''
    unmapped = params.unmapped ? "--unmapped" : ''
    mismatches = params.relaxMismatches ? "--score_min L,0,-${params.numMismatches}" : ''
    if (single) {
        """
        bismark --bam $pbat $non_directional $unmapped $mismatches $index $trimmed_reads
        """
    } else {
        """
        bismark --bam --dovetail $pbat $non_directional $unmapped $mismatches $index -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]}
        """
    }
}

/*
 * STEP 4 - Bismark deduplicate
 */
if (params.nodedup) {
    bam_dedup = bam
} else {
    process bismark_deduplicate {
        tag "$bam"
        publishDir "${params.outdir}/bismark/deduplicated", mode: 'copy'
        
        input:
        file bam
        
        output:
        file '*deduplicated.bam' into bam_dedup
        file '*.deduplication_report.txt' into bismark_dedup_log_1, bismark_dedup_log_2, bismark_dedup_log_3
        file '*.{png,gz}' into bismark_dedup_results
        
        script:
        if (single) {
            """
            deduplicate_bismark -s --bam $bam
            """
        } else {
            """
            deduplicate_bismark -p --bam $bam
            """
        }
    }
}

/*
 * STEP 5 - Bismark methylation extraction
 */
process bismark_methXtract {
    tag "$bam_dedup"
    publishDir "${params.outdir}/bismark/methylation", mode: 'copy'
    
    input:
    file bam_dedup
    
    output:
    file '*.splitting_report.txt' into bismark_splitting_report_1, bismark_splitting_report_2, bismark_splitting_report_3
    file '*.M-bias.txt' into bismark_mbias_1, bismark_mbias_3, bismark_mbias_3
    file '*.{png,gz}' into bismark_methXtract_results
    
    script:
    if (single) {
        """
        bismark_methylation_extractor \\
            --multi ${task.cpus} \\
            --buffer_size ${task.memory} \\
            --bedGraph \\
            --counts \\
            --gzip \\
            -s \\
            --report \\
            $bam_dedup
        """
    } else {
        """
        bismark_methylation_extractor \\
            --multi ${task.cpus} \\
            --buffer_size ${task.memory} \\
            --ignore_r2 2 \\
            --ignore_3prime_r2 2 \\
            --bedGraph \\
            --counts \\
            --gzip \\
            -p \\
            --no_overlap \\
            --report \\
            $bam_dedup
        """
    }
}


/*
 * STEP 6 - Bismark Sample Report
 */
process bismark_report {
    publishDir "${params.outdir}/bismark/summaries", mode: 'copy'
    
    input:
    file bismark_align_log_1
    file bismark_dedup_log_1
    file bismark_splitting_report_1
    file bismark_mbias_1
    
    output:
    file '*{html,txt}' into bismark_reports_results
    
    script:
    """
    bismark2report \\
        --alignment_report $bismark_align_log_1 \\
        --dedup_report $bismark_dedup_log_1 \\
        --splitting_report $bismark_splitting_report_1 \\
        --mbias_report $bismark_mbias_1
    """
}

/*
 * STEP 7 - Bismark Summary Report
 */
process bismark_summary {
    publishDir "${params.outdir}/bismark", mode: 'copy'
    
    input:
    file bam_2.flatten().toList()
    file bismark_align_log_2.flatten().toList()
    file bismark_dedup_log_2.flatten().toList()
    file bismark_splitting_report_2.flatten().toList()
    file bismark_mbias_2.flatten().toList()
    
    output:
    file '*{html,txt}' into bismark_summary_results
    
    script:
    """
    bismark2summary .
    """
}

/*
 * STEP 7 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    
    input:
    file ('fastqc/*') from fastqc_results.flatten().toList()
    file ('trimgalore/*') from trimgalore_results.flatten().toList()
    file ('bismark/*') from bismark_align_log_3.flatten().toList()
    file ('bismark/*') from bismark_dedup_log_3.flatten().toList()
    file ('bismark/*') from bismark_splitting_report_3.flatten().toList()
    file ('bismark/*') from bismark_mbias_3.flatten().toList()
    file ('bismark/*') from bismark_reports_results.flatten().toList()
    file ('bismark/*') from bismark_summary_results.flatten().toList()
    
    output:
    file '*multiqc_report.html'
    file '*multiqc_data'
    
    script:
    """
    multiqc -f .
    """
}
