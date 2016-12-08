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
params.project = false
params.emailAddress = false
params.genome = false
params.bismark_index = params.genome ? params.genomes[ params.genome ].bismark ?: false : false
params.saveReference = false
params.reads = "data/*_R{1,2}.fastq.gz"
params.outdir = './results'
params.notrim = false
params.nodedup = false
params.relaxMismatches = false
params.numMismatches = 0.6
// 0.6 will allow a penalty of bp * -0.6
// For 100bp reads, this is -60. Mismatches cost -6, gap opening -5 and gap extension -2
// Sp -60 would allow 10 mismatches or ~ 8 x 1-2bp indels
// Bismark default is 0.2 (L,0,-0.2), Bowtie2 default is 0.6 (L,0,-0.6)

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

log.info "=================================================="
log.info " NGI-MethylSeq : Bisulfite-Seq Best Practice v${version}"
log.info "=================================================="
log.info "Reads          : ${params.reads}"
log.info "Genome         : ${params.genome}"
log.info "Bismark Index  : ${params.bismark_index}"
log.info "Current home   : $HOME"
log.info "Current user   : $USER"
log.info "Current path   : $PWD"
log.info "Script dir     : $baseDir"
log.info "Working dir    : $workDir"
log.info "Output dir     : ${params.outdir}"
log.info "---------------------------------------------------"
log.info "Deduplication  : ${params.nodedup ? 'No' : 'Yes'}"
if(params.rrbs){            log.info "RRBS Mode      : On" }
if(params.relaxMismatches){ log.info "Mismatch Func  : L,0,-${params.numMismatches} (Bismark default = L,0,-0.2)" }
log.info "---------------------------------------------------"
if(params.notrim){      log.info "Trimming Step  : Skipped" }
if(params.pbat){        log.info "Trim Profile   : PBAT" }
if(params.single_cell){ log.info "Trim Profile   : Single Cell" }
if(params.epignome){    log.info "Trim Profile   : Epignome" }
if(params.accel){       log.info "Trim Profile   : Accel" }
if(params.cegx){        log.info "Trim Profile   : CEGX" }
if(params.clip_r1 > 0)  log.info "Trim R1        : ${params.clip_r1}"
if(params.clip_r2 > 0)  log.info "Trim R2        : ${params.clip_r2}"
if(params.three_prime_clip_r1 > 0) log.info "Trim 3' R1     : ${params.three_prime_clip_r1}"
if(params.three_prime_clip_r2 > 0) log.info "Trim 3' R2     : ${params.three_prime_clip_r2}"
log.info "---------------------------------------------------"
log.info "Config Profile : ${workflow.profile}"
if(params.project) log.info "UPPMAX Project : ${params.project}"
log.info "=================================================="

// Validate inputs
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: -1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
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
if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = []
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy'
        
        input:
        set val(name), file(reads) from read_files_trimming
        
        output:
        set val(name), file('*fq.gz') into trimmed_reads
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
}

/*
 * STEP 3 - align with Bismark
 */
process bismark_align {
    tag "$name"
    publishDir "${params.outdir}/bismark_alignments", mode: 'copy'
    
    input:
    file index from bismark_index
    set val(name), file(reads) from trimmed_reads
    
    output:
    file "*.bam" into bam, bam_2
    file "*report.txt" into bismark_align_log_1, bismark_align_log_2, bismark_align_log_3
    if(params.unmapped){ file "*.fq.gz" into bismark_unmapped }
    
    script:
    pbat = params.pbat ? "--pbat" : ''
    non_directional = params.single_cell || params.non_directional ? "--non_directional" : ''
    unmapped = params.unmapped ? "--unmapped" : ''
    mismatches = params.relaxMismatches ? "--score_min L,0,-${params.numMismatches}" : ''
    if (single) {
        """
        bismark --bam $pbat $non_directional $unmapped $mismatches $index $reads
        """
    } else {
        """
        bismark \\
            --bam \\
            --dovetail \\
            $pbat $non_directional $unmapped $mismatches \\
            $index \\
            -1 ${reads[0]} \\
            -2 ${reads[1]}
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
        tag "${bam.baseName}"
        publishDir "${params.outdir}/bismark_deduplicated", mode: 'copy'
        
        input:
        file bam
        
        output:
        file "${bam.baseName}.deduplicated.bam" into bam_dedup
        file "${bam.baseName}.deduplication_report.txt" into bismark_dedup_log_1, bismark_dedup_log_2, bismark_dedup_log_3
        
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
    tag "${bam.baseName}"
    publishDir "${params.outdir}/bismark_methylation_calls", mode: 'copy'
    
    input:
    file bam from bam_dedup
    
    output:
    file "${bam.baseName}_splitting_report.txt" into bismark_splitting_report_1, bismark_splitting_report_2, bismark_splitting_report_3
    file "${bam.baseName}.M-bias.txt" into bismark_mbias_1, bismark_mbias_2, bismark_mbias_3
    file '*.{png,gz}' into bismark_methXtract_results
    
    script:
    if (single) {
        """
        bismark_methylation_extractor \\
            --multi ${task.cpus} \\
            --buffer_size ${task.memory.toGiga()}G \\
            --ignore_r2 2 \\
            --bedGraph \\
            --counts \\
            --gzip \\
            -s \\
            --report \\
            $bam
        """
    } else {
        """
        bismark_methylation_extractor \\
            --multi ${task.cpus} \\
            --buffer_size ${task.memory.toGiga()}G \\
            --ignore_r2 2 \\
            --ignore_3prime_r2 2 \\
            --bedGraph \\
            --counts \\
            --gzip \\
            -p \\
            --no_overlap \\
            --report \\
            $bam
        """
    }
}


/*
 * STEP 6 - Bismark Sample Report
 */
process bismark_report {
    tag "$name"
    publishDir "${params.outdir}/bismark_reports", mode: 'copy'
    
    input:
    file bismark_align_log_1
    file bismark_dedup_log_1
    file bismark_splitting_report_1
    file bismark_mbias_1
    
    output:
    file '*{html,txt}' into bismark_reports_results
    
    script:
    name = bismark_align_log_1.toString() - ~/(_R1)?(_trimmed|_val_1).+$/
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
    publishDir "${params.outdir}/bismark_summary", mode: 'copy'
    
    input:
    file ('*') from bam_2.flatten().toList()
    file ('*') from bismark_align_log_2.flatten().toList()
    file ('*') from bismark_dedup_log_2.flatten().toList()
    file ('*') from bismark_splitting_report_2.flatten().toList()
    file ('*') from bismark_mbias_2.flatten().toList()
    
    output:
    file '*{html,txt}' into bismark_summary_results
    
    script:
    """
    bismark2summary
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


// E-mail sent upon pipeline completion
workflow.onComplete {
    if(params.emailAddress){
        def subject = "NGI-MethylSeq pipeline completed ${ workflow.success ? 'successfully' : 'with errors' }"

        ['mail', '-s', subject, params.emailAddress].execute() << """

        NGI-MethylSeq pipeline execution summary
        ----------------------------------------
        Starterd at  : ${workflow.start}
        Completed at : ${workflow.complete}
        Duration     : ${workflow.duration}
        Success      : ${workflow.success}
        Launch Dir   : ${workflow.launchDir}
        Work Dir     : ${workflow.workDir}
        Command      : ${workflow.commandLine}
        Resumed      : ${workflow.resume}
        Profile      : ${workflow.profile == 'standard' ? 'UPPMAX' : workflow.profile}
        Nextflow v   : ${nextflow.version}
        Exit status  : ${workflow.exitStatus}
        Error msg    : ${workflow.errorMessage ?: '-'}
        Error report : ${workflow.errorReport ?: '-'}
        """
    }
}
