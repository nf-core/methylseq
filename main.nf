#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
             B S - S E Q   M E T H Y L A T I O N   B E S T - P R A C T I C E
========================================================================================
 New Methylation (BS-Seq) Best Practice Analysis Pipeline. Started June 2016.
 @Authors
 Phil Ewels <phil.ewels@scilifelab.se>
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow main.nf

 Pipeline variables can be configured with the following command line options:
 --genome [ID] (default: GRCh37)
 --index [path] (default: set by genome ID in config)
 --reads [path] (default: data/*{_1,_2}*.fastq.gz)
 --outdir [path] (default: ./results)
 --name [str] (default: BS-Seq Best Practice)

 For example:
 $ nextflow main.nf --reads 'path/to/data/sample_*_{1,2}.fq.gz' --genome GRCm38
---------------------------------------------------------------------------------------
The pipeline can determine whether the input data is single or paired end. This relies on
specifying the input files correctly. For paired en data us the example above, i.e.
'sample_*_{1,2}.fastq.gz'. Without the glob {1,2} (or similiar) the data will be treated
as single end.
----------------------------------------------------------------------------------------
 Pipeline overview:
 - FastQC - read quility control
 - Trim Galore! - trimming
 - Bismark - align
 - Bismark - deduplication
 - Bismark - methylation extraction
 - Bismark - sample report
 - Bismark - summary report
 - MultiQC
----------------------------------------------------------------------------------------
*/



/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1

// Configurable variables
params.genome = 'GRCh37'
params.index = params.genomes[ params.genome ].bismark
params.reads = "data/*{_1,_2}*.fastq.gz"
params.outdir = './results'

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
log.info "Reads        : ${params.reads}"
log.info "Genome       : ${params.genome}"
log.info "Index        : ${params.index}"
log.info "Current home : $HOME"
log.info "Current user : $USER"
log.info "Current path : $PWD"
log.info "Script dir   : $baseDir"
log.info "Working dir  : $workDir"
log.info "Output dir   : ${params.outdir}"
log.info "===================================="
if(params.pbat){        log.info "Trim Profile : PBAT" }
if(params.single_cell){ log.info "Trim Profile : Single Cell" }
if(params.epignome){    log.info "Trim Profile : Epignome" }
if(params.accel){       log.info "Trim Profile : Accel" }
if(params.cegx){        log.info "Trim Profile : CEGX" }
log.info "Output dir   : ${params.outdir}"
log.info "Trim R1      : ${params.clip_r1}"
log.info "Trim R2      : ${params.clip_r2}"
log.info "Trim 3' R1   : ${params.three_prime_clip_r1}"
log.info "Trim 3' R2   : ${params.three_prime_clip_r2}"
log.info "===================================="

// Validate inputs
index = file(params.index)
if( !index.exists() ) exit 1, "Missing Bismark index: $index"

/*
 * Create a channel for read files - groups based on shared prefixes
 */
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { path ->
        def prefix = readPrefix(path, params.reads)
        tuple(prefix, path)
    }
    .groupTuple(sort: true)
    .set { read_files }

read_files.into { read_files_fastqc; read_files_trimming }


/*
 * STEP 1 - FastQC
 */

process fastqc {
    tag "$prefix"
    
    module 'bioinfo-tools'
    module 'FastQC'
    
    memory { 2.GB * task.attempt }
    time { 4.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    set val(prefix), file(reads:'*') from read_files_fastqc
    
    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    
    """
    fastqc $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */

process trim_galore {
    tag "$prefix"
    
    module 'bioinfo-tools'
    module 'TrimGalore'
    
    cpus 3
    memory { 3.GB * task.attempt }
    time { 16.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/trim_galore", mode: 'copy'
    
    input:
    set val(prefix), file(reads:'*') from read_files_trimming
    
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
        trim_galore --gzip $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
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
    
    module 'bioinfo-tools'
    module 'samtools'
    module 'bismark'
    
    cpus 6
    memory { 32.GB * task.attempt }
    time  { 36.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/bismark/aligned", mode: 'copy'
    
    input:
    file index
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
    if (single) {
        """
        bismark --bam $pbat $non_directional $unmapped $index $trimmed_reads
        """
    } else {
        """
        bismark --bam $pbat $non_directional $unmapped $index -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]}
        """
    }
}


/*
 * STEP 4 - Bismark deduplicate
 */

process bismark_deduplicate {
    tag "$bam"
    
    module 'bioinfo-tools'
    module 'samtools'
    module 'bismark'
    
    memory { 32.GB * task.attempt }
    time  { 12.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'
   
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

/*
 * STEP 5 - Bismark methylation extraction
 */

process bismark_methXtract {
    tag "$bam_dedup"
    
    module 'bioinfo-tools'
    module 'samtools'
    module 'bismark'
    
    cpus 4
    memory { 8.GB * task.attempt }
    time  { 8.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'
    
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
    module 'bioinfo-tools'
    module 'bismark'
    
    memory '2GB'
    time '1h'
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/bismark/summaries", mode: 'copy'
    
    input:
    file bismark_align_log_1
    file bismark_dedup_log_1
    file bismark_splitting_report_1
    file bismark_mbias_1
    
    output:
    file '*{html,txt}' into bismark_reports_results
    
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
    module 'bioinfo-tools'
    module 'bismark'
    
    memory '2GB'
    time '1h'
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/bismark", mode: 'copy'
    
    input:
    file bam_2.toList()
    file bismark_align_log_2.toList()
    file bismark_dedup_log_2.toList()
    file bismark_splitting_report_2.toList()
    file bismark_mbias_2.toList()
    
    output:
    file '*{html,txt}' into bismark_summary_results
    
    """
    bismark2summary .
    """
}


/*
 * STEP 7 - MultiQC
 */

process multiqc {
    module 'bioinfo-tools'
    // Don't load MultiQC module here as overwrites environment installation.
    // Load env module in process instead if multiqc command isn't found.
    
    memory '4GB'
    time '2h'
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    
    input:
    file ('fastqc/*') from fastqc_results.toList()
    file ('trimgalore/*') from trimgalore_results.toList()
    file ('bismark/*') from bismark_align_log_3.toList()
    file ('bismark/*') from bismark_dedup_log_3.toList()
    file ('bismark/*') from bismark_splitting_report_3.toList()
    file ('bismark/*') from bismark_mbias_3.toList()
    file ('bismark/*') from bismark_reports_results.toList()
    file ('bismark/*') from bismark_summary_results.toList()
    
    output:
    file '*multiqc_report.html'
    file '*multiqc_data'
   
    """
    # Load MultiQC with environment module if not already in PATH
    type multiqc >/dev/null 2>&1 || { module load MultiQC; };
    multiqc -f .
    """
}


/*
 * Helper function, given a file Path
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 *
 * For example:
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 *
 * Returns:
 *   'file_alpha'
 */
def readPrefix( Path actual, template ) {

    final fileName = actual.getFileName().toString()

    def filePattern = template.toString()
    int p = filePattern.lastIndexOf('/')
    if( p != -1 ) filePattern = filePattern.substring(p+1)
    if( !filePattern.contains('*') && !filePattern.contains('?') )
        filePattern = '*' + filePattern

    def regex = filePattern
        .replace('.','\\.')
        .replace('*','(.*)')
        .replace('?','(.?)')
        .replace('{','(?:')
        .replace('}',')')
        .replace(',','|')

    def matcher = (fileName =~ /$regex/)
    if( matcher.matches() ) {
        def end = matcher.end(matcher.groupCount() )
        def prefix = fileName.substring(0,end)
        while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') )
            prefix=prefix[0..-2]
        return prefix
    }
    return fileName
}
