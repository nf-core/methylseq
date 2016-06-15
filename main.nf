#!/usr/bin/env nextflow

/*
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
params.name = "BS-Seq Best practice"
params.reads = "data/*{_1,_2}*.fastq.gz"
params.outdir = './results'

params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

if(params.pbat){
    params.clip_r1 = 6
    params.clip_r2 = 6
}
if(params.single_cell){
    params.clip_r1 = 9
    params.clip_r2 = 9
}
if(params.epignome){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 6
    params.three_prime_clip_r2 = 6
}
if(params.accel){
    params.clip_r1 = 10
    params.clip_r2 = 15
    params.three_prime_clip_r1 = 10
    params.three_prime_clip_r2 = 10
}
if(params.cegx){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 2
    params.three_prime_clip_r2 = 2
}

single = 'null'

log.info "===================================="
log.info " RNAbp : RNA-Seq Best Practice v${version}"
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
if( !index.exists() ) exit 1, "Missing Bismark index: ${index}"

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
    module 'java'
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
    module 'java'
    module 'FastQC'
    module 'cutadapt'
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
    
    if (single) {
        """
        trim_galore --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    } else {
        """
        trim_galore --paired --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
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
    
    cpus 8
    memory { 32.GB * task.attempt }
    time  { 36.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/bismark", mode: 'copy'
    
    input:
    file index
    file trimmed_reads
    
    output:
    file '*.bam' into bam
    file '*.{txt,png,gz}' into bismark_align_results_1, bismark_align_results_2
    
    script:
    if (single) {
        """
        bismark --bam $index $trimmed_reads
        """
    } else {
        """
        bismark --bam $index -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]}
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
   
    publishDir "${params.outdir}/bismark", mode: 'copy'
    
    input:
    file bam
    
    output:
    file '*deduplicated.bam' into bam_dedup
    file '*.{txt,png,gz}' into bismark_dedup_results_1, bismark_dedup_results_2
    
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
    
    publishDir "${params.outdir}/bismark", mode: 'copy'
    
    input:
    file bam_dedup
    
    output:
    file '*.{txt,png,gz}' into bismark_methXtract_results_1, bismark_methXtract_results_2
    
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
            --ignore_r2 1 \\
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
 * STEP 6 - Bismark Summary
 */

process bismark_summary {
    module 'bioinfo-tools'
    module 'samtools'
    module 'bismark'
    
    memory '2GB'
    time '1h'
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/bismark", mode: 'copy'
    
    input:
    file bismark_align_results_1.toList()
    file bismark_dedup_results_1.toList()
    file bismark_methXtract_results_1.toList()
    
    output:
    file '*{html,txt}' into bismark_summary_results
    
    """
    bismark2summary $PWD/results/bismark
    """
}


/*
 * STEP 7 - MultiQC
 */

process multiqc {
    module 'bioinfo-tools'
    module 'MultiQC'
    
    memory '4GB'
    time '2h'
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    
    input:
    file ('fastqc/*') from fastqc_results.toList()
    file ('trimgalore/*') from trimgalore_results.toList()
    file ('bismark/*') from bismark_align_results_2.toList()
    file ('bismark/*') from bismark_dedup_results_2.toList()
    file ('bismark/*') from bismark_methXtract_results_2.toList()
    file ('bismark/*') from bismark_summary_results.toList()
    
    output:
    file '*multiqc_report.html'
    file '*multiqc_data'
   
    """
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
