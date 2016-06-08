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

single='null'

log.info "===================================="
log.info " RNAbp : RNA-Seq Best Practice v${version}"
log.info "===================================="
log.info "Reads        : ${params.reads}"
log.info "Genome       : ${params.genome}"
log.info "Index        : ${params.index}"
log.info "Annotation   : ${params.gtf}"
log.info "Current home : $HOME"
log.info "Current user : $USER"
log.info "Current path : $PWD"
log.info "Script dir   : $baseDir"
log.info "Working dir  : $workDir"
log.info "Output dir   : ${params.outdir}"
log.info "===================================="

// Validate inputs
index = file(params.index)
if( !index.exists() ) exit 1, "Missing Bismark index: ${index}"

/*
 * Create a channel for read files 
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
 
read_files.into  { read_files_fastqc; read_files_trimming }


/*
 * STEP 1 - FastQC
 */

process fastqc {
    tag "$name"

    module 'bioinfo-tools'
    module 'FastQC'

    memory { 2.GB * task.attempt }
    time { 1.h * task.attempt }

    errorStrategy { task.exitStatus == 140 ? 'retry' : 'warning' }
    maxRetries 3
  
    publishDir "$params.outdir/fastqc"

    input:
    set val(name), file(reads:'*') from read_files_fastqc

    output:
    file '*_fastqc.html' into fastqc_html
    file '*_fastqc.zip' into fastqc_zip

    """
    fastqc -q ${reads}
    """
}


/*
 * STEP 2 - Trim Galore!
 */

process trim_galore {
    tag "$name"

    module 'bioinfo-tools'
    module 'FastQC'
    module 'cutadapt'
    module 'TrimGalore'

    cpus 3
    memory { 3.GB * task.attempt }
    time { 4.h * task.attempt }
    
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
    maxRetries 3
    publishDir "$params.outdir/trim_galore"

    input:
    set val(name), file(reads:'*') from read_files_trimming
    
    output:
    file '*fq.gz' into trimmed_reads
    file '*trimming_report.txt' 

    script:
    single = reads instanceof Path
    if( single ) {
        """
        trim_galore --gzip --fastqc_args "-q" $reads
        """
    } else {
        """
        trim_galore --paired --gzip --fastqc_args "-q" $reads
        """
    }
}

/*
 * STEP 3 - align with Bismark
 */

process bismark_align {
    tag "$trimmed_reads"
    
    module 'bioinfo-tools'
    module 'bismark'
    
    cpus 8
    memory { 32.GB * task.attempt }
    time  { 8.h * task.attempt }
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    publishDir "$params.outdir/bismark"
    
    input:
    file index
    file trimmed_reads
    
    output:
    file '*.bam' into bam
    file '*.{txt,png,gz}'
    
    script:
    if( single ) {
        """
        bismark --bam $index $trimmed_reads
        """
    } else {
        """
        bismark --bam $index -1 $trimmed_reads[0] -2 $trimmed_reads[1]
        """
    }
}


/*
 * STEP 4 - Bismark deduplicate
 */

process bismark_deduplicate {
    tag "$bam"
    
    module 'bioinfo-tools'
    module 'bismark'
    memory { 32.GB * task.attempt }
    time  {4.h * task.attempt }
    
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'warning' }
    maxRetries 3
   
    publishDir "$params.outdir/bismark" 
    
    input:
    file bam
    
    output:
    file '*deduplicated.bam' into bam_dedup
    file '*.{txt,png,gz}'
 
    script:
    if( single ) {
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
    module 'bismark'
    cpus 4
    memory { 8.GB * task.attempt }
    time  {6.h * task.attempt }
    
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'warning' }
    maxRetries 3
   
    publishDir "$params.outdir/bismark" 
    
    input:
    file bam_dedup
    
    output:
    file '*splitting_report.txt' into methXtract_done
    file '*.{txt,png,gz}'
 
    script:
    if( single ) {
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
    module 'bismark'
    
    memory '2GB'   
    time '1h'

    publishDir "$params.outdir/bismark"    
  
    errorStrategy 'ignore'
 
    input:
    file methXtract_done
    
    output:
    file '*{html,txt}' into bismark_summary_done
   
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
    time '4h'

    publishDir "$params.outdir/MultiQC"    
  
    errorStrategy 'ignore'
 
    input:
    file bismark_summary_done
    
    output:
    file 'multiqc_report.html'
    file '*multiqc_data'
   
    """
    multiqc -f $PWD/results
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
    println(fileName) 
    return fileName
}
