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
params.name = false
params.project = false
params.email = false
params.genome = false
params.bismark_index = params.genome ? params.genomes[ params.genome ].bismark ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.saveReference = false
params.saveTrimmed = false
params.saveAlignedIntermediates = false
params.reads = "data/*_R{1,2}.fastq.gz"
params.outdir = './results'
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.notrim = false
params.nodedup = false
params.unmapped = false
params.non_directional = false
params.comprehensive = false
params.relaxMismatches = false
params.numMismatches = 0.6
// 0.6 will allow a penalty of bp * -0.6
// For 100bp reads, this is -60. Mismatches cost -6, gap opening -5 and gap extension -2
// So -60 would allow 10 mismatches or ~ 8 x 1-2bp indels
// Bismark default is 0.2 (L,0,-0.2), Bowtie2 default is 0.6 (L,0,-0.6)

// Validate inputs
if( params.bismark_index ){
    bismark_index = Channel
        .fromPath(params.bismark_index)
        .ifEmpty { exit 1, "Bismark index not found: ${params.bismark_index}" }
}
else if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
else {
    exit 1, "No reference genome specified! Please use --genome, --bismark_index or --fasta"
}
multiqc_config = file(params.multiqc_config)

// Validate inputs
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

params.rrbs = false
params.pbat = false
params.single_cell = false
params.epignome = false
params.accel = false
params.zymo = false
params.cegx = false
if(params.pbat){
    params.clip_r1 = 6
    params.clip_r2 = 9
    params.three_prime_clip_r1 = 6
    params.three_prime_clip_r2 = 9
} else if(params.single_cell){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 6
    params.three_prime_clip_r2 = 6
} else if(params.epignome){
    params.clip_r1 = 8
    params.clip_r2 = 8
    params.three_prime_clip_r1 = 8
    params.three_prime_clip_r2 = 8
} else if(params.accel || params.zymo){
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

/*
 * Create a channel for input read files
 */
params.singleEnd = false
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_trimming }

log.info "=================================================="
log.info " NGI-MethylSeq : Bisulfite-Seq Best Practice v${version}"
log.info "=================================================="
def summary = [:]
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Data Type']      = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']         = params.genome
if(params.bismark_index) summary['Bismark Index'] = params.bismark_index
else if(params.fasta)    summary['Fasta Ref'] = params.fasta
if(params.rrbs) summary['RRBS Mode'] = 'On'
if(params.relaxMismatches) summary['Mismatch Func'] = "L,0,-${params.numMismatches} (Bismark default = L,0,-0.2)"
if(params.notrim)       summary['Trimming Step'] = 'Skipped'
if(params.pbat)         summary['Trim Profile'] = 'PBAT'
if(params.single_cell)  summary['Trim Profile'] = 'Single Cell'
if(params.epignome)     summary['Trim Profile'] = 'Epignome'
if(params.accel)        summary['Trim Profile'] = 'Accel'
if(params.cegx)         summary['Trim Profile'] = 'CEGX'
summary['Trim R1'] = params.clip_r1
summary['Trim R2'] = params.clip_r2
summary["Trim 3' R1"] = params.three_prime_clip_r1
summary["Trim 3' R2"] = params.three_prime_clip_r2
summary['Deduplication']  = params.nodedup ? 'No' : 'Yes'
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
summary['Save Unmapped']  = params.unmapped ? 'Yes' : 'No'
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Directional Mode'] = params.non_directional ? 'No' : 'Yes'
summary['All C Contexts'] = params.comprehensive ? 'Yes' : 'No'
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = (workflow.profile == 'standard' ? 'UPPMAX' : workflow.profile)
if(params.project) summary['UPPMAX Project'] = params.project
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "========================================="


/*
 * PREPROCESSING - Build Bismark index
 */
if(!params.bismark_index && fasta){
    process makeBismarkIndex {
        tag fasta
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "Bisulfite_Genome" into bismark_index

        script:
        """
        bismark_genome_preparation ./
        """
    }
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

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
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(name), file(reads) from read_files_trimming

        output:
        set val(name), file('*fq.gz') into trimmed_reads
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        rrbs = params.rrbs ? "--rrbs" : ''
        non_directional = params.rrbs && params.non_directional ? "--non_directional" : ''
        if (params.singleEnd) {
            """
            trim_galore --fastqc --gzip $rrbs $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}

/*
 * STEP 3 - align with Bismark
 */
process bismark_align {
    tag "$name"
    publishDir "${params.outdir}/bismark_alignments", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fq.gz") > 0) "unmapped/$filename"
            else if (filename.indexOf(".bam") == -1) "logs/$filename"
            else params.saveAlignedIntermediates || params.nodedup || params.rrbs ? filename : null
        }

    input:
    file index from bismark_index
    set val(name), file(reads) from trimmed_reads

    output:
    file "*.bam" into bam, bam_2
    file "*report.txt" into bismark_align_log_1, bismark_align_log_2, bismark_align_log_3
    if(params.unmapped){ file "*.fq.gz" into bismark_unmapped }

    script:
    pbat = params.pbat ? "--pbat" : ''
    non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
    unmapped = params.unmapped ? "--unmapped" : ''
    mismatches = params.relaxMismatches ? "--score_min L,0,-${params.numMismatches}" : ''
    if (params.singleEnd) {
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
if (params.nodedup || params.rrbs) {
    bam_dedup = bam
} else {
    process bismark_deduplicate {
        tag "${bam.baseName}"
        publishDir "${params.outdir}/bismark_deduplicated", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

        input:
        file bam

        output:
        file "${bam.baseName}.deduplicated.bam" into bam_dedup, bam_dedup_qualimap
        file "${bam.baseName}.deduplication_report.txt" into bismark_dedup_log_1, bismark_dedup_log_2, bismark_dedup_log_3

        script:
        if (params.singleEnd) {
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
    publishDir "${params.outdir}/bismark_methylation_calls", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("splitting_report.txt") > 0) "logs/$filename"
            else if (filename.indexOf("M-bias") > 0) "m-bias/$filename"
            else if (filename.indexOf(".cov") > 0) "methylation_coverage/$filename"
            else if (filename.indexOf("bedGraph") > 0) "bedGraph/$filename"
            else "methylation_calls/$filename"
        }

    input:
    file bam from bam_dedup

    output:
    file "${bam.baseName}_splitting_report.txt" into bismark_splitting_report_1, bismark_splitting_report_2, bismark_splitting_report_3
    file "${bam.baseName}.M-bias.txt" into bismark_mbias_1, bismark_mbias_2, bismark_mbias_3
    file '*.{png,gz}' into bismark_methXtract_results

    script:
    ignore_r2 = params.rrbs ? "--ignore_r2 2" : ''
    comprehensive = params.comprehensive ? '--comprehensive --merge_non_CpG' : ''
    if (params.singleEnd) {
        """
        bismark_methylation_extractor $comprehensive \\
            --multi ${task.cpus} \\
            --buffer_size ${task.memory.toGiga()}G \\
            $ignore_r2 \\
            --bedGraph \\
            --counts \\
            --gzip \\
            -s \\
            --report \\
            $bam
        """
    } else {
        """
        bismark_methylation_extractor $comprehensive \\
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
    file ('*') from bam_2.collect()
    file ('*') from bismark_align_log_2.collect()
    file ('*') from bismark_dedup_log_2.collect()
    file ('*') from bismark_splitting_report_2.collect()
    file ('*') from bismark_mbias_2.collect()

    output:
    file '*{html,txt}' into bismark_summary_results

    script:
    """
    bismark2summary
    """
}

/*
 * STEP 8 - Qualimap
 */
process qualimap {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    input:
    file bam from bam_dedup_qualimap

    output:
    file "${bam.baseName}_qualimap" into qualimap_results

    script:
    gcref = params.genome == 'GRCh37' ? '-gd HUMAN' : ''
    gcref = params.genome == 'GRCm38' ? '-gd MOUSE' : ''
    """
    samtools sort $bam -o ${bam.baseName}.sorted.bam
    qualimap bamqc $gcref \\
        -bam ${bam.baseName}.sorted.bam \\
        -outdir ${bam.baseName}_qualimap \\
        --collect-overlap-pairs \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}
    """
}

/*
 * STEP 9 - MultiQC
 */
process multiqc {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    echo true

    input:
    file multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('bismark/*') from bismark_align_log_3.collect()
    file ('bismark/*') from bismark_dedup_log_3.collect()
    file ('bismark/*') from bismark_splitting_report_3.collect()
    file ('bismark/*') from bismark_mbias_3.collect()
    file ('bismark/*') from bismark_reports_results.collect()
    file ('bismark/*') from bismark_summary_results.collect()
    file ('qualimap/*') from qualimap_results.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*multiqc_data"

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    """
    multiqc -f -c $multiqc_config . 2>&1
    """
}




/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Build the e-mail subject and header
    def subject = "NGI-MethylSeq Pipeline Complete: $workflow.runName"
    subject += "\nContent-Type: text/html"

    // Set up the e-mail variables
    def email_fields = [:]
    email_fields['version'] = version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container

    // Render the e-mail HTML template
    def f = new File("$baseDir/assets/summary_email.html")
    def engine = new groovy.text.GStringTemplateEngine()
    def template = engine.createTemplate(f).make(email_fields)
    def email_html = template.toString()

    // Send the HTML e-mail
    if (params.email) {
        [ 'mail', '-s', subject, params.email ].execute() << email_html
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_f = new File( output_d, "pipeline_report.html" )
    output_f.withWriter { w ->
        w << email_html
    }

}
