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

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.1'
try {
  if( ! nextflow.version.matches(">= $nf_required_version") ){
    throw GroovyException('Nextflow version too old')
  }
} catch (all) {
  log.error "====================================================\n" +
            "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

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
    makeBismarkIndex_stderr = Channel.create()
    makeBismarkIndex_done = Channel.create()
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
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
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
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "BismarkIndex" into bismark_index
        file '.command.err' into makeBismarkIndex_stderr
        val 'done' into makeBismarkIndex_done

        script:
        """
        mkdir BismarkIndex
        cp $fasta BismarkIndex/
        bismark_genome_preparation BismarkIndex
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
    stdout fastqc_stdout
    val 'done' into fastqc_done

    script:
    """
    fastqc -q $reads
    fastqc --version
    """
}

/*
 * STEP 2 - Trim Galore!
 */
if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = Channel.create()
    trimgalore_logs = Channel.create()
    trimgalore_logs_done = Channel.create()
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
        file "*trimming_report.txt" into trimgalore_results, trimgalore_logs
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
        val 'done' into trimgalore_logs_done

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
    file "*report.txt" into bismark_align_log_1, bismark_align_log_2, bismark_align_log_3, bismark_align_log_4
    if(params.unmapped){ file "*.fq.gz" into bismark_unmapped }
    val 'done' into bismark_align_done

    script:
    pbat = params.pbat ? "--pbat" : ''
    non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
    unmapped = params.unmapped ? "--unmapped" : ''
    mismatches = params.relaxMismatches ? "--score_min L,0,-${params.numMismatches}" : ''
    if (params.singleEnd) {
        """
        bismark \\
            --bam $pbat $non_directional $unmapped $mismatches \\
            --genome $index \\
            $reads
        """
    } else {
        """
        bismark \\
            --bam $pbat $non_directional $unmapped $mismatches \\
            --genome $index \\
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
        stdout bismark_deduplicate_stdout
        val 'done' into bismark_dedup_done

        script:
        if (params.singleEnd) {
            """
            deduplicate_bismark -s --bam $bam
            deduplicate_bismark --version
            """
        } else {
            """
            deduplicate_bismark -p --bam $bam
            deduplicate_bismark --version
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
    file '.command.err' into bismark_methXtract_stderr
    val 'done' into bismark_methXtract_done

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
    stdout bismark_report_stdout
    val 'done' into bismark_reports_done

    script:
    name = bismark_align_log_1.toString() - ~/(_R1)?(_trimmed|_val_1).+$/
    """
    bismark2report \\
        --alignment_report $bismark_align_log_1 \\
        --dedup_report $bismark_dedup_log_1 \\
        --splitting_report $bismark_splitting_report_1 \\
        --mbias_report $bismark_mbias_1
    bismark2report --version
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
    stdout bismark_summary_stdout
    val 'done' into bismark_summary_done

    script:
    """
    bismark2summary
    bismark2summary --version
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
    stdout qualimap_stdout
    val 'done' into qualimap_done

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
 * Parse software version numbers
 */
software_versions = [
  'Bismark genomePrep': null, 'FastQC': null, 'Trim Galore!': null, 'Bismark': null, 'Bismark Deduplication': null,
  'Bismark methXtract': null, 'Bismark Report': null, 'Bismark Summary': null, 'Qualimap': null
]
makeBismarkIndex_stderr.subscribe { stdout -> software_versions['Bismark genomePrep'] = stdout.getText().find(/Bisulfite Genome Indexer version v(\S+)/) { match, version -> "v$version" } }
fastqc_stdout.subscribe { stdout -> software_versions['FastQC'] = stdout.find(/FastQC v(\S+)/) { match, version -> "v$version" } }
trimgalore_logs.subscribe { stdout -> software_versions['Trim Galore!'] = stdout.getText().find(/Trim Galore version: (\S+)/) { match, version -> "v$version" } }
bismark_align_log_4.subscribe { logfile -> software_versions['Bismark'] = logfile.getText().find(/Bismark report for: .* \(version: v(.+)\)/) { match, version -> "v$version" } }
bismark_deduplicate_stdout.subscribe { stdout -> software_versions['Bismark Deduplication'] = stdout.find(/Deduplicator Version: v(\S+)/) { match, version -> "v$version" } }
bismark_methXtract_stderr.subscribe { stdout -> software_versions['Bismark methXtract'] = stdout.getText().find(/Bismark methylation extractor version v(\S+)/) { match, version -> "v$version" } }
bismark_report_stdout.subscribe { stdout -> software_versions['Bismark Report'] = stdout.find(/bismark2report version: v(\S+)/) { match, version -> "v$version" } }
bismark_summary_stdout.subscribe { stdout -> software_versions['Bismark Summary'] = stdout.find(/bismark2summary version: (\S+)/) { match, version -> "v$version" } }
qualimap_stdout.subscribe { stdout -> software_versions['Qualimap'] = stdout.find(/QualiMap v.(\S+)/) { match, version -> "v$version" } }

process get_software_versions {
    input:
    file makeBismarkIndex from makeBismarkIndex_done.collect()
    file fastqc from fastqc_done.collect()
    file trimgalore from trimgalore_logs_done.collect()
    file bismark_align from bismark_align_done.collect()
    file bismark_deduplicate from bismark_dedup_done.collect()
    file bismark_methXtract from bismark_methXtract_done.collect()
    file bismark_report from bismark_reports_done.collect()
    file bismark_summary from bismark_summary_done.collect()
    file qualimap from qualimap_done.collect()

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo "
    id: 'ngi-rnaseq'
    section_name: 'NGI-RNAseq Software Versions'
    plot_type: 'html'
    description: 'are collected at run time from the software output.'
    data: |
        <dl class=\\"dl-horizontal\\">
${software_versions.collect{ k,v -> "            <dt>$k</dt><dd>${v ?: '<span style=\\"color:#CCCCCC;\\">N/A</a>'}</dd>" }.join("\n")}
        </dl>
    " > software_versions_mqc.yaml
    """
}


/*
 * STEP 9 - MultiQC
 */
process multiqc {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

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
    file ('software_versions/*') from software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*multiqc_data"
    file '.command.err' into multiqc_stderr

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    """
    multiqc -f -c $multiqc_config .
    """
}
multiqc_stderr.subscribe { stdout -> software_versions['MultiQC'] = stdout.getText().find(/This is MultiQC v(\S+)/) { match, version -> "v$version" } }

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[NGI-MethylSeq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[NGI-MethylSeq] FAILED: $workflow.runName"
    }
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
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['software_versions'] = software_versions
    email_fields['software_versions']['Nextflow Version'] = "v$workflow.nextflow.version"
    email_fields['software_versions']['Nextflow Build'] = workflow.nextflow.build
    email_fields['software_versions']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.debug "[NGI-MethylSeq] Sent summary e-mail using sendmail"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.debug "[NGI-MethylSeq] Sendmail failed, failing back to sending summary e-mail using mail"
        }
        log.info "[NGI-MethylSeq] Sent summary e-mail to $params.email"
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[NGI-MethylSeq] Pipeline Complete"
}
