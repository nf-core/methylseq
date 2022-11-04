/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMethylseq.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.known_splices ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Trimming / kit presets
if(params.pbat){
    params.clip_r1 = 9
    params.clip_r2 = 9
    params.three_prime_clip_r1 = 9
    params.three_prime_clip_r2 = 9
}
else if( params.single_cell ){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 6
    params.three_prime_clip_r2 = 6
}
else if( params.epignome ){
    params.clip_r1 = 8
    params.clip_r2 = 8
    params.three_prime_clip_r1 = 8
    params.three_prime_clip_r2 = 8
}
else if( params.accel || params.zymo ){
    params.clip_r1 = 10
    params.clip_r2 = 15
    params.three_prime_clip_r1 = 10
    params.three_prime_clip_r2 = 10
}
else if( params.cegx ){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 2
    params.three_prime_clip_r2 = 2
}
else if( params.em_seq ){
    params.maxins = 1000
    params.clip_r1 = 8
    params.clip_r2 = 8
    params.three_prime_clip_r1 = 8
    params.three_prime_clip_r2 = 8
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

// Aligner: bismark or bismark_hisat
if( params.aligner =~ /bismark/ ){
    include { BISMARK as ALIGNER } from '../subworkflows/local/bismark'
}
// Aligner: bwameth
else if ( params.aligner == 'bwameth' ){
    include { BWAMETH as ALIGNER } from '../subworkflows/local/bwameth'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { TRIMGALORE      } from '../modules/nf-core/trimgalore/main'
include { QUALIMAP_BAMQC  } from '../modules/nf-core/qualimap/bamqc/main'
include { PRESEQ_LCEXTRAP } from '../modules/nf-core/preseq/lcextrap/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METHYLSEQ {

    versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    versions = versions.mix(INPUT_CHECK.out.versions)



    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    versions = versions.mix(FASTQC.out.versions.first())

    /*
     * MODULE: Run TrimGalore!
     */
    if (!params.skip_trimming) {
        TRIMGALORE(INPUT_CHECK.out.reads)
        reads = TRIMGALORE.out.reads
        versions = versions.mix(TRIMGALORE.out.versions.first())
    } else {
        reads = INPUT_CHECK.out.reads
    }



    /*
     * SUBWORKFLOW: Align reads, deduplicate and extract methylation with Bismark
     */
    ALIGNER (reads)
    versions = versions.mix(ALIGNER.out.versions.unique{ it.baseName })

    /*
     * MODULE: Qualimap BamQC
     */
    QUALIMAP_BAMQC (
        ALIGNER.out.dedup,
        []
    )
    versions = versions.mix(QUALIMAP_BAMQC.out.versions.first())

    /*
     * MODULE: Run Preseq
     */
    PRESEQ_LCEXTRAP (
        ALIGNER.out.bam
    )
    versions = versions.mix(PRESEQ_LCEXTRAP.out.versions.first().ifEmpty(null))

    CUSTOM_DUMPSOFTWAREVERSIONS (
        versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowMethylseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowMethylseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
        ch_methods_description = Channel.value(methods_description)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{ it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.log.collect{ it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ALIGNER.out.mqc.ifEmpty([]))
        if (!params.skip_trimming) {
            ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{ it[1] })
        }
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] }.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.collect().ifEmpty([]),
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
        versions    = versions.mix(MULTIQC.out.versions)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
