//
// Subworkflow with functionality specific to the nf-core/methylseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { paramsHelp              } from 'plugin/nf-schema'
include { completionEmail         } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version // boolean: Display version and exit
    validate_params // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir //  string: The output directory where the results will be saved
    input //  string: Path to input samplesheet
    help // boolean: Display help message and exit
    help_full // boolean: Show the full help message
    show_hidden // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )


    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/methylseq ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/', '')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/methylseq/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )


    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2, genome ->
            if (!fastq_2) {
                return [meta.id, meta + [single_end: true], [fastq_1]]
            }
            else {
                return [meta.id, meta + [single_end: false], [fastq_1, fastq_2]]
            }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .set { ch_samplesheet }
    ch_samplesheet.dump(tag: "ch_samplesheet")

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email //  string: email address
    email_on_fail //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.isBound() ? multiqc_reports.getVal() : [],
            )
        }

        completionSummary(monochrome_logs)

        // Merge topic channel versions into the main versions file
        // This is a workaround until all modules use topic channels
        def mqcVersionsFile = file("${outdir}/pipeline_info/nf_core_methylseq_software_mqc_versions.yml")
        def topicVersionsFile = file("${outdir}/pipeline_info/nf_core_methylseq_topic_versions.yml")
        if (topicVersionsFile.exists() && mqcVersionsFile.exists()) {
            mqcVersionsFile.append(topicVersionsFile.text)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()

    // --combined_index (opt-in, nf-core/methylseq#615) launch-time sanity warnings.
    if (params.combined_index) {
        if (!params.aligner.startsWith('bismark')) {
            log.warn("--combined_index only applies to the bismark aligners; it is ignored for --aligner ${params.aligner}.")
        }
        else if (params.local_alignment) {
            log.warn("--combined_index is ignored when --local_alignment is set (combined-index alignment does not support soft-clipping); using the faithful per-strand path.")
        }
        else if (params.bismark_index) {
            log.warn("--combined_index with a pre-built --bismark_index: the index must have been built with 'bismark_genome_preparation --combined_genome', otherwise alignment will fail (pre-built iGenomes indexes are not compatible).")
        }
    }
    // Building a bwa-mem3 `--meth` index needs far more memory than aligning against it (~100 GB for hg38):
    // point users with large genomes at a pre-built --bwamem3_index.
    // (`params` here is the copy taken at include time, so test for --bwamem3_index with containsKey)
    if (params.aligner == 'bwamem3' && !(params.containsKey('bwamem3_index') && params.bwamem3_index) && params.fasta) {
        def fasta_bytes = 0L
        try { fasta_bytes = file(params.fasta).size() } catch (Exception _e) { }
        if (fasta_bytes > 300_000_000L) {
            log.warn(
                "No --bwamem3_index given, so the pipeline will build a bwa-mem3 --meth index from --fasta (${Math.round(fasta_bytes / 1e6)} MB). " +
                "For human-sized genomes this needs roughly 100 GB of memory: consider building it once with `bwa-mem3 index --meth` and passing it with --bwamem3_index " +
                "(see https://nf-co.re/methylseq/usage). Use --save_reference to keep the index built by this run."
            )
        }
    }
    // --pbat, --non_directional and the non-directional presets (--zymo, --single_cell) only change alignment for Bismark
    def non_directional_options = ['pbat', 'non_directional', 'zymo', 'single_cell'].findAll { option -> params[option] }.collect { option -> "--${option}" }
    if (non_directional_options) {
        // bwa-mem3 `--meth` seeds read 1 as C->T and read 2 as G->A, so it only handles directional libraries
        if (params.aligner == 'bwamem3') {
            error(
                "--aligner bwamem3 only supports directional libraries: ${non_directional_options.join(', ')} not supported. " +
                "Use --aligner bismark for PBAT or non-directional libraries."
            )
        }
        if (params.aligner == 'bwameth') {
            log.warn("bwa-meth aligns directional libraries only: ${non_directional_options.join(', ')} only change trimming with --aligner bwameth. Use --aligner bismark for PBAT or non-directional libraries.")
        }
        if ((params.taps || params.aligner == 'bwamem') && params.taps_caller == 'rastair') {
            log.warn("rastair 0.8 (--taps_caller rastair) assumes directional libraries: ${non_directional_options.join(', ')} are ignored when calling. --taps_caller rastair2 guesses each read's orientation instead.")
        }
    }
    // TAPS calling needs a FASTA index and an aligner without bisulfite methylation extraction
    if (params.taps && params.aligner.startsWith('bismark')) {
        error("--taps is not supported with --aligner ${params.aligner}: use --aligner bwamem3 (recommended), bwameth or bwamem.")
    }
    // MethylDackel writes methylKit files instead of bedGraphs, which the targeted analysis reads
    if (params.methyl_kit && params.run_targeted_sequencing && params.aligner in ['bwameth', 'bwamem3'] && !params.taps) {
        error("--methyl_kit can't be combined with --run_targeted_sequencing for --aligner ${params.aligner}: the targeted analysis reads MethylDackel bedGraphs, which --methyl_kit replaces.")
    }
    // Options that only apply to some steps: say so instead of ignoring them
    def taps_run = params.taps || params.aligner == 'bwamem'
    def ignored_options = [
        [!taps_run && (params.trim_OT || params.trim_OB), '--trim_OT/--trim_OB only apply to TAPS data (--taps or --aligner bwamem); they are ignored'],
        [params.rastair2_mbias_region && !(taps_run && params.taps_caller == 'rastair2'), '--rastair2_mbias_region only applies with --taps_caller rastair2 on TAPS data; it is ignored'],
        [params.collecthsmetrics && !params.run_targeted_sequencing, '--collecthsmetrics only applies with --run_targeted_sequencing; it is ignored'],
        [params.bamqc_regions_file && !params.run_qualimap, '--bamqc_regions_file only applies with --run_qualimap; it is ignored'],
        [params.run_targeted_sequencing && taps_run, 'TAPS runs write no bedGraphs, so --run_targeted_sequencing only runs Picard CollectHsMetrics (with --collecthsmetrics)'],
        [params.slamseq && !params.aligner.startsWith('bismark'), '--slamseq only applies to the bismark aligners; it is ignored'],
        [params.slamseq && params.aligner.startsWith('bismark') && params.containsKey('bismark_index') && params.bismark_index, '--slamseq only changes the index the pipeline builds, so it has no effect with --bismark_index'],
    ]
    ignored_options.findAll { entry -> entry[0] }.each { entry -> log.warn("${entry[1]}.") }
    // rastair 2 and newer is distributed under a non-commercial license: require an explicit opt-in.
    if (params.taps_caller == 'rastair2') {
        if (!(params.taps || params.aligner == 'bwamem')) {
            log.warn("--taps_caller rastair2 only applies to TAPS data (--taps or --aligner bwamem); it is ignored.")
        }
        else if (!params.accept_rastair2_license) {
            error(
                "--taps_caller rastair2 uses rastair version 2, which is free for non-commercial research use only. " +
                "Commercial use requires a license from Oxford University Innovation (enquiries@innovation.ox.ac.uk, reference 24811). " +
                "Run `rastair license` for the full terms. Set --accept_rastair2_license to confirm your use complies, or use --taps_caller rastair."
            )
        }
        else {
            log.warn(
                "rastair version 2 is free for non-commercial research use only; commercial use requires a license from Oxford University Innovation " +
                "(enquiries@innovation.ox.ac.uk, reference 24811). Publications using its results must include the notice \"Copyright © 2025, University of Oxford\"."
            )
        }
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [metas[0], fastqs]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
//
// bwa-mem3 `--meth` index as [ meta, index ]: a directory or .tar.gz, or an index prefix
// (e.g. s3://bucket/genomes/hg38.fa) whose <prefix>.* and <prefix>.meth.* files are staged together
//
def bwamem3IndexChannel(index) {
    if (!index) {
        return []
    }
    def index_path = file(index)
    // Only archives go to UNTAR: an index prefix may itself end in .gz (bwa-mem3 index genome.fa.gz writes genome.fa.gz.amb, ...)
    if (index.toString().endsWith('.tar.gz') || index_path.isDirectory()) {
        return channel.of([[id: index_path.baseName], index_path])
    }
    // `mem --meth` reads only .amb/.ann/.pac from the plain prefix; it seeds from the .meth index
    def index_files = ['amb', 'ann', 'pac', 'meth.amb', 'meth.ann', 'meth.bwt.2bit.64', 'meth.pac'].collect { ext ->
        file("${index}.${ext}", checkIfExists: true)
    }
    return channel.of([[id: index_path.name], index_files])
}

def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" + "  Currently, the available genome keys are:\n" + "  ${params.genomes.keySet().join(", ")}\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // FIXME Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    def citation_text = ["Tools used in the workflow included:", "FastQC (Andrews 2010),", "Trim Galore! (Krueger)", "Bismark (Krueger 2011)", "bwa-meth (Pedersen 2014)", "Picard (Broad Institute 2019)", "Qualimap (Okonechnikov 2015)", "Preseq (Daley 2013)", "MultiQC (Ewels et al. 2016)", "."].join(' ').trim()
}

def toolBibliographyText() {
    // FIXME Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).",
        "Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354",
        "https://www.bioinformatics.babraham.ac.uk/projects/trim_galore",
        "Felix Krueger, Simon R. Andrews, Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications, Bioinformatics, Volume 27, Issue 11, 1 June 2011, Pages 1571–1572, doi: /10.1093/bioinformatics/btr167",
        "Pedersen, Brent S. and Eyring, Kenneth and De, Subhajyoti and Yang, Ivana V. and Schwartz, David A. Fast and accurate alignment of long bisulfite-seq reads, arXiv:1401.1129, doi: 10.48550/arXiv.1401.1129",
        "Picard Tools, Broad Institute. <http://broadinstitute.github.io/picard/>",
        "Konstantin Okonechnikov, Ana Conesa, Fernando García-Alcalde, Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data, Bioinformatics, Volume 32, Issue 2, 15 January 2016, Pages 292–294, doi: 10.1093/bioinformatics/btv566",
        "Daley, T., Smith, A. Predicting the molecular complexity of sequencing libraries. Nat Methods 10, 325–327 (2013). doi: 10.1038/nmeth.2375</li>",
    ].join('</li> <li>').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
