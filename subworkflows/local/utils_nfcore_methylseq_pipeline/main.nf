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
    if (params.slamseq) {
        error("--slamseq was removed in 5.0.0: Bismark 3.x has no SLAM-seq mode.")
    }
    genomeExistsError()

    // `params` here lacks the --genome attributes main.nf sets later; resolve the FASTA as main.nf does,
    // where a fasta param that is set, even to null, wins over --genome
    def fasta = params.containsKey('fasta') ? params.fasta : getGenomeAttribute('fasta')
    if (!fasta) {
        error("No reference FASTA: give --fasta, or a --genome with one (a fasta set to null in a config or params file overrides --genome). Every aligner needs it, even with a pre-built index such as --bismark_index or --bwameth_index.")
    }
    // Bismark finds the FASTA in its genome folder by these extensions only
    def fasta_name = fasta.toString().tokenize('/').last()
    if (bismarkBwamem3() && !(fasta_name.replaceFirst(/\.gz$/, '') ==~ /.+\.(fa|fasta)$/)) {
        error("Bismark's 5-Base/TAPS mode (--aligner bismark) needs a FASTA named *.fa or *.fasta (optionally .gz): rename ${fasta_name}.")
    }

    // --relax_mismatches sets an L,... --score_min, which Bowtie 2 --local rejects
    if (params.aligner == 'bismark' && params.local_alignment && params.relax_mismatches) {
        error("--relax_mismatches can't be combined with --local_alignment for --aligner bismark: Bowtie 2 --local needs a G,... --score_min. Drop --relax_mismatches (Bowtie 2's local default G,20,8 applies) or use --aligner bismark_hisat.")
    }

    // TAPS and Illumina 5-Base read methylated C as T; with --aligner bismark, Bismark aligns with bwa-mem3 and calls
    if (params.taps && params.illumina_5base) {
        error("--taps and --illumina_5base can't be combined: pick the one that matches the library. Both run the same alignment.")
    }
    if (isMcToT() && params.aligner == 'bismark_hisat') {
        error(
            "${params.taps ? '--taps' : '--illumina_5base'} is not supported with --aligner bismark_hisat: Bismark's 5-Base/TAPS mode aligns with bwa-mem3. " +
            "Use --aligner bismark for Bismark calls, or --aligner bwamem3 for rastair calls."
        )
    }
    if (bismarkBwamem3()) {
        def unsupported = [
            'non_directional', 'pbat', 'zymo', 'single_cell', 'rrbs', 'combined_index', 'nomeseq', 'em_seq', 'accel',
            'local_alignment', 'relax_mismatches', 'minins', 'maxins', 'known_splices', 'run_methurator', 'with_umi',
        ].findAll { option -> params[option] }.collect { option -> "--${option}" }
        if (unsupported) {
            error(
                "${unsupported.join(', ')} can't be combined with ${params.taps ? '--taps' : '--illumina_5base'} --aligner bismark: Bismark's 5-Base/TAPS mode " +
                "takes paired-end, directional libraries and aligns with bwa-mem3. " +
                "The rastair route (--aligner bwamem3) is the alternative where it supports the option."
            )
        }
        if (params.taps) {
            log.warn("--taps with --aligner bismark: Bismark aligns with bwa-mem3 and makes the calls; rastair is not run. Use --aligner bwamem3 for rastair calls (a config or profile that sets `taps = true` leaves the aligner at bismark).")
        }
        else {
            log.info("--illumina_5base with --aligner bismark: Bismark aligns with bwa-mem3 and makes the calls. Use --aligner bwamem3 for rastair calls.")
        }
    }
    else if (params.illumina_5base) {
        log.info("--illumina_5base with --aligner ${params.aligner}: aligned and called as --taps would be (${params.aligner == 'bwamem3' ? 'bwa-mem3 --meth=taps' : params.aligner}, then --taps_caller ${params.taps_caller}).")
    }
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
    // Building a bwa-mem3 `--meth` index needs far more memory than aligning against it:
    // point users with large genomes at a pre-built --bwamem3_index.
    def bwamem3_index = params.containsKey('bwamem3_index') ? params.bwamem3_index : getGenomeAttribute('bwamem3')
    if ((params.aligner == 'bwamem3' || bismarkBwamem3()) && !bwamem3_index && fasta) {
        def fasta_bytes = 0L
        try { fasta_bytes = file(fasta).size() } catch (Exception _e) { }
        // A gzipped FASTA is about a third of its uncompressed size, which is what the index build scales with
        def gzipped = fasta.toString().endsWith('.gz')
        def uncompressed_bytes = gzipped ? fasta_bytes * 3 : fasta_bytes
        if (uncompressed_bytes > 300_000_000L) {
            def size_text = gzipped ? "about ${Math.round(uncompressed_bytes / 1e6)} MB uncompressed" : "${Math.round(fasta_bytes / 1e6)} MB"
            log.warn(
                "No --bwamem3_index given, so the pipeline will build a bwa-mem3 --meth index from the FASTA (${size_text}). " +
                "The index build requests 560 MB of memory per 10 MB of uncompressed FASTA, about 175 GB for hg38, and twice that on a retry: " +
                "consider building it once with `bwa-mem3 index --meth --meth-only` and passing it with --bwamem3_index (see https://nf-co.re/methylseq/usage). " +
                "Use --save_reference to keep the index built by this run."
            )
        }
    }
    // --pbat, --non_directional and the non-directional presets (--zymo, --single_cell) only change alignment for Bismark
    def non_directional_options = ['pbat', 'non_directional', 'zymo', 'single_cell'].findAll { option -> params[option] }.collect { option -> "--${option}" }
    if (non_directional_options) {
        // bwa-mem3 `--meth` seeds read 1 as C->T and read 2 as G->A, so it only handles directional libraries
        if (params.aligner == 'bwamem3') {
            error(
                "--aligner bwamem3 only supports directional libraries: ${non_directional_options.join(', ')} not supported." +
                (isMcToT() ? '' : " Use --aligner bismark for PBAT or non-directional libraries.")
            )
        }
        if (params.aligner == 'bwameth') {
            log.warn("bwa-meth aligns directional libraries only: ${non_directional_options.join(', ')} only change trimming with --aligner bwameth.${isMcToT() ? '' : ' Use --aligner bismark for PBAT or non-directional libraries.'}")
        }
        if (isTaps() && params.taps_caller == 'rastair') {
            log.warn("rastair 0.8 (--taps_caller rastair) assumes directional libraries: ${non_directional_options.join(', ')} are ignored when calling. --taps_caller rastair2 guesses each read's orientation instead.")
        }
    }
    // MethylDackel writes methylKit files instead of bedGraphs, which the targeted analysis reads
    if (params.methyl_kit && params.run_targeted_sequencing && params.aligner in ['bwameth', 'bwamem3'] && !isTaps()) {
        error("--methyl_kit can't be combined with --run_targeted_sequencing for --aligner ${params.aligner}: the targeted analysis reads MethylDackel bedGraphs, which --methyl_kit replaces.")
    }
    // UMIs are handled by fgumi runall, which drives bwa-mem3 itself
    if (params.with_umi && params.aligner != 'bwamem3') {
        error("--with_umi is only supported with --aligner bwamem3 (got --aligner ${params.aligner}).")
    }
    // The UMIs are extracted after trimming: fgumi runall aligns straight after extraction
    def five_prime_clip = !params.skip_trimming && (params.clip_r1 > 0 || params.clip_r2 > 0 ||
        (!params.skip_trimming_presets && (params.pbat || params.single_cell || params.accel || params.zymo || params.em_seq)))
    if (params.with_umi && params.umi_read_structure && five_prime_clip) {
        error(
            "--with_umi takes inline UMIs from the reads after trimming, so 5' clipping (--clip_r1, --clip_r2 or a library preset) would remove them. " +
            "Use --skip_trimming_presets without --clip_r1/--clip_r2, and skip the clipped bases in --umi_read_structure instead, e.g. '8M10S+T 8M10S+T'."
        )
    }
    if (params.with_umi && params.skip_deduplication) {
        log.warn("--skip_deduplication with --with_umi: the UMIs are moved to the RX tag, but no duplicates are marked.")
    }
    if (params.run_targeted_sequencing && !params.target_regions_file) {
        error("--target_regions_file must be specified when using --run_targeted_sequencing.")
    }
    if (params.run_targeted_sequencing && params.target_regions_file) {
        def panels = targetPanels()
        def repeated = panels.countBy { panel -> panel.name }.findAll { _name, count -> count > 1 }.keySet()
        if (repeated) {
            error("--target_regions_file: each BED needs its own file name, which (without .bed or .bed.gz) labels its results. Used more than once: ${repeated.join(', ')}.")
        }
        // The panel name ends up in file names, MultiQC anchors and MultiQC path filters, where a dot would make one panel's filter match another's files
        def unsafe = panels.size() > 1 ? panels.findAll { panel -> !(panel.name ==~ /[A-Za-z0-9_-]+/) } : []
        if (unsafe) {
            error("--target_regions_file: with several BEDs, each file name without .bed or .bed.gz may only hold letters, digits, '_' and '-'. Rename: ${unsafe.collect { panel -> panel.bed }.join(', ')}.")
        }
        // The per-panel methylation summary IDs are lower case, with '-' as '_'
        def summary_clashes = panels.groupBy { panel -> panel.name.toLowerCase().replace('-', '_') }.findAll { _id, group -> group.size() > 1 }
        if ((params.run_methylation_summary || params.control_contigs) && summary_clashes) {
            error("--target_regions_file: with --run_methylation_summary or --control_contigs, panel names must still differ when lower-cased and with '-' read as '_'. Rename one of: ${summary_clashes.values().collect { group -> group.collect { panel -> panel.name }.join(' and ') }.join('; ')}.")
        }
    }
    if (params.run_methurator && isTaps()) {
        error("--run_methurator is not supported with TAPS data or Illumina 5-Base data (--taps, --illumina_5base or --aligner bwamem): methurator relies on MethylDackel.")
    }
    // Options that only apply to some steps: say so instead of ignoring them
    def methylation_summary = params.run_methylation_summary || params.control_contigs
    def ignored_options = [
        [params.control_fasta && !params.aligner.startsWith('bismark'), "--control_fasta currently only applies with --aligner bismark or bismark_hisat; it is ignored for --aligner ${params.aligner}"],
        [params.control_fasta && bismarkBwamem3(), "--control_fasta is not supported in Bismark's 5-Base/TAPS mode (--illumina_5base or --taps with --aligner bismark): bwa-mem3 aligns genome reads to the separate controls, which inflates their methylation; it is ignored. Put the controls in the reference and use --control_contigs instead"],
        [params.control_reads != schemaDefault('control_reads') && !params.control_fasta, '--control_reads only applies with --control_fasta; it is ignored'],
        [!isTaps() && !bismarkBwamem3() && (params.trim_OT || params.trim_OB), '--trim_OT/--trim_OB only apply to rastair calls (TAPS or 5-Base data with --aligner bwamem3, bwameth or bwamem, or --aligner bwamem); they are ignored'],
        [bismarkBwamem3() && (params.trim_OT || params.trim_OB), '--trim_OT/--trim_OB only apply to rastair calls; with --aligner bismark trim read ends with --ignore_r1, --ignore_r2, --ignore_3prime_r1 or --ignore_3prime_r2. They are ignored'],
        [bismarkBwamem3() && params.taps_caller != 'rastair', "--taps_caller only applies to rastair calls (--aligner bwamem3, bwameth or bwamem); with --aligner bismark, Bismark makes the calls. It is ignored"],
        [bismarkBwamem3() && params.accept_rastair2_license, '--accept_rastair2_license only applies with --taps_caller rastair2 on the rastair route; it is ignored'],
        [params.five_base_baseq && !bismarkBwamem3(), "--five_base_baseq only applies to Bismark's 5-Base/TAPS mode (--illumina_5base or --taps with --aligner bismark)${isTaps() ? '; rastair applies its own base-quality filter' : ''}. It is ignored"],
        [params.rastair2_mbias_region && !(isTaps() && params.taps_caller == 'rastair2'), '--rastair2_mbias_region only applies with --taps_caller rastair2 on the rastair route (--aligner bwamem3, bwameth or bwamem); it is ignored'],
        [(methylation_summary || params.control_fasta) && params.aligner.startsWith('bismark') && params.nomeseq, 'Bismark coverage files mix the cytosine contexts with --nomeseq, so --run_methylation_summary, --control_contigs and --control_fasta are ignored'],
        [methylation_summary && !isTaps() && !params.aligner.startsWith('bismark') && params.methyl_kit, 'MethylDackel writes methylKit files instead of bedGraphs with --methyl_kit, so --run_methylation_summary and --control_contigs are ignored'],
        [params.run_methylsieve && !isTaps(), '--run_methylsieve only applies on the rastair route (TAPS or 5-Base data with --aligner bwamem3, bwameth or bwamem, or --aligner bwamem); it is ignored'],
        [params.umi_read_structure && !params.with_umi, '--umi_read_structure only applies with --with_umi; it is ignored'],
        [params.collecthsmetrics && !params.run_targeted_sequencing, '--collecthsmetrics only applies with --run_targeted_sequencing; it is ignored'],
        [params.bamqc_regions_file && !params.run_qualimap, '--bamqc_regions_file only applies with --run_qualimap; it is ignored'],
        [userParam('bismark_index') && !params.aligner.startsWith('bismark'), "--bismark_index only applies to the bismark aligners; it is ignored for --aligner ${params.aligner}"],
        [userParam('bismark_index') && bismarkBwamem3(), "--bismark_index is not used by Bismark's 5-Base/TAPS mode, which aligns with bwa-mem3 against the FASTA (--bwamem3_index); it is ignored"],
        [userParam('bwameth_index') && params.aligner != 'bwameth', "--bwameth_index only applies to --aligner bwameth; it is ignored for --aligner ${params.aligner}"],
        [userParam('bwamem_index') && params.aligner != 'bwamem', "--bwamem_index only applies to --aligner bwamem${params.aligner == 'bwamem3' ? ' (use --bwamem3_index for bwa-mem3)' : ''}; it is ignored"],
        [userParam('bwamem3_index') && params.aligner != 'bwamem3' && !bismarkBwamem3(), "--bwamem3_index only applies to --aligner bwamem3 and to Bismark's 5-Base/TAPS mode; it is ignored for --aligner ${params.aligner}"],
    ]
    ignored_options.findAll { entry -> entry[0] }.each { entry -> log.warn("${entry[1]}.") }
    // rastair 2 and newer is distributed under a non-commercial license: require an explicit opt-in.
    if (params.taps_caller == 'rastair2') {
        if (!isTaps()) {
            if (!bismarkBwamem3()) {
                log.warn("--taps_caller rastair2 only applies to rastair calls (TAPS or 5-Base data with --aligner bwamem3, bwameth or bwamem, or --aligner bwamem); it is ignored.")
            }
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
// Stop on --control_contigs names that the FASTA index doesn't list: a typo would otherwise drop the control silently
//
def checkControlContigs(fai) {
    def contigs = fai.readLines().collect { line -> line.tokenize('\t')[0] } as Set
    def missing = params.control_contigs.tokenize(',').findAll { contig -> !(contig in contigs) }
    if (missing) {
        error("--control_contigs: ${missing.join(', ')} not in the reference FASTA (${fai.name} lists ${contigs.size()} contigs). Use the FASTA header names, up to the first space, and add the controls to the FASTA the indexes are built from.")
    }
    return true
}

//
// Fail before alignment when a sample can't give fgumi its UMIs
//
def checkUmiInput(meta, reads) {
    def fastqs = reads instanceof List ? reads : [reads]
    if (params.umi_read_structure) {
        def structures = params.umi_read_structure.tokenize(' ')
        if (structures.size() != 2 && structures.size() != fastqs.size()) {
            error("--umi_read_structure '${params.umi_read_structure}' has ${structures.size()} read structure(s), but sample ${meta.id} has ${fastqs.size()} FASTQ file(s): give one per FASTQ file, e.g. '8M+T 8M+T' (single-end samples use the first of two).")
        }
        def used = structures.take(fastqs.size())
        // Without an M segment the sample gets no UMI, and fgumi dedup then drops every read
        if (!used.any { structure -> structure.contains('M') }) {
            error("--umi_read_structure '${params.umi_read_structure}' gives sample ${meta.id} no UMI: its ${fastqs.size() == 1 ? 'single-end reads use' : 'reads use'} '${used.join(' ')}', which has no M segment.")
        }
        // Trim Galore --rrbs clips 2 bp from the 5' end of read 2 of paired-end reads, before the UMIs are extracted
        if (params.rrbs && !params.skip_trimming && fastqs.size() == 2 && used[1].contains('M')) {
            error("--rrbs trims 2 bp from the 5' end of read 2 for paired-end samples such as ${meta.id}, which would remove the start of an inline read 2 UMI ('${used[1]}'). Use a UMI on read 1 only, UMIs in the read names, or --skip_trimming.")
        }
        return
    }
    // fgumi leaves reads without a UMI untagged, and fgumi dedup then drops them: check the first read name
    def header = ''
    fastqs[0].withInputStream { stream ->
        header = new java.util.zip.GZIPInputStream(stream).withReader { reader -> reader.readLine() } ?: ''
    }
    def fields = header.tokenize(' \t')[0].tokenize(':')
    if (fields.size() < 8 || !(fields[-1] ==~ /(?i)[ACGTNR+-]+/)) {
        error(
            "--with_umi without --umi_read_structure reads the UMI from the last colon-separated field of the read name, " +
            "but the first read of sample ${meta.id} has none ('${header}'). " +
            "Set --umi_read_structure for inline UMIs, or drop --with_umi."
        )
    }
}

//
// TAPS or Illumina 5-Base data: methylated C reads as T.
// Keep in sync with conf/modules/bwamem3_mem.config and conf/subworkflows/fastq_align_dedup_umi_bwamem3.config
//
def isMcToT() {
    return params.taps || params.illumina_5base
}

//
// Bismark's 5-Base/TAPS mode: Bismark aligns with bwa-mem3 and calls.
// Keep in sync with conf/modules/bismark_align.config
//
def bismarkBwamem3() {
    return params.aligner == 'bismark' && isMcToT()
}

//
// The rastair route for TAPS and 5-Base data.
// Keep in sync with conf/modules/bedtools_intersect.config
//
def isTaps() {
    return (isMcToT() || params.aligner == 'bwamem') && !params.aligner.startsWith('bismark')
}

//
// A parameter's default in nextflow_schema.json
//
def schemaDefault(name) {
    def schema = new groovy.json.JsonSlurper().parse(file("${projectDir}/nextflow_schema.json"))
    return schema['$defs'].values().findResult { group -> group['properties']?.get(name)?.get('default') }
}

//
// One panel per --target_regions_file BED (a comma-separated list), named after the file without .bed or .bed.gz
//
def targetPanels() {
    return params.target_regions_file.toString().tokenize(',').collect { bed -> bed.trim() }.findAll { bed -> bed }.collect { bed ->
        [name: bed.tokenize('/').last().replaceFirst(/\.bed(\.gz)?$/, ''), bed: bed]
    }
}

//
// One Picard HsMetrics section per panel (Picard names samples after the BAM), spliced into the last module_order set.
// The panels' General Statistics columns share titles and a column ID, so they start hidden.
//
def targetPanelsMultiqcConfig(panel_names, multiqc_configs) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def orders = multiqc_configs.collect { config -> yaml.load(config.text) }.findAll { config -> config instanceof Map && config.module_order }.collect { config -> config.module_order }
    def module_order = orders ? orders.last() : []
    def hsmetrics = 'CollectHsMetrics.coverage_metrics'
    def panel_entries = panel_names.collect { name ->
        [picard: [
            name: "HsMetrics: ${name}".toString(),
            anchor: "picard_hsmetrics_${name}".toString(),
            info: "Picard CollectHsMetrics on the targets in ${name}.".toString(),
            path_filters: ["*.${name}.${hsmetrics}".toString()],
        ]]
    }
    def is_picard = { entry -> entry == 'picard' || (entry instanceof Map && entry.containsKey('picard')) }
    def other_picard = { entry ->
        def config = entry instanceof Map ? new LinkedHashMap(entry.picard ?: [:]) : [:]
        config.path_filters_exclude = (config.path_filters_exclude ?: []) + ["*.${hsmetrics}".toString()]
        [picard: config]
    }
    def position = module_order.findIndexOf(is_picard)
    def new_order = position < 0
        ? module_order + [other_picard.call('picard')] + panel_entries
        : module_order.take(position) + [other_picard.call(module_order[position])] + panel_entries + module_order.drop(position + 1).findAll { entry -> !is_picard.call(entry) }
    def hidden_columns = panel_names.collectEntries { name -> ["HsMetrics: ${name}: HsMetrics".toString(), false] }
    return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([module_order: new_order, table_columns_visible: hidden_columns]))
}

//
// A parameter the user set: --genome attributes are added to `params` after this file's copy of it is taken
//
def userParam(name) {
    return params.containsKey(name) && params[name]
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

    if (bismarkBwamem3() && metas[0].single_end) {
        error("Please check input samplesheet -> sample ${metas[0].id} is single-end: Bismark's 5-Base/TAPS mode (${params.taps ? '--taps' : '--illumina_5base'} --aligner bismark) needs paired-end reads.")
    }

    // Trim Galore's --fastqc_raw stages each input as <sample>[_1|_2]_raw.fastq.gz and refuses to overwrite it
    if (!params.skip_trimming && !params.skip_fastqc && fastqs.size() == 1) {
        def id = metas[0].id
        def aliases = metas[0].single_end ? ["${id}_raw.fastq.gz".toString()] : ["${id}_1_raw.fastq.gz".toString(), "${id}_2_raw.fastq.gz".toString()]
        def clashing = fastqs[0].collect { fastq -> file(fastq).name }.findAll { name -> name in aliases }
        if (clashing) {
            error("Please check input samplesheet -> sample ${id}: rename ${clashing.join(' and ')}, the name Trim Galore uses for its raw FastQC report, or trimming fails.")
        }
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
// bwa-mem3 `--meth` index as [ meta, index ]: a directory or .tar.gz, or an index prefix
// (e.g. s3://bucket/genomes/hg38.fa) whose <prefix>.* and <prefix>.meth.* files are staged together
//
def bwamem3IndexChannel(index) {
    if (!index) {
        return []
    }
    def index_name = index.toString()
    def usage = "--bwamem3_index must be a directory, a .tar.gz archive or an index prefix (see https://nf-co.re/methylseq/usage)"
    if (index_name.contains('?')) {
        error("${usage}: '${index_name}' has a query string, which hides its file name. Use a URL without one.")
    }
    def index_path = file(index)
    // Only archives go to UNTAR: an index prefix may itself end in .gz (bwa-mem3 index genome.fa.gz writes genome.fa.gz.amb, ...)
    if (index_name.endsWith('.tar.gz') || index_path.isDirectory()) {
        return channel.of([[id: index_path.baseName], index_path])
    }
    if (index_name =~ /\.(tgz|tar|tar\.bz2|tbz2|tar\.xz|txz|zip)$/) {
        error("${usage}: '${index_name}' is an archive format the pipeline can't unpack. Repack it as .tar.gz or unpack it and pass the directory.")
    }
    // `mem --meth` reads only .amb/.ann/.pac from the plain prefix; it seeds from the .meth index
    def index_files = ['amb', 'ann', 'pac', 'meth.amb', 'meth.ann', 'meth.bwt.2bit.64', 'meth.pac'].collect { ext -> file("${index}.${ext}") }
    def missing = index_files.findAll { index_file -> !index_file.exists() }
    if (missing) {
        error("${usage}: treated '${index_name}' as an index prefix, but these files are missing: ${missing.collect { index_file -> index_file.name }.join(', ')}.")
    }
    return channel.of([[id: index_path.name], index_files])
}

//
// Exit pipeline if incorrect --genome key provided
//
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
    def citation_text = ["Tools used in the workflow included:", "FastQC (Andrews 2010),", "Trim Galore! (Krueger)", "Bismark (Krueger 2011)", "bwa-meth (Pedersen 2014)", "Picard (Broad Institute 2019)", params.with_umi ? "fgumi (Fulcrum Genomics)," : "", "Qualimap (Okonechnikov 2015)", "Preseq (Daley 2013)", "MultiQC (Ewels et al. 2016)", "."].findAll { citation -> citation }.join(' ').trim()
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
        params.with_umi ? "fgumi, Fulcrum Genomics. <https://github.com/fulcrumgenomics/fgumi>" : "",
        "Konstantin Okonechnikov, Ana Conesa, Fernando García-Alcalde, Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data, Bioinformatics, Volume 32, Issue 2, 15 January 2016, Pages 292–294, doi: 10.1093/bioinformatics/btv566",
        "Daley, T., Smith, A. Predicting the molecular complexity of sequencing libraries. Nat Methods 10, 325–327 (2013). doi: 10.1038/nmeth.2375</li>",
    ].findAll { reference -> reference }.join('</li> <li>').trim()

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
