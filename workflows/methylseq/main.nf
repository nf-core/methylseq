/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.preview.types = true

include { paramsSummaryMap           } from 'plugin/nf-schema'
include { WRITE_FILE                 } from '../../modules/local/writefile/main'
include { FASTQC                     } from '../../modules/nf-core/fastqc/main'
include { TRIMGALORE                 } from '../../modules/nf-core/trimgalore/main'
include { QUALIMAP_BAMQC             } from '../../modules/nf-core/qualimap/bamqc/main'
include { PRESEQ_LCEXTRAP            } from '../../modules/nf-core/preseq/lcextrap/main'
include { MULTIQC                    } from '../../modules/nf-core/multiqc/main'
include { CAT_FASTQ                  } from '../../modules/nf-core/cat/fastq/main'
include { FASTQ_ALIGN_DEDUP_BISMARK  } from '../../subworkflows/nf-core/fastq_align_dedup_bismark/main'
include { FASTQ_ALIGN_DEDUP_BWAMETH  } from '../../subworkflows/nf-core/fastq_align_dedup_bwameth/main'
include { paramsSummaryMultiqc       } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML     } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { TARGETED_SEQUENCING        } from '../../subworkflows/local/targeted_sequencing'
include { methodsDescriptionText     } from '../../subworkflows/local/utils_nfcore_methylseq_pipeline'
include { validateInputSamplesheet   } from '../../subworkflows/local/utils_nfcore_methylseq_pipeline'

include { Sample } from '../../utils/types.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLSEQ {

    take:
    ch_samples: Channel<Sample>
    val_fasta: Value<Path>?
    val_fasta_index: Value<Path>?
    val_bismark_index: Value<Path>?
    val_bwameth_index: Value<Path>?
    params: MethylseqParams

    main:

    ch_multiqc_files = channel.empty()

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    ch_samples_merged = CAT_FASTQ (
        ch_samples.filter { r -> !r.meta.single_end }
    )
    ch_samples_single = ch_samples.filter { r -> r.meta.single_end }

    ch_fastq = ch_samples_merged.mix(ch_samples_single)

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        ch_fastqc = FASTQC( ch_fastq )
    } else {
        ch_fastqc = channel.empty()
    }

    //
    // MODULE: Run TrimGalore!
    //
    if (!params.skip_trimming) {
        ch_trimmed_fastq = TRIMGALORE( ch_fastq )
        ch_reads = ch_trimmed_fastq.map { r -> r + record(reads: r.trim_reads) }
    } else {
        ch_trimmed_fastq = channel.empty()
        ch_reads = ch_fastq
    }

    //
    // SUBWORKFLOW: Align reads, deduplicate and extract methylation with Bismark
    //

    val_bismark_summary = null

    // Aligner: bismark or bismark_hisat
    if ( params.aligner =~ /bismark/ && val_fasta && val_bismark_index ) {
        //
        // Run Bismark alignment + downstream processing
        //
        bismark = FASTQ_ALIGN_DEDUP_BISMARK (
            ch_reads,
            val_fasta,
            val_bismark_index,
            params.skip_deduplication || params.rrbs,
            params.cytosine_report || params.nomeseq
        )
        ch_alignment  = bismark.results
        val_bismark_summary = bismark.bismark_summary
        ch_aligner_mqc = bismark.multiqc
    }
    // Aligner: bwameth
    else if ( params.aligner == 'bwameth' && val_fasta && val_fasta_index && val_bwameth_index ) {

        bwameth = FASTQ_ALIGN_DEDUP_BWAMETH (
            ch_reads,
            val_fasta,
            val_fasta_index,
            val_bwameth_index,
            params.skip_deduplication || params.rrbs,
            workflow.profile.tokenize(',').intersect(['gpu']).size() >= 1
        )
        ch_alignment  = bwameth.results
        ch_aligner_mqc = bwameth.multiqc
    }
    else {
        error "ERROR: Invalid aligner '${params.aligner}'. Valid options are: 'bismark', 'bismark_hisat', or 'bwameth'"
    }

    //
    // MODULE: Qualimap BamQC
    // skipped by default. to use run with `--run_qualimap` param.
    //
    if(params.run_qualimap) {
        bamqc_regions_file = params.bamqc_regions_file ? file( params.bamqc_regions_file, checkIfExists: true ) : null
        ch_qualimap = QUALIMAP_BAMQC (
            ch_alignment.combine(gff: bamqc_regions_file)
        )
    } else {
        ch_qualimap = channel.empty()
    }

    //
    // MODULE: Targeted sequencing analysis
    // skipped by default. to use run with `--run_targeted_sequencing` param.
    //
    if (params.run_targeted_sequencing){
        if (!params.target_regions_file) {
            error "ERROR: --target_regions_file must be specified when using --run_targeted_sequencing"
        }
        targeted_sequencing = TARGETED_SEQUENCING (
            ch_alignment,
            channel.value(file(params.target_regions_file, checkIfExists: true)),
            val_fasta,
            val_fasta_index,
            params.collecthsmetrics
        )
        ch_targeted_sequencing = targeted_sequencing.results
        val_reference_dict = targeted_sequencing.reference_dict
        val_intervallist = targeted_sequencing.intervallist
    } else {
        ch_targeted_sequencing = channel.empty()
        val_reference_dict = null
        val_intervallist = null
    }

    //
    // MODULE: Preseq LCEXTRAP
    // skipped by default. to use run with `--run_preseq` param.
    //
    if(params.run_preseq) {
        ch_preseq = PRESEQ_LCEXTRAP ( ch_alignment )
    } else {
        ch_preseq = channel.empty()
    }

    ch_results = ch_alignment
        .join(ch_fastqc, by: 'id', remainder: true)
        .join(ch_trimmed_fastq, by: 'id', remainder: true)
        .join(ch_qualimap, by: 'id', remainder: true)
        .join(ch_targeted_sequencing, by: 'id', remainder: true)
        .join(ch_preseq, by: 'id', remainder: true)

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        multiqc_custom_config = params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : null
        multiqc_logo          = params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : null

        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        workflow_summary = record(
            name: 'workflow_summary_mqc.yaml',
            items: [paramsSummaryMultiqc(summary_params)]
        )

        multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        methods_description = record(
            name: 'methods_description_mqc.yaml',
            items: [methodsDescriptionText(multiqc_custom_methods_description)]
        )

        val_versions = softwareVersionsToYAML( channel.topic('versions') )
            .collect()
            .map { items ->
                record(name: 'nf_core_methylseq_software_mqc_versions.yml', items: items.toSorted(), newLine: true)
            }

        ch_collected_files = WRITE_FILE(
            channel.of(workflow_summary, methods_description).mix(val_versions)
        )
        ch_multiqc_files = ch_multiqc_files.mix(ch_collected_files)

        if(params.run_qualimap) {
            ch_multiqc_files = ch_multiqc_files.mix(ch_qualimap.map { r -> r.qualimap_bamqc })
        }
        if (params.run_preseq) {
            ch_multiqc_files = ch_multiqc_files.mix(ch_preseq.map { r -> r.lc_log })
        }
        ch_multiqc_files = ch_multiqc_files.mix(ch_aligner_mqc)
        if (!params.skip_trimming) {
            ch_multiqc_files = ch_multiqc_files.mix(ch_reads.map { r -> r.trim_log })
        }
        if (params.run_targeted_sequencing && params.collecthsmetrics) {
            ch_multiqc_files = ch_multiqc_files.mix(ch_targeted_sequencing.map { r -> r.picard_hsmetrics })
        }
        if (!params.skip_fastqc) {
            ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc.map { r -> r.fastqc_zip })
        }

        val_multiqc_inputs = ch_multiqc_files
            .flatMap()
            .collect()
            .map { multiqc_files ->
                record(
                    multiqc_files: multiqc_files.toSet(),
                    multiqc_config: multiqc_config,
                    extra_multiqc_config: multiqc_custom_config,
                    multiqc_logo: multiqc_logo
                )
            }
        val_multiqc_report = MULTIQC ( val_multiqc_inputs ).map { r -> r.report }
    } else {
        val_multiqc_report = null
    }

    emit:
    results         : Channel<MethylseqResult> = ch_results
    bismark_summary : Value<Record>? = val_bismark_summary
    reference_dict  : Value<Record>? = val_reference_dict
    intervallist    : Value<Record>? = val_intervallist
    multiqc_report  : Value<Path>? = val_multiqc_report
}

record MethylseqParams {
    skip_fastqc: Boolean
    skip_trimming: Boolean
    aligner: String
    skip_deduplication: Boolean
    rrbs: Boolean
    cytosine_report: Boolean
    nomeseq: Boolean
    run_qualimap: Boolean
    bamqc_regions_file: String
    run_targeted_sequencing: Boolean
    target_regions_file: String
    collecthsmetrics: Boolean
    run_preseq: Boolean
    outdir: String
    skip_multiqc: Boolean
    multiqc_config: String
    multiqc_logo: String
    multiqc_methods_description: String
}

record MethylseqResult {
    id: String
    single_end: Boolean

    // fastqc
    fastqc_html: Set<Path>
    fastqc_zip: Set<Path>

    // trimgalore
    trim_reads: List<Path>
    trim_log: List<Path>
    trim_unpaired: List<Path>
    trim_html: List<Path>
    trim_zip: List<Path>

    // alignment (bismark / bwameth)
    bam: Path
    bai: Path

    // bismark
    align_report: Path
    unmapped: Path?
    dedup_report: Path
    coverage2cytosine_coverage: Path
    coverage2cytosine_report: Path
    coverage2cytosine_summary: Path
    methylation_bedgraph: Path
    methylation_calls: Path
    methylation_coverage: Path
    methylation_report: Path
    methylation_mbias: Path
    bismark_report: Path

    // bwameth
    samtools_flagstat: Path
    samtools_stats: Path
    methyldackel_extract_bedgraph: Path
    methyldackel_extract_methylkit: Path
    methyldackel_mbias: Path
    picard_metrics: Path

    // qualimap
    qualimap_bamqc: Path?

    // targeted sequencing
    bedgraph_intersect: Path
    picard_hsmetrics: Path

    // preseq
    lc_extrap: Path?
    lc_log: Path?

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
