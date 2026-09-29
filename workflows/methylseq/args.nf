nextflow.enable.types = true

include { SampleMeta             } from '../../utils/types.nf'
include { BismarkArgs            } from '../../subworkflows/nf-core/fastq_align_dedup_bismark'
include { BwamethArgs            } from '../../subworkflows/nf-core/fastq_align_dedup_bwameth'
include { BwamemArgs             } from '../../subworkflows/nf-core/fastq_align_dedup_bwamem'
include { MethyldackelArgs       } from '../../subworkflows/nf-core/bam_methyldackel'
include { TargetedSequencingArgs } from '../../subworkflows/local/targeted_sequencing'

/*
 * Resolve tool args:
 * `<tool>_args` samplesheet column > pipeline defaults merged with `--opts.<tool>.<option>`
 */
def toolArgs(tool: String, s: SampleMeta, opts: Map<String,Map<String,?>>, defaults: Map<String,?>) -> String {
    return s.tool_args[tool] ?: runArgs(tool, opts, defaults)
}

/*
 * Resolve args for a run-level tool (not per sample): pipeline defaults merged with `--opts.<tool>.<option>`
 */
def runArgs(tool: String, opts: Map<String,Map<String,?>>, defaults: Map<String,?>) -> String {
    return cli(defaults + cliOpts(opts[tool] ?: [:]))
}

/*
 * Render tool options as CLI args. Boolean -> bare flag or omitted, null -> omitted,
 * any other value -> `<flag> <value>`. The flag is `-k` for a 1-char key, `--key` otherwise,
 * or the key itself if it starts with `-`.
 */
def cli(opts: Map<String,?>) -> String {
    return opts.keySet()
        .collect { k ->
            def v = opts[k]
            def flag = k.startsWith('-') ? k : (k.length() == 1 ? "-${k}" : "--${k}")
            v instanceof Boolean ? (v ? flag : '') : v != null ? "${flag} ${v}" : ''
        }
        .findAll { a -> a != '' }
        .join(' ')
}

// CLI values arrive as strings: 'true'/'false' mean flag on/off
def cliOpts(opts: Map<String,?>) -> Map<String,?> {
    return opts.keySet().inject([:]) { acc, k ->
        def v = "${opts[k]}"
        acc + [(k): opts[k]] + (v == 'true' ? [(k): true] : v == 'false' ? [(k): false] : [:])
    }
}

def trimgaloreOpts(s: SampleMeta, p: TrimgaloreParams) -> Map<String,?> {
    // Clip presets per protocol: [clip_r1, clip_r2, three_prime_clip_r1, three_prime_clip_r2], 0 = none
    def preset: List<Integer> =
        p.skip_trimming_presets ? [0, 0, 0, 0] :
        p.pbat                  ? [8, 8, 8, 8] :
        p.single_cell           ? [6, 6, 6, 6] :
        p.zymo || p.em_seq      ? [10, 10, 10, 10] :
        p.accel                 ? [10, 15, 10, 10] :
                                  [0, 0, 0, 0]
    return [
        fastqc: true,
        rrbs: p.rrbs,
        nextseq: p.nextseq_trim > 0 ? p.nextseq_trim : null,
        length: p.length_trim ?: null,
        clip_r1: p.clip_r1 > 0 ? p.clip_r1 : preset[0] ?: null,
        clip_r2: s.single_end ? null : p.clip_r2 > 0 ? p.clip_r2 : preset[1] ?: null,
        three_prime_clip_r1: p.three_prime_clip_r1 > 0 ? p.three_prime_clip_r1 : preset[2] ?: null,
        three_prime_clip_r2: s.single_end ? null : p.three_prime_clip_r2 > 0 ? p.three_prime_clip_r2 : preset[3] ?: null
    ]
}

def bismarkArgs(s: SampleMeta, p: BismarkParams) -> BismarkArgs {
    return record(
        align: toolArgs('bismark_align', s, p.opts, bismarkAlignOpts(s, p)),
        deduplicate: toolArgs('bismark_deduplicate', s, p.opts, [:]),
        methylation_extractor: toolArgs('bismark_methylationextractor', s, p.opts, bismarkMethylationExtractorOpts(s, p)),
        coverage2cytosine: toolArgs('bismark_coverage2cytosine', s, p.opts, ['nome-seq': p.nomeseq]),
        report: toolArgs('bismark_report', s, p.opts, [:])
    )
}

def bismarkAlignOpts(s: SampleMeta, p: BismarkParams) -> Map<String,?> {
    // Combined-index alignment is incompatible with --local_alignment, so gated off there
    def hisat = p.aligner == 'bismark_hisat'
    def non_directional = p.single_cell || p.non_directional || p.zymo
    def use_combined = p.aligner.startsWith('bismark') && p.combined_index && !p.local_alignment
    return [
        hisat2: hisat,
        bowtie2: !hisat,
        'known-splicesite-infile': hisat && p.known_splices ? "<(hisat2_extract_splice_sites.py ${p.known_splices})" : null,
        pbat: p.pbat,
        non_directional: non_directional,
        combined_index: use_combined,
        combined_index_sequential: use_combined && non_directional,
        unmapped: p.unmapped,
        score_min: p.relax_mismatches ? "L,0,-${p.num_mismatches}" : null,
        local: p.local_alignment,
        minins: s.single_end ? null : p.minins,
        maxins: s.single_end ? null : p.maxins ?: (p.em_seq ? 1000 : null)
    ]
}

def bismarkMethylationExtractorOpts(s: SampleMeta, p: BismarkParams) -> Map<String,?> {
    def pe = !s.single_end
    return [
        comprehensive: p.comprehensive,
        cutoff: p.meth_cutoff,
        CX: p.nomeseq,
        ignore: p.ignore_r1 > 0 ? p.ignore_r1 : null,
        ignore_3prime: p.ignore_3prime_r1 > 0 ? p.ignore_3prime_r1 : null,
        no_overlap: pe && p.no_overlap,
        include_overlap: pe && !p.no_overlap,
        ignore_r2: pe && p.ignore_r2 > 0 ? p.ignore_r2 : null,
        ignore_3prime_r2: pe && p.ignore_3prime_r2 > 0 ? p.ignore_3prime_r2 : null
    ]
}

def bismarkGenomePreparationArgs(p: BismarkParams) -> String {
    def hisat = p.aligner == 'bismark_hisat'
    return runArgs('bismark_genomepreparation', p.opts, [
        hisat2: hisat,
        bowtie2: !hisat,
        slam: p.slamseq,
        // Combined-index build is skipped for --local_alignment (combined rejects --local)
        combined_genome: p.combined_index && !p.local_alignment && p.aligner.startsWith('bismark')
    ])
}

def bwamethArgs(s: SampleMeta, opts: Map<String,Map<String,?>>) -> BwamethArgs {
    return record(
        align: toolArgs('bwameth_align', s, opts, [:]),
        fq2bammeth: toolArgs('parabricks_fq2bammeth', s, opts, ['low-memory': true]),
        markduplicates: toolArgs('picard_markduplicates', s, opts, markduplicatesOpts())
    )
}

def bwamemArgs(s: SampleMeta, opts: Map<String,Map<String,?>>) -> BwamemArgs {
    return record(
        align: toolArgs('bwa_mem', s, opts, [:]),
        fq2bam: toolArgs('parabricks_fq2bam', s, opts, [:]),
        addorreplacereadgroups: toolArgs('picard_addorreplacereadgroups', s, opts, [RGID: 1, RGLB: 'lib1', RGPL: 'illumina', RGPU: 'unit1', RGSM: 'sample1']),
        markduplicates: toolArgs('picard_markduplicates', s, opts, markduplicatesOpts())
    )
}

// Picard takes boolean values as strings
def markduplicatesOpts() -> Map<String,?> {
    return [ASSUME_SORTED: 'true', REMOVE_DUPLICATES: 'false', VALIDATION_STRINGENCY: 'LENIENT', PROGRAM_RECORD_ID: "'null'", TMP_DIR: 'tmp']
}

def methyldackelArgs(s: SampleMeta, p: MethyldackelParams) -> MethyldackelArgs {
    def extract: Map<String,?> = [
        CHG: p.all_contexts,
        CHH: p.all_contexts,
        mergeContext: p.merge_context,
        ignoreFlags: p.ignore_flags,
        methylKit: p.methyl_kit,
        minDepth: p.min_depth > 0 ? p.min_depth : null
    ]
    def mbias: Map<String,?> = [
        CHG: p.all_contexts,
        CHH: p.all_contexts,
        ignoreFlags: p.ignore_flags
    ]
    return record(
        extract: toolArgs('methyldackel_extract', s, p.opts, extract),
        mbias: toolArgs('methyldackel_mbias', s, p.opts, mbias)
    )
}

def targetedSequencingArgs(s: SampleMeta, opts: Map<String,Map<String,?>>) -> TargetedSequencingArgs {
    return record(
        intersect_cov: toolArgs('bedtools_intersect_cov', s, opts, [:]),
        collecthsmetrics: toolArgs('picard_collecthsmetrics', s, opts, [MINIMUM_MAPPING_QUALITY: 20, COVERAGE_CAP: 1000, NEAR_DISTANCE: 500])
    )
}

def qualimapOpts(p: QualimapParams) -> Map<String,?> {
    return [
        '-gd': p.genome?.startsWith('GRCh') ? 'HUMAN' : p.genome?.startsWith('GRCm') ? 'MOUSE' : null
    ]
}

def methuratorOpts(p: MethuratorParams) -> Map<String,?> {
    return [
        compute_ci: p.methurator_compute_ci,
        rrbs: p.rrbs,
        'minimum-coverage': p.methurator_minimum_coverage,
        't-max': p.methurator_t_max
    ]
}

def multiqcArgs(p: MultiqcParams) -> String {
    return runArgs('multiqc', p.opts, [title: p.multiqc_title ? "\"${p.multiqc_title}\"" : null])
}

record TrimgaloreParams {
    rrbs: Boolean
    nextseq_trim: Integer
    length_trim: Integer?
    clip_r1: Integer
    clip_r2: Integer
    three_prime_clip_r1: Integer
    three_prime_clip_r2: Integer
    skip_trimming_presets: Boolean
    pbat: Boolean
    single_cell: Boolean
    accel: Boolean
    zymo: Boolean
    em_seq: Boolean
}

record BismarkParams {
    opts: Map<String,Map<String,?>>
    aligner: String
    known_splices: Path?
    pbat: Boolean
    single_cell: Boolean
    non_directional: Boolean
    zymo: Boolean
    em_seq: Boolean
    slamseq: Boolean
    combined_index: Boolean
    local_alignment: Boolean
    unmapped: Boolean
    relax_mismatches: Boolean
    num_mismatches: Float
    minins: Integer?
    maxins: Integer?
    comprehensive: Boolean
    meth_cutoff: Integer?
    nomeseq: Boolean
    ignore_r1: Integer
    ignore_3prime_r1: Integer
    no_overlap: Boolean
    ignore_r2: Integer
    ignore_3prime_r2: Integer
}

record MethyldackelParams {
    opts: Map<String,Map<String,?>>
    all_contexts: Boolean
    merge_context: Boolean
    ignore_flags: Boolean
    methyl_kit: Boolean
    min_depth: Integer
}

record QualimapParams {
    genome: String?
}

record MethuratorParams {
    methurator_compute_ci: Boolean
    rrbs: Boolean
    methurator_minimum_coverage: String?
    methurator_t_max: Integer?
}

record MultiqcParams {
    opts: Map<String,Map<String,?>>
    multiqc_title: String?
}
