nextflow.enable.types = true

include { BismarkArgs            } from '../../subworkflows/nf-core/fastq_align_dedup_bismark'
include { MethyldackelArgs       } from '../../subworkflows/nf-core/bam_methyldackel'

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

def trimgaloreArgs(meta: Record, p: TrimgaloreParams) -> String {
    // Clip presets per protocol: [clip_r1, clip_r2, three_prime_clip_r1, three_prime_clip_r2], 0 = none
    def preset: List<Integer> =
        p.skip_trimming_presets ? [0, 0, 0, 0] :
        p.pbat                  ? [8, 8, 8, 8] :
        p.single_cell           ? [6, 6, 6, 6] :
        p.zymo || p.em_seq      ? [10, 10, 10, 10] :
        p.accel                 ? [10, 15, 10, 10] :
                                  [0, 0, 0, 0]
    return cli([
        fastqc: true,
        rrbs: p.rrbs,
        nextseq: p.nextseq_trim > 0 ? p.nextseq_trim : null,
        length: p.length_trim ?: null,
        clip_r1: p.clip_r1 > 0 ? p.clip_r1 : preset[0] ?: null,
        clip_r2: meta.single_end ? null : p.clip_r2 > 0 ? p.clip_r2 : preset[1] ?: null,
        three_prime_clip_r1: p.three_prime_clip_r1 > 0 ? p.three_prime_clip_r1 : preset[2] ?: null,
        three_prime_clip_r2: meta.single_end ? null : p.three_prime_clip_r2 > 0 ? p.three_prime_clip_r2 : preset[3] ?: null
    ])
}

def bismarkArgs(meta: Record, p: BismarkParams) -> BismarkArgs {
    return record(
        align: cli(bismarkAlignOpts(meta, p)),
        methylation_extractor: cli(bismarkMethylationExtractorOpts(meta, p)),
        coverage2cytosine: cli(['nome-seq': p.nomeseq])
    )
}

def bismarkAlignOpts(meta: Record, p: BismarkParams) -> Map<String,?> {
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
        minins: meta.single_end ? null : p.minins,
        maxins: meta.single_end ? null : p.maxins ?: (p.em_seq ? 1000 : null)
    ]
}

def bismarkMethylationExtractorOpts(meta: Record, p: BismarkParams) -> Map<String,?> {
    def pe = !meta.single_end
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
    return cli([
        hisat2: hisat,
        bowtie2: !hisat,
        slam: p.slamseq,
        // Combined-index build is skipped for --local_alignment (combined rejects --local)
        combined_genome: p.combined_index && !p.local_alignment && p.aligner.startsWith('bismark')
    ])
}

// Picard takes boolean values as strings
def markduplicatesArgs() -> String {
    return cli([ASSUME_SORTED: 'true', REMOVE_DUPLICATES: 'false', VALIDATION_STRINGENCY: 'LENIENT', PROGRAM_RECORD_ID: "'null'", TMP_DIR: 'tmp'])
}

def methyldackelArgs(p: MethyldackelParams) -> MethyldackelArgs {
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
        extract: cli(extract),
        mbias: cli(mbias)
    )
}

def qualimapArgs(p: QualimapParams) -> String {
    return cli([
        '-gd': p.genome?.startsWith('GRCh') ? 'HUMAN' : p.genome?.startsWith('GRCm') ? 'MOUSE' : null
    ])
}

def methuratorArgs(p: MethuratorParams) -> String {
    return cli([
        compute_ci: p.methurator_compute_ci,
        rrbs: p.rrbs,
        'minimum-coverage': p.methurator_minimum_coverage,
        't-max': p.methurator_t_max
    ])
}

def multiqcArgs(p: MultiqcParams) -> String {
    return cli([title: p.multiqc_title ? "\"${p.multiqc_title}\"" : null])
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
    multiqc_title: String?
}
