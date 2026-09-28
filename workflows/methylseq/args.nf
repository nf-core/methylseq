nextflow.enable.types = true

include { SampleMeta             } from '../../utils/types.nf'
include { BismarkArgs            } from '../../subworkflows/nf-core/fastq_align_dedup_bismark'
include { BwamethArgs            } from '../../subworkflows/nf-core/fastq_align_dedup_bwameth'
include { BwamemArgs             } from '../../subworkflows/nf-core/fastq_align_dedup_bwamem'
include { MethyldackelArgs       } from '../../subworkflows/nf-core/bam_methyldackel'
include { TargetedSequencingArgs } from '../../subworkflows/local/targeted_sequencing'

/*
 * Resolve tool args:
 * `<tool>_args` samplesheet column > `--args.<tool>` > pipeline default
 */
def toolArgs(tool: String, s: SampleMeta, args: Map<String,String>, fallback: String) -> String {
    return s.tool_args[tool] ?: args[tool] ?: fallback
}

/*
 * Resolve args for a run-level tool (not per sample): `--args.<tool>` > pipeline default
 */
def runArgs(tool: String, args: Map<String,String>, fallback: String) -> String {
    return args[tool] ?: fallback
}

def trimgaloreArgs(s: SampleMeta, p: TrimgaloreParams) -> String {
    return [
        // Static args
        '--fastqc',

        // Special flags
        p.rrbs ? '--rrbs' : '',
        p.nextseq_trim > 0 ? "--nextseq ${p.nextseq_trim}" : '',
        p.length_trim ? "--length ${p.length_trim}" : '',

        // Trimming - R1
        p.clip_r1 > 0 ? "--clip_r1 ${p.clip_r1}" : (
            p.skip_trimming_presets ? '' : (
                p.pbat ? "--clip_r1 8" : (
                    p.single_cell ? "--clip_r1 6" : (
                        (p.accel || p.zymo || p.em_seq) ? "--clip_r1 10" : ''
                    )
                )
            )
        ),

        // Trimming - R2
        s.single_end ? '' : (
            p.clip_r2 > 0 ? "--clip_r2 ${p.clip_r2}" : (
                p.skip_trimming_presets ? '' : (
                    p.pbat ? "--clip_r2 8" : (
                        p.single_cell ? "--clip_r2 6" : (
                            (p.zymo || p.em_seq) ? "--clip_r2 10" : (
                                p.accel ? "--clip_r2 15" : ''
                            )
                        )
                    )
                )
            )
        ),

        // Trimming - 3' R1
        p.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${p.three_prime_clip_r1}" : (
            p.skip_trimming_presets ? '' : (
                p.pbat ? "--three_prime_clip_r1 8" : (
                    p.single_cell ? "--three_prime_clip_r1 6" : (
                        (p.accel || p.zymo || p.em_seq) ? "--three_prime_clip_r1 10" : ''
                    )
                )
            )
        ),

        // Trimming - 3' R2
        s.single_end ? '' : (
            p.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${p.three_prime_clip_r2}" : (
                p.skip_trimming_presets ? '' : (
                    p.pbat ? "--three_prime_clip_r2 8" : (
                        p.single_cell ? "--three_prime_clip_r2 6" : (
                            (p.accel || p.zymo || p.em_seq) ? "--three_prime_clip_r2 10" : ''
                        )
                    )
                )
            )
        ),
    ].join(' ').trim()
}

def bismarkArgs(s: SampleMeta, p: BismarkParams) -> BismarkArgs {
    return record(
        align: toolArgs('bismark_align', s, p.args, bismarkAlignArgs(s, p)),
        deduplicate: toolArgs('bismark_deduplicate', s, p.args, ''),
        methylation_extractor: toolArgs('bismark_methylationextractor', s, p.args, bismarkMethylationExtractorArgs(s, p)),
        coverage2cytosine: toolArgs('bismark_coverage2cytosine', s, p.args, p.nomeseq ? "--nome-seq" : ""),
        report: toolArgs('bismark_report', s, p.args, '')
    )
}

def bismarkAlignArgs(s: SampleMeta, p: BismarkParams) -> String {
    // Combined-index alignment is incompatible with --local_alignment, so gated off there
    def non_directional = p.single_cell || p.non_directional || p.zymo
    def use_combined = p.aligner.startsWith('bismark') && p.combined_index && !p.local_alignment
    return [
        (p.aligner == 'bismark_hisat') ? ' --hisat2' : ' --bowtie2',
        (p.aligner == 'bismark_hisat' && p.known_splices) ? " --known-splicesite-infile <(hisat2_extract_splice_sites.py ${p.known_splices})" : '',
        p.pbat ? ' --pbat' : '',
        non_directional ? ' --non_directional' : '',
        use_combined ? ' --combined_index' : '',
        (use_combined && non_directional) ? ' --combined_index_sequential' : '',
        p.unmapped ? ' --unmapped' : '',
        p.relax_mismatches ? " --score_min L,0,-${p.num_mismatches}" : '',
        p.local_alignment ? " --local" : '',
        !s.single_end && p.minins ? " --minins ${p.minins}" : '',
        s.single_end ? '' : (
            p.maxins ? " --maxins ${p.maxins}" : (
                p.em_seq ? " --maxins 1000" : ''
            )
        )
    ].join(' ').trim()
}

def bismarkMethylationExtractorArgs(s: SampleMeta, p: BismarkParams) -> String {
    return [
        p.comprehensive   ? ' --comprehensive' : '',
        p.meth_cutoff     ? " --cutoff ${p.meth_cutoff}" : '',
        p.nomeseq         ? '--CX' : '',
        p.ignore_r1 > 0   ? "--ignore ${p.ignore_r1}" : '',
        p.ignore_3prime_r1 > 0   ? "--ignore_3prime ${p.ignore_3prime_r1}" : '',
        s.single_end ? '' : (p.no_overlap           ? ' --no_overlap'                         : '--include_overlap'),
        s.single_end ? '' : (p.ignore_r2        > 0 ? "--ignore_r2 ${p.ignore_r2}"       : ""),
        s.single_end ? '' : (p.ignore_3prime_r2 > 0 ? "--ignore_3prime_r2 ${p.ignore_3prime_r2}": "")
    ].join(' ').trim()
}

def bismarkGenomePreparationArgs(p: BismarkParams) -> String {
    return runArgs('bismark_genomepreparation', p.args, [
        (p.aligner == 'bismark_hisat') ? ' --hisat2' : ' --bowtie2',
        p.slamseq ? ' --slam' : '',
        // Combined-index build is skipped for --local_alignment (combined rejects --local)
        (p.combined_index && !p.local_alignment && p.aligner.startsWith('bismark')) ? ' --combined_genome' : ''
    ].join(' ').trim())
}

def bwamethArgs(s: SampleMeta, args: Map<String,String>) -> BwamethArgs {
    return record(
        align: toolArgs('bwameth_align', s, args, ''),
        fq2bammeth: toolArgs('parabricks_fq2bammeth', s, args, '--low-memory'),
        markduplicates: toolArgs('picard_markduplicates', s, args, markduplicatesArgs())
    )
}

def bwamemArgs(s: SampleMeta, args: Map<String,String>) -> BwamemArgs {
    return record(
        align: toolArgs('bwa_mem', s, args, ''),
        fq2bam: toolArgs('parabricks_fq2bam', s, args, ''),
        addorreplacereadgroups: toolArgs('picard_addorreplacereadgroups', s, args, "--RGID 1 --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM sample1"),
        markduplicates: toolArgs('picard_markduplicates', s, args, markduplicatesArgs())
    )
}

def markduplicatesArgs() -> String {
    return "--ASSUME_SORTED true --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT --PROGRAM_RECORD_ID 'null' --TMP_DIR tmp"
}

def methyldackelArgs(s: SampleMeta, p: MethyldackelParams) -> MethyldackelArgs {
    def extract = [
        p.all_contexts ? ' --CHG --CHH' : '',
        p.merge_context ? ' --mergeContext' : '',
        p.ignore_flags ? " --ignoreFlags" : '',
        p.methyl_kit ? " --methylKit" : '',
        p.min_depth > 0 ? " --minDepth ${p.min_depth}" : ''
    ].join(" ").trim()
    def mbias = [
        p.all_contexts ? ' --CHG --CHH' : '',
        p.ignore_flags ? " --ignoreFlags" : ''
    ].join(" ").trim()
    return record(
        extract: toolArgs('methyldackel_extract', s, p.args, extract),
        mbias: toolArgs('methyldackel_mbias', s, p.args, mbias)
    )
}

def targetedSequencingArgs(s: SampleMeta, args: Map<String,String>) -> TargetedSequencingArgs {
    return record(
        intersect_cov: toolArgs('bedtools_intersect_cov', s, args, ''),
        collecthsmetrics: toolArgs('picard_collecthsmetrics', s, args, "--MINIMUM_MAPPING_QUALITY 20 --COVERAGE_CAP 1000  --NEAR_DISTANCE 500")
    )
}

def qualimapArgs(p: QualimapParams) -> String {
    return [
        p.genome?.startsWith('GRCh') ? '-gd HUMAN' : '',
        p.genome?.startsWith('GRCm') ? '-gd MOUSE' : ''
    ].join(" ").trim()
}

def methuratorArgs(p: MethuratorParams) -> String {
    return [
        p.methurator_compute_ci ? ' --compute_ci' : '',
        p.rrbs ? ' --rrbs' : '',
        p.methurator_minimum_coverage ? " --minimum-coverage ${p.methurator_minimum_coverage}" : "",
        p.methurator_t_max ? " --t-max ${p.methurator_t_max}" : ""
    ].join(" ").trim()
}

def multiqcArgs(p: MultiqcParams) -> String {
    return runArgs('multiqc', p.args, p.multiqc_title ? "--title \"${p.multiqc_title}\"" : '')
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
    args: Map<String,String>
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
    args: Map<String,String>
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
    args: Map<String,String>
    multiqc_title: String?
}
