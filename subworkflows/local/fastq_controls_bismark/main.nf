//
// Align the first reads of each sample to a separate spike-in control FASTA with Bismark, and call every context
//

include { GUNZIP                       } from '../../../modules/nf-core/gunzip/main'
include { BISMARK_GENOMEPREPARATION    } from '../../../modules/nf-core/bismark/genomepreparation/main'
include { BWAMEM3_INDEX                } from '../../../modules/nf-core/bwamem3/index/main'
include { BISMARK_ALIGN                } from '../../../modules/nf-core/bismark/align/main'
include { BISMARK_METHYLATIONEXTRACTOR } from '../../../modules/nf-core/bismark/methylationextractor/main'
include { BISMARK_COVERAGE2CYTOSINE    } from '../../../modules/nf-core/bismark/coverage2cytosine/main'

workflow FASTQ_CONTROLS_BISMARK {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    control_fasta // channel: [ val(meta), path(fasta) ], one FASTA, plain or gzipped
    bwamem3 // boolean: Bismark's 5-Base/TAPS mode, which aligns with bwa-mem3

    main:
    ch_fasta = control_fasta.branch { _meta, fasta ->
        gzipped: fasta.name.endsWith('.gz')
        plain: true
    }
    GUNZIP(ch_fasta.gzipped)
    ch_fasta = GUNZIP.out.gunzip
        .mix(ch_fasta.plain)
        .map { meta, fasta ->
            checkControlFasta(fasta)
            [meta, fasta]
        }
        .first()

    // 5-Base/TAPS: bwa-mem3 aligns, and the extractor and coverage2cytosine need no Bismark index
    if (bwamem3) {
        BWAMEM3_INDEX(ch_fasta)
        ch_index = BWAMEM3_INDEX.out.index
        ch_bismark_index = channel.value([[:], []])
    }
    else {
        BISMARK_GENOMEPREPARATION(ch_fasta)
        ch_index = BISMARK_GENOMEPREPARATION.out.index
        ch_bismark_index = ch_index
    }

    ch_inputs = reads
        .combine(ch_fasta)
        .combine(ch_index)
        .multiMap { meta, fastq, meta_fasta, fasta, meta_index, index ->
            // meta.controls selects the control options in the shared Bismark configs
            reads: [meta + [controls: true], fastq]
            fasta: [meta_fasta, fasta]
            bismark_index: bwamem3 ? [[:], []] : [meta_index, index]
            bwamem3_index: bwamem3 ? [meta_index, index] : [[:], []]
        }

    // No deduplication: the controls are covered so deeply that most duplicates are not PCR copies
    BISMARK_ALIGN(
        ch_inputs.reads,
        ch_inputs.fasta,
        ch_inputs.bismark_index,
        ch_inputs.bwamem3_index,
    )

    BISMARK_METHYLATIONEXTRACTOR(
        BISMARK_ALIGN.out.bam,
        ch_bismark_index,
    )

    BISMARK_COVERAGE2CYTOSINE(
        BISMARK_METHYLATIONEXTRACTOR.out.coverage,
        ch_fasta,
        ch_bismark_index,
    )

    emit:
    cx_report        = BISMARK_COVERAGE2CYTOSINE.out.report // channel: [ val(meta), path(<id>.controls.CX_report.txt.gz) ]
    align_report     = BISMARK_ALIGN.out.report // channel: [ val(meta), path(controls.*_report.txt) ]
    splitting_report = BISMARK_METHYLATIONEXTRACTOR.out.report // channel: [ val(meta), path(controls.*_splitting_report.txt) ]
}

// The control names become MultiQC columns and are matched against the reference contigs
def checkControlFasta(fasta) {
    def ids = []
    def empty = []
    def bases = 0
    fasta.eachLine { line ->
        if (line.startsWith('>')) {
            if (ids && !bases) {
                empty << ids.last()
            }
            // Bismark names a sequence by the first word after '>'
            ids << (line.substring(1).tokenize() ?: [''])[0]
            if (line.length() > 1 && Character.isWhitespace(line.charAt(1))) {
                ids[-1] = line
            }
            bases = 0
        }
        else {
            bases += line.trim().length()
        }
    }
    if (ids && !bases) {
        empty << ids.last()
    }
    if (!ids) {
        error("--control_fasta ${fasta.name} holds no sequences.")
    }
    if (empty) {
        error("--control_fasta: these records have no sequence: ${empty.join(', ')}")
    }
    def invalid = ids.findAll { id -> !(id ==~ /[A-Za-z0-9._:|+@~-]+/) }
    if (invalid) {
        error("--control_fasta: sequence names may only hold letters, digits and ._:|+@~-, got: ${invalid.collect { id -> "'${id}'" }.join(', ')}")
    }
    def duplicated = ids.countBy { id -> id }.findAll { _id, count -> count > 1 }.keySet()
    if (duplicated) {
        error("--control_fasta names these sequences more than once: ${duplicated.join(', ')}")
    }
    if (fasta.size() > 10_000_000) {
        log.warn("--control_fasta ${fasta.name} is ${Math.round(fasta.size() / 1e6)} MB: it is meant for small spike-in sequences such as lambda and pUC19.")
    }
}
