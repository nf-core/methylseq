#!/usr/bin/env nextflow
/*
========================================================================================
             B S - S E Q   M E T H Y L A T I O N   B E S T - P R A C T I C E
========================================================================================
 New Methylation (BS-Seq) Best Practice Analysis Pipeline. Started June 2016.
 #### Homepage / Documentation
 https://github.com/nf-core/methylseq
 #### Authors
 Phil Ewels <phil.ewels@scilifelab.se>
----------------------------------------------------------------------------------------
*/


/*
 * SET UP CONFIGURATION VARIABLES
 */
params.name = false
params.project = false
params.clusterOptions = false
params.email = false
params.plaintext_email = false
params.genome = false
params.bismark_index = params.genome ? params.genomes[ params.genome ].bismark ?: false : false
params.bwa_meth_index = params.genome ? params.genomes[ params.genome ].bwa_meth ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.fasta_index = params.genome ? params.genomes[ params.genome ].fasta_index ?: false : false

// Validate inputs
assert params.aligner == 'bwameth' || params.aligner == 'bismark' : "Invalid aligner option: ${params.aligner}. Valid options: 'bismark', 'bwameth'"

Channel
    .fromPath("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)
    .into { ch_wherearemyfiles_for_trimgalore; ch_wherearemyfiles_for_alignment }

if( params.aligner == 'bismark' ){
    assert params.bismark_index || params.fasta : "No reference genome index or fasta file specified"
    ch_wherearemyfiles_for_alignment.set { ch_wherearemyfiles_for_bismark_align }

    if( params.bismark_index ){
        Channel
            .fromPath(params.bismark_index, checkIfExists: true)
            .ifEmpty { exit 1, "Bismark index file not found: ${params.bismark_index}" }
            .set { ch_bismark_index_for_bismark_align }
    }
    else if( params.fasta ){
        Channel
            .fromPath(params.fasta, checkIfExists: true)
            .ifEmpty { exit 1, "fasta file not found : ${params.fasta}" }
            .set { ch_fasta_for_makeBismarkIndex }
    }
}
else if( params.aligner == 'bwameth' ){
    assert params.fasta : "No Fasta reference specified! This is required by MethylDackel."
    ch_wherearemyfiles_for_alignment.into { ch_wherearemyfiles_for_bwamem_align; ch_wherearemyfiles_for_samtools_sort_index_flagstat }

    Channel
        .fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "fasta file not found : ${params.fasta}" }
        .into { ch_fasta_for_makeBwaMemIndex; ch_fasta_for_makeFastaIndex; ch_fasta_for_methyldackel }

    if( params.bwa_meth_index ){
        Channel
            .fromPath("${params.bwa_meth_index}*", checkIfExists: true)
            .ifEmpty { exit 1, "bwa-meth index file(s) not found: ${params.bwa_meth_index}" }
            .set { ch_bwa_meth_indices_for_bwamem_align }
        ch_fasta_for_makeBwaMemIndex.close()
    }

    if( params.fasta_index ){
        Channel
            .fromPath(params.fasta_index, checkIfExists: true)
            .ifEmpty { exit 1, "fasta index file not found: ${params.fasta_index}" }
            .set { ch_fasta_index_for_methyldackel }
        ch_fasta_for_makeFastaIndex.close()
    }
}

Channel
    .fromPath(params.multiqc_config, checkIfExists: true)
    .ifEmpty { exit 1, "multiqc config file not found: ${params.multiqc_config}" }
    .set { ch_config_for_multiqc }


if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax_devel' ){
    if( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Library prep presets
params.rrbs = false
params.pbat = false
params.single_cell = false
params.epignome = false
params.accel = false
params.zymo = false
params.cegx = false
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
} else {
    params.clip_r1 = 0
    params.clip_r2 = 0
    params.three_prime_clip_r1 = 0
    params.three_prime_clip_r2 = 0
}

/*
 * Create a channel for input read files
 */

if( params.readPaths ){
    if( params.singleEnd ){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { ch_read_files_for_fastqc; ch_read_files_for_trim_galore }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { ch_read_files_for_fastqc; ch_read_files_for_trim_galore }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { ch_read_files_for_fastqc; ch_read_files_for_trim_galore }
}

log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/methylseq : Bisulfite-Seq Best Practice v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/methylseq'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Aligner']        = params.aligner
summary['Data Type']      = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']         = params.genome
if( params.bismark_index ) summary['Bismark Index'] = params.bismark_index
if( params.bwa_meth_index ) summary['BWA-Meth Index'] = "${params.bwa_meth_index}*"
if( params.fasta )    summary['Fasta Ref'] = params.fasta
if( params.fasta_index )    summary['Fasta Index'] = params.fasta_index
if( params.rrbs ) summary['RRBS Mode'] = 'On'
if( params.relaxMismatches ) summary['Mismatch Func'] = "L,0,-${params.numMismatches} (Bismark default = L,0,-0.2)"
if( params.notrim )       summary['Trimming Step'] = 'Skipped'
if( params.pbat )         summary['Trim Profile'] = 'PBAT'
if( params.single_cell )  summary['Trim Profile'] = 'Single Cell'
if( params.epignome )     summary['Trim Profile'] = 'TruSeq (EpiGnome)'
if( params.accel )        summary['Trim Profile'] = 'Accel-NGS (Swift)'
if( params.zymo )         summary['Trim Profile'] = 'Zymo Pico-Methyl'
if( params.cegx )         summary['Trim Profile'] = 'CEGX'
summary['Trim R1'] = params.clip_r1
summary['Trim R2'] = params.clip_r2
summary["Trim 3' R1"] = params.three_prime_clip_r1
summary["Trim 3' R2"] = params.three_prime_clip_r2
summary['Deduplication']  = params.nodedup || params.rrbs ? 'No' : 'Yes'
summary['Directional Mode'] = params.single_cell || params.zymo || params.non_directional ? 'No' : 'Yes'
summary['All C Contexts'] = params.comprehensive ? 'Yes' : 'No'
if( params.mindepth ) summary['Minimum Depth'] = params.mindepth
if( params.ignoreFlags ) summary['MethylDackel'] = 'Ignoring SAM Flags'
if( params.methylKit ) summary['MethylDackel'] = 'Producing methylKit output'
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
summary['Save Unmapped']  = params.unmapped ? 'Yes' : 'No'
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if( workflow.containerEngine ) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = System.getProperty("user.name")
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if( params.project ) summary['UPPMAX Project'] = params.project
if( params.email ) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "========================================="

// Show a big error message if we're running on the base config and an uppmax cluster
if( workflow.profile == 'standard' ){
    if( "hostname".execute().text.contains('.uppmax.uu.se') ) {
        log.error "====================================================\n" +
                  "  WARNING! You are running with the default 'standard'\n" +
                  "  pipeline config profile, which runs on the head node\n" +
                  "  and assumes all software is on the PATH.\n" +
                  "  ALL JOBS ARE RUNNING LOCALLY and stuff will probably break.\n" +
                  "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                  "============================================================"
    }
}


/*
 * PREPROCESSING - Build Bismark index
 */
if( !params.bismark_index && params.aligner == 'bismark' ){
    process makeBismarkIndex {
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_makeBismarkIndex

        output:
        file "BismarkIndex" into ch_bismark_index_for_bismark_align

        script:
        """
        mkdir BismarkIndex
        cp $fasta BismarkIndex/
        bismark_genome_preparation BismarkIndex
        """
    }
}

/*
 * PREPROCESSING - Build bwa-mem index
 */
if( !params.bwa_meth_index && params.aligner == 'bwameth' ){
    process makeBwaMemIndex {
        tag "$fasta"
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_makeBwaMemIndex

        output:
        file "${fasta}*" into ch_bwa_meth_indices_for_bwamem_align

        script:
        """
        bwameth.py index $fasta
        """
    }
}

/*
 * PREPROCESSING - Index Fasta file
 */
if( !params.fasta_index && params.aligner == 'bwameth' ){
    process makeFastaIndex {
        tag "$fasta"
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_makeFastaIndex

        output:
        file "${fasta}.fai" into ch_fasta_index_for_methyldackel

        script:
        """
        samtools faidx $fasta
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
    set val(name), file(reads) from ch_read_files_for_fastqc

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc_results_for_multiqc

    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 2 - Trim Galore!
 */
if( params.notrim ){
    ch_trimmed_reads_for_alignment = ch_read_files_for_trim_galore
    ch_trim_galore_results_for_multiqc = Channel.from(false)
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if( filename.indexOf("_fastqc") > 0 ) "FastQC/$filename"
                else if( filename.indexOf("trimming_report.txt" ) > 0) "logs/$filename"
                else if( !params.saveTrimmed && filename == "where_are_my_files.txt" ) filename
                else if( params.saveTrimmed && filename != "where_are_my_files.txt" ) filename
                else null
            }

        input:
        set val(name), file(reads) from ch_read_files_for_trim_galore
        file wherearemyfiles from ch_wherearemyfiles_for_trimgalore.collect()

        output:
        set val(name), file('*fq.gz') into ch_trimmed_reads_for_alignment
        file "*trimming_report.txt" into ch_trim_galore_results_for_multiqc
        file "*_fastqc.{zip,html}"
        file "where_are_my_files.txt"

        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        rrbs = params.rrbs ? "--rrbs" : ''
        if( params.singleEnd ) {
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
 * STEP 3.1 - align with Bismark
 */
if( params.aligner == 'bismark' ){
    process bismark_align {
        tag "$name"
        publishDir "${params.outdir}/bismark_alignments", mode: 'copy',
            saveAs: {filename ->
                if( filename.indexOf(".fq.gz") > 0 ) "unmapped/$filename"
                else if( filename.indexOf("report.txt") > 0 ) "logs/$filename"
                else if( (!params.saveAlignedIntermediates && !params.nodedup && !params.rrbs).every() && filename == "where_are_my_files.txt" ) filename
                else if( (params.saveAlignedIntermediates || params.nodedup || params.rrbs).any() && filename != "where_are_my_files.txt" ) filename
                else null
            }

        input:
        set val(name), file(reads) from ch_trimmed_reads_for_alignment
        file index from ch_bismark_index_for_bismark_align.collect()
        file wherearemyfiles from ch_wherearemyfiles_for_bismark_align.collect()

        output:
        set val(name), file("*.bam") into ch_bam_for_bismark_deduplicate, ch_bam_for_bismark_summary, ch_bam_for_preseq
        set val(name), file("*report.txt") into ch_bismark_align_log_for_bismark_report, ch_bismark_align_log_for_bismark_summary, ch_bismark_align_log_for_multiqc
        file "*.fq.gz" optional true
        file "where_are_my_files.txt"

        script:
        pbat = params.pbat ? "--pbat" : ''
        non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
        unmapped = params.unmapped ? "--unmapped" : ''
        mismatches = params.relaxMismatches ? "--score_min L,0,-${params.numMismatches}" : ''
        multicore = ''
        if( task.cpus ){
            // Numbers based on recommendation by Felix for a typical mouse genome
            if( params.single_cell || params.zymo || params.non_directional ){
                cpu_per_multicore = 5
                mem_per_multicore = (18.GB).toBytes()
            } else {
                cpu_per_multicore = 3
                mem_per_multicore = (13.GB).toBytes()
            }
            // How many multicore splits can we afford with the cpus we have?
            ccore = ((task.cpus as int) / cpu_per_multicore) as int
            // Check that we have enough memory, assuming 13GB memory per instance (typical for mouse alignment)
            try {
                tmem = (task.memory as nextflow.util.MemoryUnit).toBytes()
                mcore = (tmem / mem_per_multicore) as int
                ccore = Math.min(ccore, mcore)
            } catch (all) {
                log.debug "Not able to define bismark align multicore based on available memory"
            }
            if( ccore > 1 ){
              multicore = "--multicore $ccore"
            }
        }
        if( params.singleEnd ) {
            """
            bismark \\
                --bam $pbat $non_directional $unmapped $mismatches $multicore \\
                --genome $index \\
                $reads
            """
        } else {
            """
            bismark \\
                --bam $pbat $non_directional $unmapped $mismatches $multicore \\
                --genome $index \\
                -1 ${reads[0]} \\
                -2 ${reads[1]}
            """
        }
    }

    /*
     * STEP 4 - Bismark deduplicate
     */
    if( params.nodedup || params.rrbs ) {
        ch_bam_for_bismark_deduplicate.into { ch_bam_dedup_for_bismark_methXtract; ch_bam_dedup_for_qualimap }
        ch_bismark_dedup_log_for_bismark_report = Channel.from(false)
        ch_bismark_dedup_log_for_bismark_summary = Channel.from(false)
        ch_bismark_dedup_log_for_multiqc  = Channel.from(false)
    } else {
        process bismark_deduplicate {
            tag "$name"
            publishDir "${params.outdir}/bismark_deduplicated", mode: 'copy',
                saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

            input:
            set val(name), file(bam) from ch_bam_for_bismark_deduplicate

            output:
            set val(name), file("*.deduplicated.bam") into ch_bam_dedup_for_bismark_methXtract, ch_bam_dedup_for_qualimap
            set val(name), file("*.deduplication_report.txt") into ch_bismark_dedup_log_for_bismark_report, ch_bismark_dedup_log_for_bismark_summary, ch_bismark_dedup_log_for_multiqc

            script:
            if( params.singleEnd ) {
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
        tag "$name"
        publishDir "${params.outdir}/bismark_methylation_calls", mode: 'copy',
            saveAs: {filename ->
                if( filename.indexOf("splitting_report.txt" ) > 0 ) "logs/$filename"
                else if( filename.indexOf("M-bias" ) > 0) "m-bias/$filename"
                else if( filename.indexOf(".cov" ) > 0 ) "methylation_coverage/$filename"
                else if( filename.indexOf("bedGraph" ) > 0 ) "bedGraph/$filename"
                else "methylation_calls/$filename"
            }

        input:
        set val(name), file(bam) from ch_bam_dedup_for_bismark_methXtract

        output:
        set val(name), file("*splitting_report.txt") into ch_bismark_splitting_report_for_bismark_report, ch_bismark_splitting_report_for_bismark_summary, ch_bismark_splitting_report_for_multiqc
        set val(name), file("*.M-bias.txt") into ch_bismark_mbias_for_bismark_report, ch_bismark_mbias_for_bismark_summary, ch_bismark_mbias_for_multiqc
        file '*.{png,gz}'

        script:
        comprehensive = params.comprehensive ? '--comprehensive --merge_non_CpG' : ''
        multicore = ''
        if( task.cpus ){
            // Numbers based on Bismark docs
            ccore = ((task.cpus as int) / 10) as int
            if( ccore > 1 ){
              multicore = "--multicore $ccore"
            }
        }
        buffer = ''
        if( task.memory ){
            mbuffer = (task.memory as nextflow.util.MemoryUnit) - 2.GB
            // only set if we have more than 6GB available
            if( mbuffer.compareTo(4.GB) == 1 ){
              buffer = "--buffer_size ${mbuffer.toGiga()}G"
            }
        }
        if(params.singleEnd) {
            """
            bismark_methylation_extractor $comprehensive \\
                $multicore $buffer \\
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
                $multicore $buffer \\
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

    ch_bismark_align_log_for_bismark_report
     .join(ch_bismark_dedup_log_for_bismark_report)
     .join(ch_bismark_splitting_report_for_bismark_report)
     .join(ch_bismark_mbias_for_bismark_report)
     .set{ ch_bismark_logs_for_bismark_report }


    /*
     * STEP 6 - Bismark Sample Report
     */
    process bismark_report {
        tag "$name"
        publishDir "${params.outdir}/bismark_reports", mode: 'copy'

        input:
        set val(name), file(align_log), file(dedup_log), file(splitting_report), file(mbias) from ch_bismark_logs_for_bismark_report

        output:
        file '*{html,txt}' into ch_bismark_reports_results_for_multiqc

        script:
        """
        bismark2report \\
            --alignment_report $align_log \\
            --dedup_report $dedup_log \\
            --splitting_report $splitting_report \\
            --mbias_report $mbias
        """
    }

    /*
     * STEP 7 - Bismark Summary Report
     */
    process bismark_summary {
        publishDir "${params.outdir}/bismark_summary", mode: 'copy'

        input:
        file ('*') from ch_bam_for_bismark_summary.collect()
        file ('*') from ch_bismark_align_log_for_bismark_summary.collect()
        file ('*') from ch_bismark_dedup_log_for_bismark_summary.collect()
        file ('*') from ch_bismark_splitting_report_for_bismark_summary.collect()
        file ('*') from ch_bismark_mbias_for_bismark_summary.collect()

        output:
        file '*{html,txt}' into ch_bismark_summary_results_for_multiqc

        script:
        """
        bismark2summary
        """
    }
} // End of bismark processing block
else {
    ch_bismark_align_log_for_multiqc = Channel.from(false)
    ch_bismark_dedup_log_for_multiqc = Channel.from(false)
    ch_bismark_splitting_report_for_multiqc = Channel.from(false)
    ch_bismark_mbias_for_multiqc = Channel.from(false)
    ch_bismark_reports_results_for_multiqc = Channel.from(false)
    ch_bismark_summary_results_for_multiqc = Channel.from(false)
}


/*
 * Process with bwa-mem and assorted tools
 */
if( params.aligner == 'bwameth' ){
    process bwamem_align {
        tag "$name"
        publishDir "${params.outdir}/bwa-mem_alignments", mode: 'copy',
            saveAs: {filename ->
                if( !params.saveAlignedIntermediates && filename == "where_are_my_files.txt" ) filename
                else if( params.saveAlignedIntermediates && filename != "where_are_my_files.txt" ) filename
                else null
            }

        input:
        set val(name), file(reads) from ch_trimmed_reads_for_alignment
        file bwa_meth_indices from ch_bwa_meth_indices_for_bwamem_align.collect()
        file wherearemyfiles from ch_wherearemyfiles_for_bwamem_align.collect()

        output:
        set val(name), file('*.bam') into ch_bam_for_samtools_sort_index_flagstat, ch_bam_for_preseq
        file "where_are_my_files.txt"

        script:
        fasta = bwa_meth_indices[0].toString() - '.bwameth' - '.c2t' - '.amb' - '.ann' - '.bwt' - '.pac' - '.sa'
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        bwameth.py \\
            --threads ${task.cpus} \\
            --reference $fasta \\
            $reads | samtools view -bS - > ${prefix}.bam
        """
    }


    /*
     * STEP 4.- samtools flagstat on samples
     */
    process samtools_sort_index_flagstat {
        tag "$name"
        publishDir "${params.outdir}/bwa-mem_alignments", mode: 'copy',
            saveAs: {filename ->
                if(filename.indexOf("report.txt") > 0) "logs/$filename"
                else if( (!params.saveAlignedIntermediates && !params.nodedup && !params.rrbs).every() && filename == "where_are_my_files.txt") filename
                else if( (params.saveAlignedIntermediates || params.nodedup || params.rrbs).any() && filename != "where_are_my_files.txt") filename
                else null
            }

        input:
        set val(name), file(bam) from ch_bam_for_samtools_sort_index_flagstat
        file wherearemyfiles from ch_wherearemyfiles_for_samtools_sort_index_flagstat.collect()

        output:
        set val(name), file("${bam.baseName}.sorted.bam") into ch_bam_sorted_for_markDuplicates
        file "${bam.baseName}.sorted.bam.bai" into ch_bam_index
        file "${bam.baseName}_flagstat_report.txt" into ch_flagstat_results_for_multiqc
        file "${bam.baseName}_stats_report.txt" into ch_samtools_stats_results_for_multiqc
        file "where_are_my_files.txt"

        script:
        def avail_mem = task.memory ? ((task.memory.toBytes() - 6000000000) / task.cpus) : false
        def sort_mem = avail_mem && avail_mem > 2000000000 ? "-m $avail_mem" : ''
        """
        samtools sort $bam \\
            -@ ${task.cpus} $sort_mem \\
            -o ${bam.baseName}.sorted.bam
        samtools index ${bam.baseName}.sorted.bam
        samtools flagstat ${bam.baseName}.sorted.bam > ${bam.baseName}_flagstat_report.txt
        samtools stats ${bam.baseName}.sorted.bam > ${bam.baseName}_stats_report.txt
        """
    }

    /*
     * STEP 5 - Mark duplicates
     */
    if( params.nodedup || params.rrbs ) {
        ch_bam_sorted_for_markDuplicates.into { ch_bam_dedup_for_methyldackel; ch_bam_dedup_for_qualimap }
        ch_bam_index.set { ch_bam_index_for_methyldackel }
        ch_markDups_results_for_multiqc = Channel.from(false)
    } else {
        process markDuplicates {
            tag "$name"
            publishDir "${params.outdir}/bwa-mem_markDuplicates", mode: 'copy',
                saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

            input:
            set val(name), file(bam) from ch_bam_sorted_for_markDuplicates

            output:
            set val(name), file("${bam.baseName}.markDups.bam") into ch_bam_dedup_for_methyldackel, ch_bam_dedup_for_qualimap
            file "${bam.baseName}.markDups.bam.bai" into ch_bam_index_for_methyldackel //ToDo check if this correctly overrides the original channel
            file "${bam.baseName}.markDups_metrics.txt" into ch_markDups_results_for_multiqc

            script:
            if( !task.memory ){
                log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
                avail_mem = 3
            } else {
                avail_mem = task.memory.toGiga()
            }
            """
            picard -Xmx${avail_mem}g MarkDuplicates \\
                INPUT=$bam \\
                OUTPUT=${bam.baseName}.markDups.bam \\
                METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
                REMOVE_DUPLICATES=false \\
                ASSUME_SORTED=true \\
                PROGRAM_RECORD_ID='null' \\
                VALIDATION_STRINGENCY=LENIENT
            samtools index ${bam.baseName}.markDups.bam
            """
        }
    }

    /*
     * STEP 6 - extract methylation with MethylDackel
     */
    process methyldackel {
        tag "$name"
        publishDir "${params.outdir}/MethylDackel", mode: 'copy'

        input:
        set val(name), file(bam) from ch_bam_dedup_for_methyldackel
        file bam_index from ch_bam_index_for_methyldackel
        file fasta from ch_fasta_for_methyldackel
        file fasta_index from ch_fasta_index_for_methyldackel

        output:
        file "${bam.baseName}*" into ch_methyldackel_results_for_multiqc

        script:
        allcontexts = params.comprehensive ? '--CHG --CHH' : ''
        mindepth = params.mindepth > 0 ? "--minDepth ${params.mindepth}" : ''
        ignoreFlags = params.ignoreFlags ? "--ignoreFlags" : ''
        methylKit = params.methylKit ? "--methylKit" : ''
        """
        MethylDackel extract $allcontexts $ignoreFlags $methylKit $mindepth $fasta $bam
        MethylDackel mbias $allcontexts $ignoreFlags $fasta $bam ${bam.baseName} --txt > ${bam.baseName}_methyldackel.txt
        """
    }

} // end of bwa-meth if block
else {
    ch_flagstat_results_for_multiqc = Channel.from(false)
    ch_samtools_stats_results_for_multiqc = Channel.from(false)
    ch_markDups_results_for_multiqc = Channel.from(false)
    ch_methyldackel_results_for_multiqc = Channel.from(false)
}


/*
 * STEP 8 - Qualimap
 */
process qualimap {
    tag "$name"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    input:
    set val(name), file(bam) from ch_bam_dedup_for_qualimap

    output:
    file "${bam.baseName}_qualimap" into ch_qualimap_results_for_multiqc

    script:
    gcref = params.genome == 'GRCh37' ? '-gd HUMAN' : ''
    gcref = params.genome == 'GRCm38' ? '-gd MOUSE' : ''
    def avail_mem = task.memory ? ((task.memory.toBytes() - 6000000000) / task.cpus) : false
    def sort_mem = avail_mem && avail_mem > 2000000000 ? "-m $avail_mem" : ''
    """
    samtools sort $bam \\
        -@ ${task.cpus} $sort_mem \\
        -o ${bam.baseName}.sorted.bam
    qualimap bamqc $gcref \\
        -bam ${bam.baseName}.sorted.bam \\
        -outdir ${bam.baseName}_qualimap \\
        --collect-overlap-pairs \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}
    """
}

/*
 * STEP 9 - preseq
 */
process preseq {
    tag "$name"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    input:
    set val(name), file(bam) from ch_bam_for_preseq

    output:
    file "${bam.baseName}.ccurve.txt" into preseq_results

    script:
    def avail_mem = task.memory ? ((task.memory.toBytes() - 6000000000) / task.cpus) : false
    def sort_mem = avail_mem && avail_mem > 2000000000 ? "-m $avail_mem" : ''
    """
    samtools sort $bam \\
        -@ ${task.cpus} $sort_mem \\
        -o ${bam.baseName}.sorted.bam
    preseq lc_extrap -v -B ${bam.baseName}.sorted.bam -o ${bam.baseName}.ccurve.txt
    """
}

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml_for_multiqc

    script:
    """
    echo "$workflow.manifest.version" &> v_ngi_methylseq.txt
    echo "$workflow.nextflow.version" &> v_nextflow.txt
    bismark_genome_preparation --version &> v_bismark_genome_preparation.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    bismark --version &> v_bismark.txt
    deduplicate_bismark --version &> v_deduplicate_bismark.txt
    bismark_methylation_extractor --version &> v_bismark_methylation_extractor.txt
    bismark2report --version &> v_bismark2report.txt
    bismark2summary --version &> v_bismark2summary.txt
    samtools --version &> v_samtools.txt
    bwa &> v_bwa.txt 2>&1 || true
    bwameth.py --version &> v_bwameth.txt
    picard MarkDuplicates --version &> v_picard_markdups.txt 2>&1 || true
    MethylDackel --version &> v_methyldackel.txt
    qualimap --version &> v_qualimap.txt || true
    preseq &> v_preseq.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}



/*
 * STEP 10 - MultiQC
 */
process multiqc {
    tag "${params.outdir}/MultiQC/$ofilename"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_config_for_multiqc
    file ('fastqc/*') from ch_fastqc_results_for_multiqc.collect().ifEmpty([])
    file ('trimgalore/*') from ch_trim_galore_results_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_align_log_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_dedup_log_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_splitting_report_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_mbias_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_reports_results_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_summary_results_for_multiqc.collect().ifEmpty([])
    file ('samtools/*') from ch_flagstat_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('samtools/*') from ch_samtools_stats_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('picard/*') from ch_markDups_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('methyldackel/*') from ch_methyldackel_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('qualimap/*') from ch_qualimap_results_for_multiqc.collect().ifEmpty([])
    file ('preseq/*') from preseq_results.collect().ifEmpty([])
    file ('software_versions/*') from ch_software_versions_yaml_for_multiqc.collect().ifEmpty([])

    output:
    file "*_report.html" into ch_multiqc_report
    file "*_data"
    file '.command.err'

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    if(custom_runName){
      rfilename = "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
      ofilename = rfilename+'.html'
    } else {
      rfilename = ''
      ofilename = 'multiqc_report.html'
    }
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/methylseq] Successful: $workflow.runName"
    if( !workflow.success ){
      subject = "[nf-core/methylseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
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
    email_fields['summary']['Container'] = workflow.container
    if( workflow.repository ) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if( workflow.commitId ) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if( workflow.revision ) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if( workflow.success ) {
            mqc_report = ch_multiqc_report.getVal()
            if( mqc_report.getClass() == ArrayList ){
                log.warn "[nfcore/methylseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
                }
        }
    } catch (all) {
        log.warn "[nfcore/methylseq] Could not attach MultiQC report to summary email"
    }

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
    if( params.email ) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/methylseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/methylseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Switch the embedded MIME images with base64 encoded src
    ngimethylseqlogo = new File("$baseDir/assets/methylseq_logo.png").bytes.encodeBase64().toString()
    scilifelablogo = new File("$baseDir/assets/SciLifeLab_logo.png").bytes.encodeBase64().toString()
    ngilogo = new File("$baseDir/assets/NGI_logo.png").bytes.encodeBase64().toString()
    email_html = email_html.replaceAll(~/cid:ngimethylseqlogo/, "data:image/png;base64,$ngimethylseqlogo")
    email_html = email_html.replaceAll(~/cid:scilifelablogo/, "data:image/png;base64,$scilifelablogo")
    email_html = email_html.replaceAll(~/cid:ngilogo/, "data:image/png;base64,$ngilogo")

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/methylseq] Pipeline Complete"

    if( !workflow.success ){
        if( workflow.profile == 'standard' ){
            if( "hostname".execute().text.contains('.uppmax.uu.se') ) {
                log.error "====================================================\n" +
                        "  WARNING! You are running with the default 'standard'\n" +
                        "  pipeline config profile, which runs on the head node\n" +
                        "  and assumes all software is on the PATH.\n" +
                        "  This is probably why everything broke.\n" +
                        "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                        "============================================================"
            }
        }
    }

}
