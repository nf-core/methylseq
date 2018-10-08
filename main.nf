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
if (params.aligner != 'bismark' && params.aligner != 'bwameth'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'bismark', 'bwameth'"
}
if( params.bismark_index && params.aligner == 'bismark' ){
    bismark_index = Channel
        .fromPath(params.bismark_index)
        .ifEmpty { exit 1, "Bismark index not found: ${params.bismark_index}" }
}
else if( params.bwa_meth_index && params.aligner == 'bwameth' ){
    bwa_meth_indices = Channel
        .fromPath( "${params.bwa_meth_index}*" )
        .ifEmpty { exit 1, "bwa-meth index not found: ${params.bwa_meth_index}" }
}
else if( params.fasta_index && params.aligner == 'bwameth' ){
    fasta_index = file(params.fasta_index)
    if( !fasta_index.exists() ) exit 1, "Fasta index file not found: ${params.fasta_index}"
}
else if( !params.fasta ) {
    exit 1, "No reference genome index specified!"
}
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
else if( params.aligner == 'bwameth') {
    exit 1, "No Fasta reference specified! This is required by MethylDackel."
}
multiqc_config = file(params.multiqc_config)

// Validate inputs
if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax_devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
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

if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_files_fastqc; read_files_trimming }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_files_fastqc; read_files_trimming }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { read_files_fastqc; read_files_trimming }
}

log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/methylseq : Bisulfite-Seq Best Practice v${params.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/methylseq'
summary['Pipeline Version'] = params.version
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Aligner']        = params.aligner
summary['Data Type']      = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']         = params.genome
if(params.bismark_index) summary['Bismark Index'] = params.bismark_index
if(params.bwa_meth_index) summary['BWA-Meth Index'] = "${params.bwa_meth_index}*"
else if(params.fasta)    summary['Fasta Ref'] = params.fasta
if(params.rrbs) summary['RRBS Mode'] = 'On'
if(params.relaxMismatches) summary['Mismatch Func'] = "L,0,-${params.numMismatches} (Bismark default = L,0,-0.2)"
if(params.notrim)       summary['Trimming Step'] = 'Skipped'
if(params.pbat)         summary['Trim Profile'] = 'PBAT'
if(params.single_cell)  summary['Trim Profile'] = 'Single Cell'
if(params.epignome)     summary['Trim Profile'] = 'TruSeq (EpiGnome)'
if(params.accel)        summary['Trim Profile'] = 'Accel-NGS (Swift)'
if(params.zymo)         summary['Trim Profile'] = 'Zymo Pico-Methyl'
if(params.cegx)         summary['Trim Profile'] = 'CEGX'
summary['Trim R1'] = params.clip_r1
summary['Trim R2'] = params.clip_r2
summary["Trim 3' R1"] = params.three_prime_clip_r1
summary["Trim 3' R2"] = params.three_prime_clip_r2
summary['Deduplication']  = params.nodedup || params.rrbs ? 'No' : 'Yes'
summary['Directional Mode'] = params.single_cell || params.zymo || params.non_directional ? 'No' : 'Yes'
summary['All C Contexts'] = params.comprehensive ? 'Yes' : 'No'
if(params.mindepth) summary['Minimum Depth'] = params.mindepth
if(params.ignoreFlags) summary['MethylDackel'] = 'Ignoring SAM Flags'
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
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = System.getProperty("user.name")
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.project) summary['UPPMAX Project'] = params.project
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "========================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
  if( ! nextflow.version.matches(">= $params.nf_required_version") ){
    throw GroovyException('Nextflow version too old')
  }
} catch (all) {
  log.error "====================================================\n" +
            "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}
// Show a big error message if we're running on the base config and an uppmax cluster
if( workflow.profile == 'standard'){
    if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
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
if(!params.bismark_index && params.fasta && params.aligner == 'bismark'){
    process makeBismarkIndex {
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "BismarkIndex" into bismark_index

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
if(!params.bwa_meth_index && params.fasta && params.aligner == 'bwameth'){
    process makeBwaMemIndex {
        tag fasta
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "${fasta}.bwameth.c2t.bwt" into bwa_meth_index
        file "${fasta}*" into bwa_meth_indices

        script:
        """
        bwameth.py index $fasta
        """
    }
}

/*
 * PREPROCESSING - Index Fasta file
 */
if(!params.fasta_index && params.fasta && params.aligner == 'bwameth'){
    process makeFastaIndex {
        tag fasta
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta

        output:
        file "${fasta}.fai" into fasta_index

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
    set val(name), file(reads) from read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 2 - Trim Galore!
 */
if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = Channel.from(false)
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
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        rrbs = params.rrbs ? "--rrbs" : ''
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
 * STEP 3.1 - align with Bismark
 */
if(params.aligner == 'bismark'){
    process bismark_align {
        tag "$name"
        publishDir "${params.outdir}/bismark_alignments", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".fq.gz") > 0) "unmapped/$filename"
                else if (filename.indexOf(".bam") == -1) "logs/$filename"
                else params.saveAlignedIntermediates || params.nodedup || params.rrbs ? filename : null
            }

        input:
        set val(name), file(reads) from trimmed_reads
        file index from bismark_index.collect()

        output:
        file "*.bam" into bam, bam_2
        file "*report.txt" into bismark_align_log_1, bismark_align_log_2, bismark_align_log_3
        if(params.unmapped){ file "*.fq.gz" into bismark_unmapped }

        script:
        pbat = params.pbat ? "--pbat" : ''
        non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
        unmapped = params.unmapped ? "--unmapped" : ''
        mismatches = params.relaxMismatches ? "--score_min L,0,-${params.numMismatches}" : ''
        multicore = ''
        if (task.cpus){
            // Numbers based on recommendation by Felix for a typical mouse genome
            if(params.single_cell || params.zymo || params.non_directional){
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
            if(ccore > 1){
              multicore = "--multicore $ccore"
            }
        }
        if (params.singleEnd) {
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
    if (params.nodedup || params.rrbs) {
        bam.into { bam_dedup; bam_dedup_qualimap }
        bismark_dedup_log_1 = Channel.from(false)
        bismark_dedup_log_2 = Channel.from(false)
        bismark_dedup_log_3 = Channel.from(false)
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

            script:
            if (params.singleEnd) {
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

        script:
        comprehensive = params.comprehensive ? '--comprehensive --merge_non_CpG' : ''
        multicore = ''
        if (task.cpus){
            // Numbers based on Bismark docs
            ccore = ((task.cpus as int) / 10) as int
            if(ccore > 1){
              multicore = "--multicore $ccore"
            }
        }
        buffer = ''
        if (task.memory){
            mbuffer = (task.memory as nextflow.util.MemoryUnit) - 2.GB
            // only set if we have more than 6GB available
            if(mbuffer.compareTo(4.GB) == 1){
              buffer = "--buffer_size ${mbuffer.toGiga()}G"
            }
        }
        if (params.singleEnd) {
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

        script:
        name = bismark_align_log_1.toString() - ~/(_R1)?(_trimmed|_val_1).+$/
        """
        bismark2report \\
            --alignment_report $bismark_align_log_1 \\
            --dedup_report $bismark_dedup_log_1 \\
            --splitting_report $bismark_splitting_report_1 \\
            --mbias_report $bismark_mbias_1
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

        script:
        """
        bismark2summary
        """
    }
} // End of bismark processing block
else {
    bismark_align_log_3 = Channel.from(false)
    bismark_dedup_log_3 = Channel.from(false)
    bismark_splitting_report_3 = Channel.from(false)
    bismark_mbias_3 = Channel.from(false)
    bismark_reports_results = Channel.from(false)
    bismark_summary_results = Channel.from(false)
}


/*
 * Process with bwa-mem and assorted tools
 */
if(params.aligner == 'bwameth'){
    process bwamem_align {
        tag "$name"
        publishDir "${params.outdir}/bwa-mem_alignments", mode: 'copy',
            saveAs: { fn -> params.saveAlignedIntermediates ? fn : null }

        input:
        set val(name), file(reads) from trimmed_reads
        file bwa_meth_indices from bwa_meth_indices.collect()

        output:
        file '*.bam' into bam_aligned

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
        tag "${bam.baseName}"
        publishDir "${params.outdir}/bwa-mem_alignments", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".txt") > 0) "logs/$filename"
                else if (params.saveAlignedIntermediates || params.nodedup || params.rrbs) filename
                else null
            }

        input:
        file bam from bam_aligned

        output:
        file "${bam.baseName}.sorted.bam" into bam_sorted
        file "${bam.baseName}.sorted.bam.bai" into bam_index
        file "${bam.baseName}_flagstat.txt" into flagstat_results
        file "${bam.baseName}_stats.txt" into samtools_stats_results

        script:
        """
        samtools sort \\
            $bam \\
            -m ${task.memory.toBytes() / task.cpus} \\
            -@ ${task.cpus} \\
            > ${bam.baseName}.sorted.bam
        samtools index ${bam.baseName}.sorted.bam
        samtools flagstat ${bam.baseName}.sorted.bam > ${bam.baseName}_flagstat.txt
        samtools stats ${bam.baseName}.sorted.bam > ${bam.baseName}_stats.txt
        """
    }

    /*
     * STEP 5 - Mark duplicates
     */
    if (params.nodedup || params.rrbs) {
        bam_sorted.into { bam_md; bam_dedup_qualimap }
        bam_index.set { bam_md_bai }
        picard_results = Channel.from(false)
    } else {
        process markDuplicates {
            tag "${bam.baseName}"
            publishDir "${params.outdir}/bwa-mem_markDuplicates", mode: 'copy',
                saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

            input:
            file bam from bam_sorted

            output:
            file "${bam.baseName}.markDups.bam" into bam_md, bam_dedup_qualimap
            file "${bam.baseName}.markDups.bam.bai" into bam_md_bai
            file "${bam.baseName}.markDups_metrics.txt" into picard_results

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
        tag "${bam.baseName}"
        publishDir "${params.outdir}/MethylDackel", mode: 'copy'

        input:
        file bam from bam_md
        file bam_index from bam_md_bai
        file fasta from fasta
        file fasta_index from fasta_index

        output:
        file "${bam.baseName}*" into methyldackel_results

        script:
        allcontexts = params.comprehensive ? '--CHG --CHH' : ''
        mindepth = params.mindepth > 0 ? "--minDepth ${params.mindepth}" : ''
        ignoreFlags = params.ignoreFlags ? "--ignoreFlags" : ''
        """
        MethylDackel extract $allcontexts $ignoreFlags $mindepth $fasta $bam
        MethylDackel mbias $allcontexts $ignoreFlags $fasta $bam ${bam.baseName} --txt > ${bam.baseName}_methyldackel.txt
        """
    }

} // end of bwa-meth if block
else {
    flagstat_results = Channel.from(false)
    samtools_stats_results = Channel.from(false)
    picard_results = Channel.from(false)
    methyldackel_results = Channel.from(false)
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
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo "$params.version" &> v_ngi_methylseq.txt
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
    qualimap --version &> v_qualimap.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}



/*
 * STEP 9 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('trimgalore/*') from trimgalore_results.collect().ifEmpty([])
    file ('bismark/*') from bismark_align_log_3.collect().ifEmpty([])
    file ('bismark/*') from bismark_dedup_log_3.collect().ifEmpty([])
    file ('bismark/*') from bismark_splitting_report_3.collect().ifEmpty([])
    file ('bismark/*') from bismark_mbias_3.collect().ifEmpty([])
    file ('bismark/*') from bismark_reports_results.collect().ifEmpty([])
    file ('bismark/*') from bismark_summary_results.collect().ifEmpty([])
    file ('samtools/*') from flagstat_results.flatten().collect().ifEmpty([])
    file ('samtools/*') from samtools_stats_results.flatten().collect().ifEmpty([])
    file ('picard/*') from picard_results.flatten().collect().ifEmpty([])
    file ('methyldackel/*') from methyldackel_results.flatten().collect().ifEmpty([])
    file ('qualimap/*') from qualimap_results.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect().ifEmpty([])

    output:
    file "*_report.html" into multiqc_report
    file "*_data"
    file '.command.err' into multiqc_stderr

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
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
    if(!workflow.success){
      subject = "[nf-core/methylseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
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
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

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

    if(!workflow.success){
        if( workflow.profile == 'standard'){
            if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
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
