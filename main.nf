#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/methylseq
========================================================================================
 nf-core/methylseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/methylseq
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/methylseq --input '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --aligner [str]                 Alignment tool to use (default: bismark)
                                      Available: bismark, bismark_hisat, bwameth
      --input [file]                  Path to input data (must be surrounded with quotes)
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more

    Options:
     --genome [str]                     Name of iGenomes reference
     --single_end [bool]                Specifies that the input is single end reads
     --comprehensive [bool]             Output information for all cytosine contexts
     --cytosine_report [bool]           Output stranded cytosine report during Bismark's bismark_methylation_extractor step.
     --ignore_flags [bool]              Run MethylDackel with the flag to ignore SAM flags.
     --meth_cutoff [int]                Specify a minimum read coverage to report a methylation call during Bismark's bismark_methylation_extractor step.
     --min_depth [int]                  Specify a minimum read coverage for MethylDackel to report a methylation call.
     --methyl_kit [bool]                Run MethylDackel with the --methyl_kit flag to produce files suitable for use with the methylKit R package.
     --skip_deduplication [bool]        Skip deduplication step after alignment. This is turned on automatically if --rrbs is specified
     --non_directional [bool]           Run alignment against all four possible strands
     --save_align_intermeds [bool]      Save aligned intermediates to results directory
     --save_trimmed [bool]              Save trimmed reads to results directory
     --unmapped [bool]                  Save unmapped reads to fastq files
     --relax_mismatches [bool]          Turn on to relax stringency for alignment (set allowed penalty with --num_mismatches)
     --num_mismatches [float]           0.6 will allow a penalty of bp * -0.6 - for 100bp reads (bismark default is 0.2)
     --known_splices [file]             Supply a .gtf file containing known splice sites (bismark_hisat only)
     --slamseq [bool]                   Run bismark in SLAM-seq mode
     --local_alignment [bool]           Allow soft-clipping of reads (potentially useful for single-cell experiments)
     --minins [int]                     Bismark: The minimum insert size for valid paired-end alignments.
     --maxins [int]                     Bismark: The maximum insert size for valid paired-end alignments.
     --bismark_align_cpu_per_multicore [int] Specify how many CPUs are required per --multicore for bismark align (default = 3)
     --bismark_align_mem_per_multicore [str] Specify how much memory is required per --multicore for bismark align (default = 13.GB)

    References                          If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta [file]                    Path to fasta reference
      --fasta_index [path]              Path to Fasta Index
      --bismark_index [path]            Path to Bismark index
      --bwa_meth_index [path]           Path to bwameth index
      --save_reference [bool]           Save reference(s) to results directory

    Trimming options:
     --skip_trimming [bool]             Skip read trimming
     --clip_r1 [int]                    Trim the specified number of bases from the 5' end of read 1 (or single-end reads).
     --clip_r2 [int]                    Trim the specified number of bases from the 5' end of read 2 (paired-end only).
     --three_prime_clip_r1 [int]        Trim the specified number of bases from the 3' end of read 1 AFTER adapter/quality trimming
     --three_prime_clip_r2 [int]        Trim the specified number of bases from the 3' end of read 2 AFTER adapter/quality trimming
     --rrbs [bool]                      Turn on if dealing with MspI digested material.

    Trimming presets:
     --pbat [bool]
     --single_cell [bool]
     --epignome [bool]
     --accell [bool]
     --zymo [bool]
     --cegx [bool]
     --em_seq [bool]

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                  The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]                 The AWS Region for your AWS Batch job to run on
      --awscli [str]                    Path to the AWS CLI tool

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
assert params.aligner == 'bwameth' || params.aligner == 'bismark' || params.aligner == 'bismark_hisat' : "Invalid aligner option: ${params.aligner}. Valid options: 'bismark', 'bwameth', 'bismark_hisat'"

/*
 * SET UP CONFIGURATION VARIABLES
 */

// These params need to be set late, after the iGenomes config is loaded
params.bismark_index = params.genome ? params.genomes[ params.genome ].bismark ?: false : false
params.bwa_meth_index = params.genome ? params.genomes[ params.genome ].bwa_meth ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.fasta_index = params.genome ? params.genomes[ params.genome ].fasta_index ?: false : false

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

Channel
    .fromPath("$projectDir/assets/where_are_my_files.txt", checkIfExists: true)
    .into { ch_wherearemyfiles_for_trimgalore; ch_wherearemyfiles_for_alignment }

ch_splicesites_for_bismark_hisat_align = params.known_splices ? Channel.fromPath(params.known_splices, checkIfExists: true) : Channel.empty()

if( params.aligner =~ /bismark/ ){
    assert params.bismark_index || params.fasta : "No reference genome index or fasta file specified"
    ch_wherearemyfiles_for_alignment.set { ch_wherearemyfiles_for_bismark_align }

    if( params.bismark_index ){
        Channel
            .fromPath(params.bismark_index, checkIfExists: true)
            .ifEmpty { exit 1, "Bismark index file not found: ${params.bismark_index}" }
            .into { ch_bismark_index_for_bismark_align; ch_bismark_index_for_bismark_methXtract }
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

if( workflow.profile == 'uppmax' ){
    if( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Trimming / kit presets
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
bismark_minins = params.minins
bismark_maxins = params.maxins
if(params.pbat){
    clip_r1 = 9
    clip_r2 = 9
    three_prime_clip_r1 = 9
    three_prime_clip_r2 = 9
}
else if( params.single_cell ){
    clip_r1 = 6
    clip_r2 = 6
    three_prime_clip_r1 = 6
    three_prime_clip_r2 = 6
}
else if( params.epignome ){
    clip_r1 = 8
    clip_r2 = 8
    three_prime_clip_r1 = 8
    three_prime_clip_r2 = 8
}
else if( params.accel || params.zymo ){
    clip_r1 = 10
    clip_r2 = 15
    three_prime_clip_r1 = 10
    three_prime_clip_r2 = 10
}
else if( params.cegx ){
    clip_r1 = 6
    clip_r2 = 6
    three_prime_clip_r1 = 2
    three_prime_clip_r2 = 2
}
else if( params.em_seq ){
    bismark_maxins = 1000
    clip_r1 = 8
    clip_r2 = 8
    three_prime_clip_r1 = 8
    three_prime_clip_r2 = 8
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

/*
 * Create a channel for input read files
 */
if (params.input_paths) {
    if (params.single_end) {
        Channel
            .from(params.input_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.input_paths was empty - no input files supplied" }
            .into { ch_read_files_fastqc; ch_read_files_trimming }
    } else {
        Channel
            .from(params.input_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.input_paths was empty - no input files supplied" }
            .into { ch_read_files_fastqc; ch_read_files_trimming }
    }
} else {
    Channel
        .fromFilePairs(params.input, size: params.single_end ? 1 : 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .into { ch_read_files_fastqc; ch_read_files_trimming }
}

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']  = custom_runName ?: workflow.runName
summary['Input']     = params.input
summary['Aligner']   = params.aligner
summary['Data Type'] = params.single_end ? 'Single-End' : 'Paired-End'
if(params.known_splices)    summary['Spliced alignment'] =  'Yes'
if(params.slamseq)          summary['SLAM-seq'] = 'Yes'
if(params.local_alignment)  summary['Local alignment'] = 'Yes'
if(params.genome)           summary['Genome']    = params.genome
if(params.bismark_index)    summary['Bismark Index'] = params.bismark_index
if(params.bwa_meth_index)   summary['BWA-Meth Index'] = "${params.bwa_meth_index}*"
if(params.fasta)            summary['Fasta Ref'] = params.fasta
if(params.fasta_index)      summary['Fasta Index'] = params.fasta_index
if(params.rrbs)             summary['RRBS Mode'] = 'On'
if(params.relax_mismatches) summary['Mismatch Func'] = "L,0,-${params.num_mismatches} (Bismark default = L,0,-0.2)"
if(params.skip_trimming)    summary['Trimming Step'] = 'Skipped'
if(params.pbat)             summary['Trim Profile'] = 'PBAT'
if(params.single_cell)      summary['Trim Profile'] = 'Single Cell'
if(params.epignome)         summary['Trim Profile'] = 'TruSeq (EpiGnome)'
if(params.accel)            summary['Trim Profile'] = 'Accel-NGS (Swift)'
if(params.zymo)             summary['Trim Profile'] = 'Zymo Pico-Methyl'
if(params.cegx)             summary['Trim Profile'] = 'CEGX'
if(params.em_seq)           summary['Trim Profile'] = 'EM Seq'
summary['Trimming']         = "5'R1: $clip_r1 / 5'R2: $clip_r2 / 3'R1: $three_prime_clip_r1 / 3'R2: $three_prime_clip_r2"
summary['Deduplication']    = params.skip_deduplication || params.rrbs ? 'No' : 'Yes'
summary['Directional Mode'] = params.single_cell || params.zymo || params.non_directional ? 'No' : 'Yes'
summary['All C Contexts']   = params.comprehensive ? 'Yes' : 'No'
summary['Cytosine report']  = params.cytosine_report ? 'Yes' : 'No'
if(params.min_depth)        summary['Minimum Depth'] = params.min_depth
if(params.ignore_flags)     summary['MethylDackel'] = 'Ignoring SAM Flags'
if(params.methyl_kit)       summary['MethylDackel'] = 'Producing methyl_kit output'
save_intermeds = [];
if(params.save_reference)   save_intermeds.add('Reference genome build')
if(params.save_trimmed)     save_intermeds.add('Trimmed FastQ files')
if(params.unmapped)         save_intermeds.add('Unmapped reads')
if(params.save_align_intermeds) save_intermeds.add('Intermediate BAM files')
if(save_intermeds.size() > 0) summary['Save Intermediates'] = save_intermeds.join(', ')
if(params.minins)           summary['Bismark min insert size'] = bismark_minins
if(params.maxins || params.em_seq) summary['Bismark max insert size'] = bismark_maxins
if(params.bismark_align_cpu_per_multicore) summary['Bismark align CPUs per --multicore'] = params.bismark_align_cpu_per_multicore
if(params.bismark_align_mem_per_multicore) summary['Bismark align memory per --multicore'] = params.bismark_align_mem_per_multicore
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Pipeline dir']     = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
if(params.project) summary['Cluster Project'] = params.project
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-methylseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/methylseq Workflow Summary'
    section_href: 'https://github.com/nf-core/methylseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml_for_multiqc
    file "software_versions.csv"

    script:
    """
    echo "$workflow.manifest.version" &> v_pipeline.txt
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
    hisat2 --version &> v_hisat2.txt
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
 * PREPROCESSING - Build Bismark index
 */
if( !params.bismark_index && params.aligner =~ /bismark/ ){
    process makeBismarkIndex {
        publishDir path: { params.save_reference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        file fasta from ch_fasta_for_makeBismarkIndex

        output:
        file "BismarkIndex" into ch_bismark_index_for_bismark_align, ch_bismark_index_for_bismark_methXtract

        script:
        aligner = params.aligner == 'bismark_hisat' ? '--hisat2' : '--bowtie2'
        slam = params.slamseq ? '--slam' : ''
        """
        mkdir BismarkIndex
        cp $fasta BismarkIndex/
        bismark_genome_preparation $aligner $slam BismarkIndex
        """
    }
}

/*
 * PREPROCESSING - Build bwa-mem index
 */
if( !params.bwa_meth_index && params.aligner == 'bwameth' ){
    process makeBwaMemIndex {
        tag "$fasta"
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        file fasta from ch_fasta_for_makeBwaMemIndex

        output:
        file "${fasta}*" into ch_bwa_meth_indices_for_bwamem_align
        file fasta

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
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

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
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from ch_read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc_results_for_multiqc

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}

/*
 * STEP 2 - Trim Galore!
 */
if( params.skip_trimming ){
    ch_trimmed_reads_for_alignment = ch_read_files_trimming
    ch_trim_galore_results_for_multiqc = Channel.from(false)
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: params.publish_dir_mode,
            saveAs: {filename ->
                if( filename.indexOf("_fastqc") > 0 ) "FastQC/$filename"
                else if( filename.indexOf("trimming_report.txt" ) > 0) "logs/$filename"
                else if( !params.save_trimmed && filename == "where_are_my_files.txt" ) filename
                else if( params.save_trimmed && filename != "where_are_my_files.txt" ) filename
                else null
            }

        input:
        set val(name), file(reads) from ch_read_files_trimming
        file wherearemyfiles from ch_wherearemyfiles_for_trimgalore.collect()

        output:
        set val(name), file('*fq.gz') into ch_trimmed_reads_for_alignment
        file "*trimming_report.txt" into ch_trim_galore_results_for_multiqc
        file "*_fastqc.{zip,html}"
        file "where_are_my_files.txt"

        script:
        def c_r1 = clip_r1 > 0 ? "--clip_r1 $clip_r1" : ''
        def c_r2 = clip_r2 > 0 ? "--clip_r2 $clip_r2" : ''
        def tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 $three_prime_clip_r1" : ''
        def tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 $three_prime_clip_r2" : ''
        def rrbs = params.rrbs ? "--rrbs" : ''
        def cores = 1
        if(task.cpus){
            cores = (task.cpus as int) - 4
            if (params.single_end) cores = (task.cpus as int) - 3
            if (cores < 1) cores = 1
            if (cores > 4) cores = 4
        }
        if( params.single_end ) {
            """
            trim_galore --fastqc --gzip $reads \
              $rrbs $c_r1 $tpc_r1 --cores $cores
            """
        } else {
            """
            trim_galore --fastqc --gzip --paired $reads \
              $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 --cores $cores
            """
        }
    }
}

/*
 * STEP 3.1 - align with Bismark
 */
if( params.aligner =~ /bismark/ ){
    process bismark_align {
        tag "$name"
        publishDir "${params.outdir}/bismark_alignments", mode: params.publish_dir_mode,
            saveAs: {filename ->
                if( filename.indexOf(".fq.gz") > 0 ) "unmapped/$filename"
                else if( filename.indexOf("report.txt") > 0 ) "logs/$filename"
                else if( (!params.save_align_intermeds && !params.skip_deduplication && !params.rrbs).every() && filename == "where_are_my_files.txt" ) filename
                else if( (params.save_align_intermeds || params.skip_deduplication || params.rrbs).any() && filename != "where_are_my_files.txt" ) filename
                else null
            }

        input:
        set val(name), file(reads) from ch_trimmed_reads_for_alignment
        file index from ch_bismark_index_for_bismark_align.collect()
        file wherearemyfiles from ch_wherearemyfiles_for_bismark_align.collect()
        file knownsplices from ch_splicesites_for_bismark_hisat_align.collect().ifEmpty([])

        output:
        set val(name), file("*.bam") into ch_bam_for_bismark_deduplicate, ch_bam_for_bismark_summary, ch_bam_for_preseq
        set val(name), file("*report.txt") into ch_bismark_align_log_for_bismark_report, ch_bismark_align_log_for_bismark_summary, ch_bismark_align_log_for_multiqc
        file "*.fq.gz" optional true
        file "where_are_my_files.txt"

        script:
        // Paired-end or single end input files
        input = params.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"

        // Choice of read aligner
        aligner = params.aligner == "bismark_hisat" ? "--hisat2" : "--bowtie2"

        // Optional extra bismark parameters
        splicesites = params.aligner == "bismark_hisat" && knownsplices.name != 'null' ? "--known-splicesite-infile <(hisat2_extract_splice_sites.py ${knownsplices})" : ''
        pbat = params.pbat ? "--pbat" : ''
        non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
        unmapped = params.unmapped ? "--unmapped" : ''
        mismatches = params.relax_mismatches ? "--score_min L,0,-${params.num_mismatches}" : ''
        soft_clipping = params.local_alignment ? "--local" : ''
        minins = bismark_minins ? "--minins $bismark_minins" : ''
        maxins = bismark_maxins ? "--maxins $bismark_maxins" : ''

        // Try to assign sensible bismark memory units according to what the task was given
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
            // Check if the user has specified this and overwrite if so
            if(params.bismark_align_cpu_per_multicore) {
                cpu_per_multicore = (params.bismark_align_cpu_per_multicore as int)
            }
            if(params.bismark_align_mem_per_multicore) {
                mem_per_multicore = (params.bismark_align_mem_per_multicore as nextflow.util.MemoryUnit).toBytes()
            }
            // How many multicore splits can we afford with the cpus we have?
            ccore = ((task.cpus as int) / cpu_per_multicore) as int
            // Check that we have enough memory, assuming 13GB memory per instance (typical for mouse alignment)
            try {
                tmem = (task.memory as nextflow.util.MemoryUnit).toBytes()
                mcore = (tmem / mem_per_multicore) as int
                ccore = Math.min(ccore, mcore)
            } catch (all) {
                log.debug "Warning: Not able to define bismark align multicore based on available memory"
            }
            if( ccore > 1 ){
              multicore = "--multicore $ccore"
            }
        }

        // Main command
        """
        bismark $input \\
            $aligner \\
            --bam $pbat $non_directional $unmapped $mismatches $multicore $minins $maxins \\
            --genome $index \\
            $reads \\
            $soft_clipping \\
            $splicesites
        """
    }

    /*
     * STEP 4 - Bismark deduplicate
     */
    if( params.skip_deduplication || params.rrbs ) {
        ch_bam_for_bismark_deduplicate.into { ch_bam_dedup_for_bismark_methXtract; ch_bam_dedup_for_qualimap }
        ch_bismark_dedup_log_for_bismark_report = Channel.from(false)
        ch_bismark_dedup_log_for_bismark_summary = Channel.from(false)
        ch_bismark_dedup_log_for_multiqc  = Channel.from(false)
    } else {
        process bismark_deduplicate {
            tag "$name"
            publishDir "${params.outdir}/bismark_deduplicated", mode: params.publish_dir_mode,
                saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

            input:
            set val(name), file(bam) from ch_bam_for_bismark_deduplicate

            output:
            set val(name), file("*.deduplicated.bam") into ch_bam_dedup_for_bismark_methXtract, ch_bam_dedup_for_qualimap
            set val(name), file("*.deduplication_report.txt") into ch_bismark_dedup_log_for_bismark_report, ch_bismark_dedup_log_for_bismark_summary, ch_bismark_dedup_log_for_multiqc

            script:
            fq_type = params.single_end ? '-s' : '-p'
            """
            deduplicate_bismark $fq_type --bam $bam
            """
        }
    }

    /*
     * STEP 5 - Bismark methylation extraction
     */
    process bismark_methXtract {
        tag "$name"
        publishDir "${params.outdir}/bismark_methylation_calls", mode: params.publish_dir_mode,
            saveAs: {filename ->
                if( filename.indexOf("splitting_report.txt" ) > 0 ) "logs/$filename"
                else if( filename.indexOf("M-bias" ) > 0) "m-bias/$filename"
                else if( filename.indexOf(".cov" ) > 0 ) "methylation_coverage/$filename"
                else if( filename.indexOf("bedGraph" ) > 0 ) "bedGraph/$filename"
                else if( filename.indexOf("CpG_report" ) > 0 ) "stranded_CpG_report/$filename"
                else "methylation_calls/$filename"
            }

        input:
        set val(name), file(bam) from ch_bam_dedup_for_bismark_methXtract
        file index from ch_bismark_index_for_bismark_methXtract.collect()

        output:
        set val(name), file("*splitting_report.txt") into ch_bismark_splitting_report_for_bismark_report, ch_bismark_splitting_report_for_bismark_summary, ch_bismark_splitting_report_for_multiqc
        set val(name), file("*.M-bias.txt") into ch_bismark_mbias_for_bismark_report, ch_bismark_mbias_for_bismark_summary, ch_bismark_mbias_for_multiqc
        file '*.{png,gz}'

        script:
        comprehensive = params.comprehensive ? '--comprehensive --merge_non_CpG' : ''
        cytosine_report = params.cytosine_report ? "--cytosine_report --genome_folder ${index} " : ''
        meth_cutoff = params.meth_cutoff ? "--cutoff ${params.meth_cutoff}" : ''
        multicore = ''
        if( task.cpus ){
            // Numbers based on Bismark docs
            ccore = ((task.cpus as int) / 3) as int
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
        if(params.single_end) {
            """
            bismark_methylation_extractor $comprehensive $meth_cutoff \\
                $multicore $buffer $cytosine_report \\
                --bedGraph \\
                --counts \\
                --gzip \\
                -s \\
                --report \\
                $bam
            """
        } else {
            """
            bismark_methylation_extractor $comprehensive $meth_cutoff \\
                $multicore $buffer $cytosine_report \\
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
        publishDir "${params.outdir}/bismark_reports", mode: params.publish_dir_mode

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
        publishDir "${params.outdir}/bismark_summary", mode: params.publish_dir_mode

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
        publishDir "${params.outdir}/bwa-mem_alignments", mode: params.publish_dir_mode,
            saveAs: {filename ->
                if( !params.save_align_intermeds && filename == "where_are_my_files.txt" ) filename
                else if( params.save_align_intermeds && filename != "where_are_my_files.txt" ) filename
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
        publishDir "${params.outdir}/bwa-mem_alignments", mode: params.publish_dir_mode,
            saveAs: {filename ->
                if(filename.indexOf("report.txt") > 0) "logs/$filename"
                else if( (!params.save_align_intermeds && !params.skip_deduplication && !params.rrbs).every() && filename == "where_are_my_files.txt") filename
                else if( (params.save_align_intermeds || params.skip_deduplication || params.rrbs).any() && filename != "where_are_my_files.txt") filename
                else null
            }

        input:
        set val(name), file(bam) from ch_bam_for_samtools_sort_index_flagstat
        file wherearemyfiles from ch_wherearemyfiles_for_samtools_sort_index_flagstat.collect()

        output:
        set val(name), file("${bam.baseName}.sorted.bam") into ch_bam_sorted_for_markDuplicates
        set val(name), file("${bam.baseName}.sorted.bam.bai") into ch_bam_index
        file "${bam.baseName}_flagstat_report.txt" into ch_flagstat_results_for_multiqc
        file "${bam.baseName}_stats_report.txt" into ch_samtools_stats_results_for_multiqc
        file "where_are_my_files.txt"

        script:
        def avail_mem = task.memory ? ((task.memory.toGiga() - 6) / task.cpus).trunc() : false
        def sort_mem = avail_mem && avail_mem > 2 ? "-m ${avail_mem}G" : ''
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
    if( params.skip_deduplication || params.rrbs ) {
        ch_bam_sorted_for_markDuplicates.into { ch_bam_dedup_for_methyldackel; ch_bam_dedup_for_qualimap }
        ch_bam_index.set { ch_bam_index_for_methyldackel }
        ch_markDups_results_for_multiqc = Channel.from(false)
    } else {
        process markDuplicates {
            tag "$name"
            publishDir "${params.outdir}/bwa-mem_markDuplicates", mode: params.publish_dir_mode,
                saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

            input:
            set val(name), file(bam) from ch_bam_sorted_for_markDuplicates

            output:
            set val(name), file("${bam.baseName}.markDups.bam") into ch_bam_dedup_for_methyldackel, ch_bam_dedup_for_qualimap
            set val(name), file("${bam.baseName}.markDups.bam.bai") into ch_bam_index_for_methyldackel //ToDo check if this correctly overrides the original channel
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
        publishDir "${params.outdir}/MethylDackel", mode: params.publish_dir_mode

        input:
        set val(name),
            file(bam),
            file(bam_index),
            file(fasta),
            file(fasta_index) from ch_bam_dedup_for_methyldackel
            .join(ch_bam_index_for_methyldackel)
            .combine(ch_fasta_for_methyldackel)
            .combine(ch_fasta_index_for_methyldackel)


        output:
        file "${bam.baseName}*" into ch_methyldackel_results_for_multiqc

        script:
        all_contexts = params.comprehensive ? '--CHG --CHH' : ''
        min_depth = params.min_depth > 0 ? "--minDepth ${params.min_depth}" : ''
        ignore_flags = params.ignore_flags ? "--ignoreFlags" : ''
        methyl_kit = params.methyl_kit ? "--methylKit" : ''
        """
        MethylDackel extract $all_contexts $ignore_flags $methyl_kit $min_depth $fasta $bam
        MethylDackel mbias $all_contexts $ignore_flags $fasta $bam ${bam.baseName} --txt > ${bam.baseName}_methyldackel.txt
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
    publishDir "${params.outdir}/qualimap", mode: params.publish_dir_mode

    input:
    set val(name), file(bam) from ch_bam_dedup_for_qualimap

    output:
    file "${bam.baseName}_qualimap" into ch_qualimap_results_for_multiqc

    script:
    gcref = params.genome.toString().startsWith('GRCh') ? '-gd HUMAN' : ''
    gcref = params.genome.toString().startsWith('GRCm') ? '-gd MOUSE' : ''
    def avail_mem = task.memory ? ((task.memory.toGiga() - 6) / task.cpus).trunc() : false
    def sort_mem = avail_mem && avail_mem > 2 ? "-m ${avail_mem}G" : ''
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
    publishDir "${params.outdir}/preseq", mode: params.publish_dir_mode

    input:
    set val(name), file(bam) from ch_bam_for_preseq

    output:
    file "${bam.baseName}.ccurve.txt" into preseq_results

    script:
    def avail_mem = task.memory ? ((task.memory.toGiga() - 6) / task.cpus).trunc() : false
    def sort_mem = avail_mem && avail_mem > 2 ? "-m ${avail_mem}G" : ''
    """
    samtools sort $bam \\
        -@ ${task.cpus} $sort_mem \\
        -o ${bam.baseName}.sorted.bam
    preseq lc_extrap -v -B ${bam.baseName}.sorted.bam -o ${bam.baseName}.ccurve.txt
    """
}

/*
 * STEP 10 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
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
    file ('software_versions/*') from ch_software_versions_yaml_for_multiqc.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $custom_config_file . \\
        -m custom_content -m picard -m qualimap -m bismark -m samtools -m preseq -m cutadapt -m fastqc
    """
}

/*
 * STEP 11 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/methylseq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/methylseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
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
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/methylseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
                }
        }
    } catch (all) {
        log.warn "[nfcore/methylseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/methylseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/methylseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/methylseq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/methylseq]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/methylseq v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
