#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
             B S - S E Q   M E T H Y L A T I O N  :  B W A - M E T H
========================================================================================
 Methylation (BS-Seq) Analysis Pipeline using bwa-meth. Started November 2016.
 #### Homepage / Documentation
 https://github.com/SciLifeLab/NGI-MethylSeq
 #### Authors
 Phil Ewels <phil.ewels@scilifelab.se>
----------------------------------------------------------------------------------------
*/


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1

// Configurable variables
params.project = false
params.genome = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.fasta_index = params.genome ? params.genomes[ params.genome ].fasta_index ?: false : false
params.bwa_meth_index = params.genome ? params.genomes[ params.genome ].bwa_meth ?: false : false
params.saveReference = true
params.reads = "data/*_R{1,2}.fastq.gz"
params.outdir = './results'
params.notrim = false
params.nodedup = false
params.allcontexts = false
params.mindepth = 0
params.ignoreFlags = false

// Validate inputs
if( params.bwa_meth_index ){
    bwa_meth_index = file("${params.bwa_meth_index}.bwameth.c2t.bwt")
    bwa_meth_indices = Channel.fromPath( "${params.bwa_meth_index}*" ).toList()
    if( !bwa_meth_index.exists() ) exit 1, "bwa-meth index not found: ${params.bwa_meth_index}"
}
if( params.fasta_index ){
    fasta_index = file(params.fasta_index)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta_index}"
}
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
} else {
    exit 1, "No reference Fasta file specified! Please use --fasta"
}

params.rrbs = false
params.pbat = false
params.single_cell = false
params.epignome = false
params.accel = false
params.zymo = false
params.cegx = false
if(params.pbat){
    params.clip_r1 = 6
    params.clip_r2 = 9
    params.three_prime_clip_r1 = 6
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

def single

log.info "==================================================="
log.info " NGI-MethylSeq : Bisulfite-Seq BWA-Meth v${version}"
log.info "==================================================="
log.info "Reads          : ${params.reads}"
log.info "Genome         : ${params.genome}"
log.info "Bismark Index  : ${params.bismark_index}"
log.info "Current home   : $HOME"
log.info "Current user   : $USER"
log.info "Current path   : $PWD"
log.info "Script dir     : $baseDir"
log.info "Working dir    : $workDir"
log.info "Output dir     : ${params.outdir}"
log.info "---------------------------------------------------"
if(params.rrbs){        log.info "RRBS Mode      : On" }
log.info "Deduplication  : ${params.nodedup ? 'No' : 'Yes'}"
log.info "PileOMeth      : C Contexts - ${params.allcontexts ? 'All (CpG, CHG, CHH)' : 'CpG only'}"
log.info "PileOMeth      : Minimum Depth - ${params.mindepth}"
if(params.ignoreFlags){ log.info "PileOMeth:     : Ignoring SAM Flags" }
log.info "---------------------------------------------------"
if(params.notrim){      log.info "Trimming Step  : Skipped" }
if(params.pbat){        log.info "Trim Profile   : PBAT" }
if(params.single_cell){ log.info "Trim Profile   : Single Cell" }
if(params.epignome){    log.info "Trim Profile   : Epignome" }
if(params.accel){       log.info "Trim Profile   : Accel" }
if(params.cegx){        log.info "Trim Profile   : CEGX" }
if(params.clip_r1 > 0)  log.info "Trim R1        : ${params.clip_r1}"
if(params.clip_r2 > 0)  log.info "Trim R2        : ${params.clip_r2}"
if(params.three_prime_clip_r1 > 0) log.info "Trim 3' R1     : ${params.three_prime_clip_r1}"
if(params.three_prime_clip_r2 > 0) log.info "Trim 3' R2     : ${params.three_prime_clip_r2}"
log.info "---------------------------------------------------"
log.info "Config Profile : ${workflow.profile}"
if(params.project) log.info "UPPMAX Project : ${params.project}"
log.info "==================================================="

// Validate inputs
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: -1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_fastqc; read_files_trimming }



/*
 * PREPROCESSING - Build bwa-mem index
 */
if(!params.bwa_meth_index){
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
if(!params.fasta_index){
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
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
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
    trimgalore_results = []
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy'
        
        input:
        set val(name), file(reads) from read_files_trimming
        
        output:
        set val(name), file('*fq.gz') into trimmed_reads
        file '*trimming_report.txt' into trimgalore_results
        
        script:
        single = reads instanceof Path
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        rrbs = params.rrbs ? "--rrbs" : ''
        if (single) {
            """
            trim_galore --gzip $rrbs $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --gzip $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}

/*
 * STEP 3 - align with bwa-mem
 */
process bwamem_align {
    tag "$name"
    publishDir "${params.outdir}/bwa-mem_alignments", mode: 'copy'
    
    input:
    set val(name), file(reads) from trimmed_reads
    file index from bwa_meth_index.first()
    file bwa_meth_indices from bwa_meth_indices.first()
    
    output:
    file '*.bam' into bam_aligned, bam_flagstat
    
    script:
    fasta = index.toString() - '.bwameth.c2t.bwt'
    prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    set -o pipefail   # Capture exit codes from bwa-meth
    bwameth.py \\
        --threads ${task.cpus} \\
        --reference $fasta \\
        $reads | samtools view -bS - > ${prefix}.bam
    """
}

/*
 * STEP 4.1 - samtools flagstat on samples
 */
process samtools_flagstat {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/bwa-mem_alignments", mode: 'copy'
    
    input:
    file bam from bam_flagstat
    
    output:
    file "${bam.baseName}_flagstat.txt" into flagstat_results
    file "${bam.baseName}_stats.txt" into samtools_stats_results
    
    script:
    """
    samtools flagstat $bam > ${bam.baseName}_flagstat.txt
    samtools stats $bam > ${bam.baseName}_stats.txt
    """
}
/*
 * STEP 4.2 - sort and index alignments
 */
process samtools_sort {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/bwa-mem_alignments_sorted", mode: 'copy'
    
    executor 'local'
    
    input:
    file bam from bam_aligned
    
    output:
    file "${bam.baseName}.sorted.bam" into bam_sorted, bam_for_index
    
    script:
    """
    samtools sort \\
        $bam
        -m ${task.memory.toBytes() / task.cpus} \\
        -@ ${task.cpus} \\
        > ${bam.baseName}.sorted.bam
    """
}
/*
 * STEP 4.3 - sort and index alignments
 */
process samtools_index {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/bwa-mem_alignments_sorted", mode: 'copy'
    
    input:
    file bam from bam_for_index
    
    output:
    file "${bam}.bai" into bam_index
    
    script:
    """
    samtools index $bam
    """
}


/*
 * STEP 5 - Mark duplicates
 */
process markDuplicates {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/bwa-mem_markDuplicates", mode: 'copy'

    input:
    file bam from bam_sorted

    output:
    file "${bam.baseName}.markDups.bam" into bam_md, bam_md_qualimap
    file "${bam.baseName}.markDups_metrics.txt" into picard_results

    script:
    """
    java -Xmx2g -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${bam.baseName}.markDups.bam \\
        METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT

    # Print version number to standard out
    echo "File name: $bam Picard version "\$(java -Xmx2g -jar \$PICARD_HOME/picard.jar  MarkDuplicates --version 2>&1)
    """
}


/*
 * STEP 6 - extract methylation with PileOMeth
 */
process pileOMeth {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/PileOMeth", mode: 'copy'
    
    input:
    file bam from bam_md
    file fasta from fasta
    file fasta_index from fasta
    
    output:
    file '*' into pileometh_results
    
    script:
    allcontexts = params.allcontexts ? '--CHG --CHH' : ''
    mindepth = params.mindepth > 0 ? "--minDepth ${params.mindepth}" : ''
    ignoreFlags = params.ignoreFlags ? "--ignoreFlags" : ''
    """
    PileOMeth extract $allcontexts $ignoreFlags $mindepth $fasta $bam
    PileOMeth mbias $allcontexts $ignoreFlags $fasta $bam ${bam.baseName}
    """
}

/*
 * STEP 7 - Qualimap
 */
process qualimap {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/Qualimap", mode: 'copy'
    
    input:
    file bam from bam_md_qualimap
    
    output:
    file '${bam.baseName}_qualimap' into qualimap_results
    
    script:
    gcref = params.genome == 'GRCh37' ? '-gd HUMAN' : ''
    gcref = params.genome == 'GRCm38' ? '-gd MOUSE' : ''
    """
    samtools sort $bam -o ${bam.baseName}.sorted.bam
    qualimap bamqc $gcref \\
        -bam ${bam.baseName}.sorted.bam \\
        -outdir ${bam.baseName}_qualimap \\
        --skip-duplicated \\
        --collect-overlap-pairs \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}
    """
}

/*
 * STEP 8 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    
    input:
    file ('fastqc/*') from fastqc_results.flatten().toList()
    file ('trimgalore/*') from trimgalore_results.flatten().toList()
    file ('samtools/*') from flagstat_results.flatten().toList()
    file ('samtools/*') from samtools_stats_results.flatten().toList()
    file ('picard/*') from picard_results.flatten().toList()
    file ('pileometh/*') from pileometh_results.flatten().toList()
    file ('qualimap/*') from qualimap_results.flatten().toList()
    
    output:
    file '*multiqc_report.html'
    file '*multiqc_data'
    
    script:
    """
    multiqc -f .
    """
}
