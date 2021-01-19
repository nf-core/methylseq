/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include {
    SAMPLESHEET_CHECK;
    get_samplesheet_paths;
    get_genome_paths} from '../process/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    
    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .set { sample }
        
    reads = sample.map { get_samplesheet_paths(it) }
    genome = sample.map { get_genome_paths(it, params.genomes) }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
    genome // channel: [ val(meta), [ fasta ] ]
}
