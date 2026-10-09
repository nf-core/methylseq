// Gather per-sample MethylDackel bedGraphs into one HDF5 file
process COMBINED_HDF5 {
    label 'process_medium'

    container '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0'

    input:
    path bedgraph_files
    val  bedgraph_meta

    output:
    path "bedGraph.id_file.csv", emit: id_file
    path "merged.bedgraph.hdf5", emit: merged_bedgraph_hdf5
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    make_combined_hdf5.py \\
        --input_files '${bedgraph_files}' \\
        --meta '${bedgraph_meta}' \\
        --id_file bedGraph.id_file.csv \\
        --hdf5 merged.bedgraph.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        h5py: \$(python -c 'import h5py; print(h5py.__version__)' 2>/dev/null)
    END_VERSIONS
    """
}
