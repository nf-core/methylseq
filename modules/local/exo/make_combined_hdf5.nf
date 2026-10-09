

process COMBINED_HDF5 {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0' :
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.20.0' }"
	
	input:
    path(bedGraph_files)
	val(bedGraph_meta)

	output:
	path "bedGraph.id_file.csv"       , emit: id_file
	path "merged.bedgraph.hdf5", emit: merged_bedgraph_hdf5
	path "versions.yml"      , emit: versions

	script:
	"""

	make_combined_hdf5.py \\
		--input_files '$bedGraph_files' \\
		--meta '$bedGraph_meta' \\
		--id_file bedGraph.id_file.csv \\
		--hdf5 merged.bedgraph.hdf5
	

	cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        h5py: \$(python -c 'import h5py; print(h5py.__version__)')
    END_VERSIONS
	"""


}