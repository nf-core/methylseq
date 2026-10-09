

process BASIC_READ_STATISTICS {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/basic_read_statistics:0.0.1' :
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/basic_read_statistics:0.0.1' }"
	
	input:
	path(bedGraph_hdf5)
    path(mbias_files)
	val(mbias_files_meta)

	output:
	path "merged_mbias.csv", emit: merged_mbias
	path "mbias.id_file.csv", emit: mbias_id_file
	path "tab_id_value_counts_df.csv", emit: tab_id_value_counts_df
	path "rnf_counts_id_df.csv", emit: rnf_counts_id_df
	path "tab_count_by_rep.csv", emit: tab_count_by_rep
	path "lt100_tab_count_by_rep.csv", emit: lt100_tab_count_by_rep
	path "metrics_table.csv", emit: metrics_table
	path "combined_coverage_distributions.png", emit: combined_coverage_distributions_plot
	path "read_count_occurence_frequencies.png", emit: read_count_occurence_frequencies_plot
	path "read_bias_in_methylated_OT_strand.png", emit: read_bias_in_methylated_OT_strand_plot
	path "read_bias_in_methylated_OB_strand.png", emit: read_bias_in_methylated_OB_strand_plot
	path "read_bias_in_unmethylated_OT_strand.png", emit: read_bias_in_unmethylated_OT_strand_plot
	path "read_bias_in_unmethylated_OB_strand.png", emit: read_bias_in_unmethylated_OB_strand_plot
	path "versions.yml"      , emit: versions

	script:
	"""

	basic_read_statistics.py \\
		--bedGraph_hdf5 $bedGraph_hdf5 \\
		--mbias_files '$mbias_files' \\
		--meta '$mbias_files_meta' \\
		--id_file mbias.id_file.csv \\
		--merged_mbias merged_mbias.csv 
		
	

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    basic_read_statistics: \$(echo \$(basic_read_statistics.py --version 2>/dev/null) | sed 's/^.*basic_read_statistics //; s/Using.*\$//' )
	END_VERSIONS
	"""


}