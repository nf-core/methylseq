

process METH_HISTOGRAM_BY_CHROMOSOME {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/basic_read_statistics:0.0.1' :
        '062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/basic_read_statistics:0.0.1' }"
	
	input:
	path(bedGraph_files)
	val(bedGraph_meta)
	path(hg38_chrom_sizes)

	output:
	path "counts_df_long.csv", emit: counts_df_long
	path "chr_percent_df.csv", emit: chr_percent_df
	path "meth_histogram_by_chromosome_CpG_count.png"
	path "meth_histogram_by_chromosome_Percent_CpGs.png"
	path "versions.yml"      , emit: versions

	script:
	"""

	meth_histogram_by_chromosome.py \\
		--input_files '$bedGraph_files' \\
		--meta '$bedGraph_meta' \\
		--hg38_chrom_sizes $hg38_chrom_sizes \\
		--counts_df_long counts_df_long.csv \\
		--chr_percent_df chr_percent_df.csv
		
	

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    pandas: \$(python -c "import pandas; print(pandas.__version__)" 2>/dev/null)
	END_VERSIONS
	"""


}