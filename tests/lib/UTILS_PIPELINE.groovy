class UTILS_PIPELINE {
    // These are the files to exclude when we want to snapshot
    static List<String> exclusionRegexesForUnstableFileContents = [

        // To exlude Bismark reports
        /.*_report\.txt/,
        /bismark_summary_report\.html/,
        /.*_report\.html/,

        // To exclude the pipeline_software_mqc_versions.yml file that contains the Nextflow version
        /nf_core_.*_software_mqc_versions\.yml/,

        // To exclude this folder that somehow is a file on stub tests
        /multiqc_plots/,

        // To exclude FASTQC reports
        /.*_raw\.(html|zip)/,
        /.*_fastqc\.(html|zip)/,

        // To exclude from the MultiQC reports
        /bismark_alignment-(cnt|pct)\.(pdf|png|svg)/,
        /bismark_deduplication-(cnt|pct)\.(pdf|png|svg)/,
        /bismark-methylation-dp\.(pdf|png|svg)/,
        /bismark_mbias_CHG_R1\.(pdf|png|svg)/,
        /bismark_mbias_CHG_R2\.(pdf|png|svg)/,
        /bismark_mbias_CHH_R1\.(pdf|png|svg)/,
        /bismark_mbias_CHH_R2\.(pdf|png|svg)/,
        /bismark_mbias_CpG_R1\.(pdf|png|svg)/,
        /bismark_mbias_CpG_R2\.(pdf|png|svg)/,
        /bismark_strand_alignment-(cnt|pct)\.(pdf|png|svg)/,
        /cutadapt_filtered_reads_plot-(cnt|pct)\.(pdf|png|svg)/,
        /cutadapt_trimmed_sequences_plot_3_Counts\.(pdf|png|svg)/,
        /cutadapt_trimmed_sequences_plot_3_Obs_Exp\.(pdf|png|svg)/,
        /fastqc_adapter_content_plot\.(pdf|png|svg)/,
        /fastqc_overrepresented_sequences_plot(.*)?\.(pdf|png|svg)/,
        /fastqc_per_base_.*_plot(.*)?\.(pdf|png|svg)/,
        /fastqc_per_sequence_.*\.(pdf|png|svg)/,
        /fastqc_sequence_(counts|duplication_levels)_plot(.*)?\.(pdf|png|svg)/,
        /fastqc-status-check-.*\.(pdf|png|svg)/,
        /fastqc_top_overrepresented_sequences_table(.*)?\.(pdf|png|svg|txt)/,
        /fastqc_overrepresented_sequences_plot(.*)?\.(pdf|png|svg)/,
        /general_stats_table\.(pdf|png|svg)/,
        /picard_deduplication-(cnt|pct)\.(pdf|png|svg)/,
        /qualimap_coverage_histogram\.(pdf|png|svg)/,
        /qualimap_gc_content\.(pdf|png|svg)/,
        /qualimap_genome_fraction\.(pdf|png|svg)/,
        /qualimap_insert_size\.(pdf|png|svg)/,
        /samtools_alignment_plot-(cnt|pct)\.(pdf|png|svg)/,
        /samtools-flagstat-dp_Percentage_of_total\.(pdf|png|svg)/,
        /samtools-flagstat-dp_Read_counts\.(pdf|png|svg)/,
        /samtools-stats-dp\.(pdf|png|svg)/,
        /multiqc_data\.json/,
        /multiqc\.log/,
        /multiqc_sources\.txt/,
        /multiqc_fastqc\.txt/,
        /multiqc_qualimap_bamqc_genome_results\.txt/,
        /multiqc_general_stats\.txt/,

        // To exclude Picard Markduplicates metrics
        /.*\.markdup\.sorted\.MarkDuplicates\.metrics\.txt/,

        // To exclude Qualimap files
        /.*\.(css|gif|js)/,
        /images_qualimapReport/,
        /raw_data_qualimapReport/,
        /qualimapReport\.html/,

        // To exclude from samtools stats
        /.*\.sorted\.bam\.(flagstat|idxstats|stats)/,

        // To exclude markdup
        /.*\.markdup\.sorted\.bam(\.bai)?/,

        // To exclude trimgalore
        /.*\.fastq\.gz_trimming_report\.txt/
    ]
}
