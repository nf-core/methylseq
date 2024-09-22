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
        /bismark_alignment-(cnt|pct)\.(pdf|svg)/,
        /bismark_deduplication-(cnt|pct)\.(pdf|svg)/,
        /bismark-methylation-dp\.(pdf|svg|png)/,
        /bismark_mbias_CHG_R1\.(pdf|svg)/,
        /bismark_mbias_CHG_R2\.(pdf|svg)/,
        /bismark_mbias_CHH_R1\.(pdf|svg)/,
        /bismark_mbias_CHH_R2\.(pdf|svg)/,
        /bismark_mbias_CpG_R1\.(pdf|svg)/,
        /bismark_mbias_CpG_R2\.(pdf|svg)/,
        /bismark_strand_alignment-(cnt|pct)\.(pdf|svg)/,
        /cutadapt_filtered_reads_plot-(cnt|pct)\.(pdf|svg)/,
        /cutadapt_trimmed_sequences_plot_3_Counts\.(pdf|svg)/,
        /cutadapt_trimmed_sequences_plot_3_Obs_Exp\.(pdf|svg)/,
        /qualimap_coverage_histogram-(cnt|pct)\.(pdf|svg)/,
        /fastqc_adapter_content_plot\.(pdf|png|svg)/,
        /fastqc_overrepresented_sequences_plot(.*)?\.(pdf|svg)/,
        /fastqc_per_base_.*_plot(.*)?\.(pdf|png|svg)/,
        /fastqc_per_sequence_.*\.(pdf|svg)/,
        /fastqc_sequence_(counts|duplication_levels)_plot(.*)?\.(pdf|svg)/,
        /fastqc-status-check-.*\.(pdf|svg)/,
        /fastqc_top_overrepresented_sequences_table(.*)?\.(pdf|png|svg|txt)/,
        /fastqc_adapter_content_plot\.(pdf|png|svg)/,
        /fastqc_overrepresented_sequences_plot(.*)?\.(pdf|svg)/,
        /general_stats_table\.(pdf|png|svg)/,
        /qualimap_coverage_histogram\.(pdf|svg)/,
        /qualimap_gc_content\.(pdf|svg)/,
        /qualimap_genome_fraction\.(pdf|svg)/,
        /qualimap_insert_size\.(pdf|svg)/,
        /multiqc_data\.json/,
        /multiqc\.log/,
        /multiqc_sources\.txt/,
        /multiqc_fastqc\.txt/,
        /multiqc_fail_strand_check_table\.txt/,
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
