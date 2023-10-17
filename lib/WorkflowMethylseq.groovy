//
// This file holds several functions specific to the workflow/methylseq.nf in the nf-core/methylseq pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowMethylseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        genomeExistsError(params, log)


        if (!params.fasta) {
            if (!params.bismark_index || params.aligner != 'bismark') {
                Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            }
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        // FIXME Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
        def citation_text = [
                "Tools used in the workflow included:",
                "FastQC (Andrews 2010),",
                "Trim Galore! (Krueger)",
                "Bismark (Krueger 2011)",
                "bwa-meth (Pedersen 2014)",
                "Picard (Broad Institute 2019)",
                "Qualimap (Okonechnikov 2015)",
                "Preseq (Daley 2013)",
                "MultiQC (Ewels et al. 2016)",
                "."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        // FIXME Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
        def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).",
            "Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354",
            "https://www.bioinformatics.babraham.ac.uk/projects/trim_galore",
            "Felix Krueger, Simon R. Andrews, Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications, Bioinformatics, Volume 27, Issue 11, 1 June 2011, Pages 1571–1572, doi: /10.1093/bioinformatics/btr167",
            "Pedersen, Brent S. and Eyring, Kenneth and De, Subhajyoti and Yang, Ivana V. and Schwartz, David A. Fast and accurate alignment of long bisulfite-seq reads, arXiv:1401.1129, doi: 10.48550/arXiv.1401.1129",
            "Picard Tools, Broad Institute. <http://broadinstitute.github.io/picard/>",
            "Konstantin Okonechnikov, Ana Conesa, Fernando García-Alcalde, Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data, Bioinformatics, Volume 32, Issue 2, 15 January 2016, Pages 292–294, doi: 10.1093/bioinformatics/btv566",
            "Daley, T., Smith, A. Predicting the molecular complexity of sequencing libraries. Nat Methods 10, 325–327 (2013). doi: 10.1038/nmeth.2375</li>",
        ].join('</li> <li>').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }
}
