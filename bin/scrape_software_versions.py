#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'NGI-MethylSeq': ['v_ngi_methylseq.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'Bismark genomePrep': ['v_bismark_genome_preparation.txt', r"Bismark Genome Preparation Version: v(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Cutadapt': ['v_cutadapt.txt', r"(\S+)"],
    'Trim Galore!': ['v_trim_galore.txt', r"version (\S+)"],
    'Bismark': ['v_bismark.txt', r"Bismark Version: v(\S+)"],
    'Bismark Deduplication': ['v_deduplicate_bismark.txt', r"Deduplicator Version: v(\S+)"],
    'Bismark methXtract': ['v_bismark_methylation_extractor.txt', r"Bismark Extractor Version: v(\S+)"],
    'Bismark Report': ['v_bismark2report.txt', r"bismark2report version: v(\S+)"],
    'Bismark Summary': ['v_bismark2summary.txt', r"bismark2summary version: (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'Qualimap': ['v_qualimap.txt', r"QualiMap v.(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['NGI-MethylSeq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark genomePrep'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Cutadapt'] = '<span style="color:#999999;\">N/A</span>'
results['Trim Galore!'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark Deduplication'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark methXtract'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark Report'] = '<span style="color:#999999;\">N/A</span>'
results['Bismark Summary'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['Qualimap'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'NGI-MethylSeq Software Versions'
section_href: 'https://github.com/SciLifeLab/NGI-MethylSeq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
