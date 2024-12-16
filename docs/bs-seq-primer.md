# Bisulfite Sequencing & Three-Base Aligners Primer

[Bisulfite sequencing](https://github.com/nf-core/methylseq) (**BS-seq**) is a widely-used technique to investigate **DNA methylation**, a crucial epigenetic modification that regulates gene expression without altering the DNA sequence.

## **Principle of Bisulfite Sequencing**

1. **Bisulfite Treatment**:

   - During bisulfite treatment, **non-methylated cytosines** in DNA are converted to **uracils**, while **methylated cytosines** remain unaffected
   - This chemical reaction forms the foundation for distinguishing methylated and unmethylated cytosines

2. **PCR Amplification**:

   - After bisulfite treatment, the DNA undergoes PCR amplification
   - During this step, uracils are converted into **thymines**
   - This conversion results in a reduced complexity of the DNA code, introducing computational challenges in downstream analysis

3. **DNA Strands in Sequencing**:
   - For a given genomic locus, bisulfite treatment and PCR amplification produce **four distinct DNA strands**:
     - Bisulfite-converted forward strand (OT)
     - Reverse complement of OT (CTOT)
     - Bisulfite-converted reverse strand (OB)
     - Reverse complement of OB (CTOB)
   - Depending on library type all four strands may end up in a sequencing library, adding to the complexity of analysis

## **Challenges in Mapping Bisulfite-Sequenced Reads**

Mapping bisulfite-treated sequences to a reference genome presents several computational challenges:

1. **Reduced DNA Complexity**:
   - Due to the conversion of cytosines to thymines, the DNA code becomes less diverse, increasing the likelihood of ambiguous alignments
2. **Multiple DNA Strands**:
   - The presence of four possible DNA strands (and their combinations) for each genomic locus increases the alignment search space
3. **Variable Methylation States**:
   - Each read can theoretically represent any possible methylation state for a locus, further complicating the alignment process

## **Applications of Bisulfite Sequencing**

Bisulfite sequencing is employed to study:

- **Genome-wide DNA methylation patterns**: Understanding epigenetic regulation of genes
- **Differential methylation**: Investigating differences between different samples (e.g. healthy versus diseased states)
- **Epigenetic inheritance**: Studying methylation changes across generations

---

## Three-Base Aligners: Bismark and BWA-Meth

Bisulfite sequencing converts unmethylated cytosines into uracils (later read as thymines during PCR), leading to a DNA code represented predominantly by three bases (A, G, and T).

Aligning these three-base reads against a standard reference genome is non-trivial, as the original strand identity and methylation state are obscured by the conversion.

To address these challenges, specialized “three-base aligners” have been developed to accurately map bisulfite-treated reads and infer their original strand context and methylation status. Here, we primarily summarise the aligners used in the nf-core/methylseq pipeline.

### Bismark ([docs](https://felixkrueger.github.io/Bismark/); [publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC3102221/))

- Bismark resolves strand ambiguity by performing up to four parallel alignments
- First, sequencing reads are converted _in silico_ to represent both forward and reverse strand conversions (C-to-T and G-to-A), mirroring fully bisulfite-converted versions of the reference genome
- Each read set is then aligned using Bowtie2 (alternatively HISAT2 or minimap2) against equally converted references; this enables support for indels, local alignments, and bisulfite converted RNA-seq-type or long reads (e.g. EM-seq using Nanopore or Pac Bio reads)
- By comparing these (up-to) four alignments, Bismark identifies each read’s correct strand origin

This approach allows Bismark to handle directional, PBAT, amplicon, and non-directional libraries robustly and to accurately align reads that represent partially methylated cytosines.

### BWA-Meth ([docs](https://github.com/brentp/bwa-meth); [publication](https://arxiv.org/abs/1401.1129))

- BWA-Meth adapts the BWA-MEM algorithm for bisulfite data, providing efficient, flexible alignments with support for indels, local alignments, and streaming-based workflows
- By converting reads _in silico_ on the fly, BWA-Meth eliminates the need for intermediate files, reducing temporary storage requirements and simplifying the overall process

The result is a fast, resource-efficient aligner that integrates smoothly with downstream analysis tools.

## At a glance

| Feature/Attribute                     | **Bismark**                                                                                                                                              | **BWA-Meth**                                                                                                                           |
| ------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| **Core algorithm & approach**         | Uses Bowtie2/HISAT2/minimap2 for 3-letter alignments                                                                                                     | Uses BWA-MEM for 3-letter alignments                                                                                                   |
| **Computational efficiency**          | More computationally intensive                                                                                                                           | Generally faster and more efficient due to single alignment step                                                                       |
| **Output format & methylation calls** | Produces SAM/BAM, deduplication, direct methylation calls (CpG, CHG, CHH), bedGraph/coverage and cytosine context files with integrated QC and reporting | Produces standard SAM/BAM; requires external tools (e.g., MethylDackel) for methylation calling and QC                                 |
| **Downstream analysis**               | Built-in methylation calling, filtering and (context-) reporting options                                                                                 | Requires external tools for methylation calling, offering more customization but added complexity                                      |
| **Recommended when**                  | Accuracy and detailed methylation context are top priority, and ample computational resources are available                                              | Speed (including GPU support), scalability, and modular workflows are essential, or when integrating with existing BWA-based pipelines |

**References**:

- **Bismark [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3102221/) and [official docs](https://felixkrueger.github.io/Bismark/)**
- **BWA-Meth [paper](https://arxiv.org/abs/1401.1129) and [GitHub](https://github.com/brentp/bwa-meth)**
- **[NVIDIA Parabricks `fq2bam_meth`](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_fq2bam_meth.html)**

## Additional Considerations for BS-seq:

#### Incomplete conversion and controls:

- **Incomplete bisulfite conversion**: Not all cytosines may be converted to uracil (and eventually to thymine) during the bisulfite treatment step. This partial conversion can lead to overestimation of methylation levels
- **Spike-In controls**: Often, fully unmethylated lambda DNA or other known-conversion controls are included to estimate and correct for the efficiency of the bisulfite reaction

#### Quality Control (QC) Metrics:

- **M-bias Plots**: Analysis tools can generate M-bias plots to visualize methylation bias across read positions. This helps identify technical artefacts, or non-uniform conversion at biased positiones at the ends of reads
- **Adapter- and quality-trimming**: Even more than standard sequencing data, bisulfite reads benefit from adapter removal and quality filtering. Tools like Trim Galore (often bundled with Bismark) ensure higher-quality alignments and more accurate methylation calls

#### Post-alignment deduplication and bias correction:

- **PCR duplication removal**: Bisulfite sequencing libraries often contain PCR duplicates. Removing these duplicates is crucial to avoid artificially inflated methylation signals (relevant for statistics in downstream analysis)

#### Context-specific methylation:

- While the focus often lies on CpG methylation (the most studied context in mammals), accurate mapping and downstream tools can also call methylation in CHG and CHH contexts. This is particularly relevant for plant genomes where non-CpG methylation is biologically significant

#### Large-scale analysis and cloud computing:

- As whole-genome bisulfite sequencing (WGBS) datasets grow larger, computational efficiency, scalability, and memory usage become critical factors
- **Cloud-based computing**: Tools that integrate into cloud-based workflows (e.g., AWS, GCP) and employ parallelization or GPU-acceleration (e.g., NVIDIA Parabricks for `fq2bam_meth`) can significantly reduce runtime for large projects

#### Long-read bisulfite sequencing:

- Although short-read Illumina sequencing predominates, emerging methods such as enzymatically converted methylation sequencing (EM-seq) are possible using long-read platforms (e.g., PacBio or Oxford Nanopore)
- Aligning and calling methylation from longer reads introduces different challenges and opportunities, including better resolution of repetitive regions and phasing of haplotypes

#### Experimental Design and Biological Replicates:

- Careful experimental design, including replicates and controls, is essential. BS-seq can be influenced by technical variability, and replication ensures statistical robustness in methylation calling and differential methylation analyses
