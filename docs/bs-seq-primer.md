# Bisulfite Sequencing & Three-Base Aligners Primer

Bisulfite sequencing (**BS-Seq**) is a widely-used technique to investigate **DNA methylation**, a crucial epigenetic modification that regulates gene expression without altering the DNA sequence.

## **Principle of Bisulfite Sequencing**

1. **Bisulfite Treatment**:

   - During bisulfite treatment, **non-methylated cytosines** in DNA are converted to **uracils**, while **methylated cytosines** remain unaffected.
   - This chemical reaction forms the foundation for distinguishing methylated and unmethylated cytosines.

2. **PCR Amplification**:

   - After bisulfite treatment, the DNA undergoes PCR amplification.
   - During this step, uracils are converted into **thymines**.
   - This conversion results in a reduced complexity of the DNA code, introducing computational challenges in downstream analysis.

3. **DNA Strands in Sequencing**:
   - For a given genomic locus, bisulfite treatment and PCR amplification produce **four distinct DNA strands**:
     - Forward strand (original sequence)
     - Reverse complement strand
     - Bisulfite-converted forward strand
     - Bisulfite-converted reverse complement strand
   - All four strands may end up in a sequencing library, adding to the complexity of analysis.

## **Challenges in Mapping Bisulfite-Sequenced Reads**

Mapping bisulfite-treated sequences to a reference genome presents several computational challenges:

1. **Reduced DNA Complexity**:

   - Due to the conversion of cytosines to thymines, the DNA code becomes less diverse, increasing the likelihood of ambiguous alignments.

2. **Multiple DNA Strands**:

   - The presence of four possible DNA strands (and their combinations) for each genomic locus increases the alignment search space.

3. **Variable Methylation States**:
   - Each read can theoretically represent any possible methylation state for a locus, further complicating the alignment process.

## **Applications of Bisulfite Sequencing**

Bisulfite sequencing is employed to study:

- **Genome-wide DNA methylation patterns**: Understanding epigenetic regulation of genes.
- **Differential methylation**: Investigating differences between healthy and diseased states (e.g., cancer).
- **Epigenetic inheritance**: Studying methylation changes across generations.

---

## Three-Base Aligners: Bismark and BWA-Meth

Bisulfite sequencing (BS-Seq) converts unmethylated cytosines into uracils (later read as thymines during PCR), leading to a DNA code represented predominantly by three bases (A, G, and T).

Aligning these three-base reads against a standard reference genome is non-trivial, as the original strand identity and methylation state are obscured by the conversion.

To address these challenges, specialized “three-base aligners” have been developed to accurately map bisulfite-treated reads and infer their original strand context and methylation status. Here, we primarily highlight the aligner used in the nf-core/methylseq pipeline

### [Bismark](https://pmc.ncbi.nlm.nih.gov/articles/PMC3102221/)

- Bismark resolves strand ambiguity by performing four parallel alignments.

- First, sequencing reads are converted in silico to represent both forward and reverse strand conversions (C-to-T and G-to-A), mirroring fully bisulfite-converted versions of the reference genome.

- Each read set is then aligned using Bowtie (or Bowtie2/HISAT2) against equally converted references.

- By comparing these four alignments, Bismark identifies each read’s correct strand origin

This approach allows Bismark to handle both directional and non-directional libraries robustly and to accurately align reads that represent partially methylated cytosines.

### [BWA-Meth](https://arxiv.org/abs/1401.1129)

- BWA-Meth adapts the BWA-MEM algorithm for bisulfite data, providing efficient, flexible alignments with support for indels, local alignments, and streaming-based workflows.

- By converting reads _in silico_ on the fly, BWA-Meth eliminates the need for intermediate files, reducing storage requirements and simplifying the overall process.

Its robustness minimizes the need for aggressive quality trimming, maintaining data integrity and streamlining the alignment pipeline. The result is a fast, resource-efficient aligner that integrates smoothly with downstream analysis tools.

| Feature/Attribute                     | **Bismark**                                                                                                 | **BWA-Meth**                                                                                                                           |
| ------------------------------------- | ----------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| **Core Algorithm & Approach**         | Uses Bowtie2/HISAT2 for a 3-way alignment (C-to-T, G-to-A, and original) tailored for bisulfite data        | Adapts BWA-MEM for bisulfite alignment, aligning reads directly to a pre-converted reference genome                                    |
| **Computational Efficiency**          | More computationally intensive, multiple alignments per read, higher memory/time cost                       | Generally faster, more efficient, and lower memory usage due to single-step alignment to the converted genome                          |
| **Output Format & Methylation Calls** | Produces SAM/BAM and direct methylation calls (CpG, CHG, CHH) with integrated QC and reporting              | Produces standard SAM/BAM; requires external tools (e.g., MethylDackel) for methylation calling and QC                                 |
| **Complexity Handling**               | Handles up to four strands (original and converted), providing robust methylation context coverage          | Primarily relies on the converted reference, slightly less complexity in handling multiple strands                                     |
| **Use Cases & Scale**                 | Ideal for small/medium-scale projects or when detailed methylation contexts (CpG/CHG/CHH) matter            | Suitable for large-scale/complex projects where computational efficiency and scalability are key                                       |
| **Error Tolerance & Sensitivity**     | High sensitivity to bisulfite-induced mismatches, potentially more accurate for diverse states              | Good sensitivity, though slightly less tolerant to certain bisulfite mismatches compared to Bismark                                    |
| **Downstream Analysis**               | Built-in methylation calling and reporting simplifies pipeline development                                  | Requires external tools for methylation calling, offering more customization but added complexity                                      |
| **Recommended When**                  | Accuracy and detailed methylation context are top priority, and ample computational resources are available | Speed (including GPU support), scalability, and modular workflows are essential, or when integrating with existing BWA-based pipelines |

**References**:

- **Bismark [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3102221/) and [official docs](https://felixkrueger.github.io/Bismark/)**
- **BWA-Meth [paper](https://arxiv.org/abs/1401.1129) and [github](https://github.com/brentp/bwa-meth)**
- **[NVIDIA Parabricks fq2bam_meth](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_fq2bam_meth.html)**

## Additional Considerations in Bisulfite Sequencing:

1. Incomplete Conversion and Controls:

- Incomplete Bisulfite Conversion: Not all cytosines may be converted to uracil (and eventually to thymine) during the bisulfite treatment step. This partial conversion can lead to overestimation of methylation levels.
- Spike-In Controls: Often, fully unmethylated lambda DNA or other known-conversion controls are included to estimate and correct for the efficiency of the bisulfite reaction.

2. Quality Control (QC) Metrics:

- M-bias Plots: Analysis tools can generate M-bias plots to visualize methylation bias across read positions. This helps identify technical artifacts such as non-uniform conversion or biased coverage at the ends of reads.
- Adapter and Quality Trimming: Like standard sequencing data, bisulfite reads require adapter removal and quality filtering. Tools like Trim Galore! (often bundled with Bismark) ensure higher-quality alignments and more accurate methylation calls.

3. Post-Alignment Deduplication and Bias Correction:

- PCR Duplication Removal: Bisulfite sequencing libraries often contain PCR duplicates. Removing these duplicates is crucial to avoid artificially inflated methylation signals.
- Strand-Specific Bias: Some methylation calling methods consider strand-specific bias. Ensuring that duplicates are removed and that data are balanced across strands is important for accurate methylation quantification.

4. Context-Specific Methylation:

- While the focus often lies on CpG methylation (the most studied context in animals), accurate mapping and downstream tools can also call methylation in CHG and CHH contexts. This is particularly relevant in plant genomes and certain vertebrate tissues where non-CpG methylation is biologically significant.

5. Large-Scale Analysis and Cloud Computing:

- As whole-genome bisulfite sequencing (WGBS) datasets grow larger, computational efficiency, scalability, and memory usage become critical factors.
- Cloud-Based Computing: Tools that integrate into cloud-based workflows (e.g., AWS, GCP) and employ parallelization or GPU-acceleration (e.g., NVIDIA Parabricks for fq2bam_meth) can significantly reduce runtime for large projects.

6. Long-Read Bisulfite Sequencing:

- Although short-read Illumina sequencing predominates, emerging methods using long-read platforms (e.g., PacBio or Oxford Nanopore) also support bisulfite-treated samples.
- Aligning and calling methylation from longer reads introduces different challenges and opportunities, including better resolution of repetitive regions and phasing of haplotypes.

7. Integration with Methylation-Specific Analysis Tools:

- Post-alignment, specialized tools like Bismark, MethylDackel, methyKit can refine methylation calls, produce reports, and integrate with downstream differential methylation analysis software.

8. Experimental Design and Biological Replicates:

- Careful experimental design, including replicates and controls, is essential. Bisulfite sequencing can be influenced by technical variability, and replication ensures statistical robustness in methylation calling and differential methylation analyses.
