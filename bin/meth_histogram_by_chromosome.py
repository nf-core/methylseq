#!/usr/bin/env python

"""
Per-chromosome CpG counts and methylation from MethylDackel bedGraph files.

A CpG counts as methylated when its methylation level is above 0 (partially or fully methylated).

Outputs:
    --counts_df_long  long table of methylated CpG counts: index (chromosome), replicate, count
    --chr_percent_df  long table: Sample, Chromosome, Methylated CpG Count, CpG Count, Percent CpGs Methylated
    meth_histogram_by_chromosome_CpG_count.png, meth_histogram_by_chromosome_Percent_CpGs.png
"""

import argparse
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

KEEP_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
CONTROL_CHROMOSOMES = ["phage_lambda", "plasmid_puc19c"]


def main():
    parser = argparse.ArgumentParser(description="Per-chromosome CpG methylation from bedGraph files")
    parser.add_argument("--input_files", required=True, help="Space-separated bedGraph files, in the same order as --meta")
    parser.add_argument("--meta", required=True, help="Sample metadata, e.g. [[id:healthy_1, group:healthy], [id:disease_1, group:disease]]")
    parser.add_argument("--hg38_chrom_sizes", required=True, help="Chromosome sizes file (two columns: name, length)")
    parser.add_argument("--counts_df_long", required=True, help="Output CSV of methylated CpG counts per chromosome")
    parser.add_argument("--chr_percent_df", required=True, help="Output CSV of CpG counts and percent methylated per chromosome")
    args = parser.parse_args()

    input_files = args.input_files.split()
    sample_ids = re.findall(r"id:([^,\]]+)", args.meta)
    if len(sample_ids) != len(input_files):
        raise SystemExit(f"Got {len(input_files)} bedGraph files but {len(sample_ids)} sample ids")

    # Only chromosomes present in the sizes file are reported (plus the spike-in controls)
    chrom_sizes = pd.read_csv(args.hg38_chrom_sizes, sep="\t", header=None, usecols=[0], names=["chromosome"])
    chromosomes = [c for c in KEEP_CHROMOSOMES if c in set(chrom_sizes["chromosome"])] + CONTROL_CHROMOSOMES

    chr_percent_df = pd.concat(
        [chromosome_stats(f, sample, chromosomes) for sample, f in zip(sample_ids, input_files)],
        ignore_index=True,
    )
    chr_percent_df.to_csv(args.chr_percent_df)

    counts_df_long = chr_percent_df.rename(
        columns={"Chromosome": "index", "Sample": "replicate", "Methylated CpG Count": "count"}
    )[["index", "replicate", "count"]]
    counts_df_long.to_csv(args.counts_df_long)

    plot_chromosome_histograms(chr_percent_df)


def chromosome_stats(bedgraph_file: str, sample: str, chromosomes: list) -> pd.DataFrame:
    try:
        tab = pd.read_table(bedgraph_file, header=None, skiprows=1, usecols=[0, 3], names=["chr", "percent_methylation"])
    except pd.errors.EmptyDataError:
        tab = pd.DataFrame({"chr": pd.Series(dtype=str), "percent_methylation": pd.Series(dtype=float)})
    cpg_count = tab.groupby("chr").size().reindex(chromosomes, fill_value=0)
    meth_count = (tab["percent_methylation"] > 0).groupby(tab["chr"]).sum().reindex(chromosomes, fill_value=0).astype(int)
    percent = (100 * meth_count / cpg_count.where(cpg_count > 0)).fillna(0)
    return pd.DataFrame({
        "Sample": sample,
        "Chromosome": chromosomes,
        "Methylated CpG Count": meth_count.to_numpy(),
        "CpG Count": cpg_count.to_numpy(),
        "Percent CpGs Methylated": percent.to_numpy(),
    })


def plot_chromosome_histograms(chr_percent_df: pd.DataFrame, dims=(8, 4), rotation=45) -> None:
    for column, filename in [
        ("CpG Count", "meth_histogram_by_chromosome_CpG_count.png"),
        ("Percent CpGs Methylated", "meth_histogram_by_chromosome_Percent_CpGs.png"),
    ]:
        f, ax = plt.subplots(figsize=dims)
        sns.barplot(data=chr_percent_df, y=column, x="Chromosome", hue="Sample", palette="Set2", ax=ax)
        sns.despine()
        plt.xticks(rotation=rotation, ha="right")
        ax.legend(bbox_to_anchor=(1, 1))
        f.savefig(filename, bbox_inches="tight")
        plt.close(f)


if __name__ == "__main__":
    main()
