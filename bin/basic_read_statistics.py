#!/usr/bin/env python

"""
Per-sample read statistics from the merged MethylDackel bedGraph HDF5 and M-bias files.

Outputs (current directory):
    merged_mbias.csv            MethylDackel M-bias tables for all samples, with a 'replicate' column
    tab_id_value_counts_df.csv  Distribution of per-CpG methylation (%) per sample
    rnf_counts_id_df.csv        Distribution of per-CpG read depth per sample
    metrics_table.csv           One row of summary metrics per sample
    *.png                       Static plots of the above
"""

import argparse
import re

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

__version__ = "0.0.2"

# CpGs with fewer reads than this are counted as low depth in metrics_table.csv
LOW_DEPTH_THRESHOLD = 10


def main():
    parser = argparse.ArgumentParser(description="Per-sample read statistics from bedGraph HDF5 and M-bias files")
    parser.add_argument("--bedGraph_hdf5", required=True, help="Single HDF5 file with bedGraph data from all samples")
    parser.add_argument("--mbias_files", required=True, help="Space-separated MethylDackel M-bias files, in the same order as --meta")
    parser.add_argument("--meta", required=True, help="Sample metadata, e.g. [[id:healthy_1, group:healthy], [id:disease_1, group:disease]]")
    parser.add_argument("--id_file", required=True, help="Output CSV mapping sample id to M-bias file")
    parser.add_argument("--merged_mbias", required=True, help="Output merged M-bias CSV file")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    mbias_files = args.mbias_files.split()
    sample_ids = re.findall(r"id:([^,\]]+)", args.meta)
    if len(sample_ids) != len(mbias_files):
        raise SystemExit(f"Got {len(mbias_files)} M-bias files but {len(sample_ids)} sample ids")

    file_id_df = pd.DataFrame({"id": sample_ids, "file": mbias_files})
    file_id_df.to_csv(args.id_file, index=False)

    mbias_df = pd.concat(
        [pd.read_csv(f, sep="\t").assign(replicate=sample) for sample, f in zip(sample_ids, mbias_files)],
        ignore_index=True,
    )
    mbias_df.to_csv(args.merged_mbias, index=False)

    meth_distL = []
    depth_distL = []
    metricsL = []
    with h5py.File(args.bedGraph_hdf5, "r") as hdf5:
        for sample in sample_ids:
            # data columns: start, end, percent methylation, methylated reads, unmethylated reads
            data = hdf5[sample + "/data"][:]
            meth_reads = data[:, 3].astype(float)
            depth = meth_reads + data[:, 4]
            covered = depth > 0
            pct_meth = np.divide(100 * meth_reads, depth, out=np.zeros_like(depth), where=covered)

            values, counts = np.unique(np.rint(pct_meth[covered]).astype(int), return_counts=True)
            meth_distL.append(pd.DataFrame({"percent_methylation": values, "count": counts, "replicate": sample}))

            values, counts = np.unique(depth[covered].astype(int), return_counts=True)
            depth_distL.append(pd.DataFrame({"read count": values, "number of read count occurrences": counts, "replicate": sample}))

            n_cpg = int(covered.sum())
            metricsL.append({
                "sample": sample,
                "CpGs covered": n_cpg,
                "median read #": float(np.median(depth[covered])) if n_cpg else 0.0,
                "mean read #": float(depth[covered].mean()) if n_cpg else 0.0,
                f"% CpGs with < {LOW_DEPTH_THRESHOLD} reads": pct(np.sum(depth[covered] < LOW_DEPTH_THRESHOLD), n_cpg),
                "mean CpG methylation (%)": float(pct_meth[covered].mean()) if n_cpg else 0.0,
                "% CpGs partially methylated": pct(np.sum((pct_meth > 0) & (pct_meth < 100) & covered), n_cpg),
            })

    meth_dist_df = pd.concat(meth_distL, ignore_index=True)
    meth_dist_df.to_csv("tab_id_value_counts_df.csv", index=False)

    depth_dist_df = pd.concat(depth_distL, ignore_index=True)
    depth_dist_df.to_csv("rnf_counts_id_df.csv", index=False)

    metrics_table = pd.DataFrame(metricsL).set_index("sample").round(3)
    metrics_table.to_csv("metrics_table.csv")

    plot_basic_read_statistics(depth_dist_df, mbias_df, meth_dist_df)


def pct(numerator, denominator) -> float:
    return float(100 * numerator / denominator) if denominator else 0.0


def plot_basic_read_statistics(depth_dist_df: pd.DataFrame, mbias_df: pd.DataFrame, meth_dist_df: pd.DataFrame) -> None:
    f, ax = plt.subplots()
    sns.lineplot(data=meth_dist_df, x="percent_methylation", y="count", hue="replicate", palette="Set2", ax=ax)
    ax.set_yscale("log")
    ax.set_xlabel("CpG methylation (%)")
    ax.set_ylabel("Number of CpGs")
    sns.despine(left=True, bottom=True)
    f.savefig("combined_coverage_distributions.png", bbox_inches="tight")
    plt.close(f)

    g = sns.FacetGrid(depth_dist_df, col="replicate", height=4, aspect=2, col_wrap=2)
    g.map(sns.lineplot, "read count", "number of read count occurrences")
    g.map(sns.scatterplot, "read count", "number of read count occurrences", s=20, color="gray")
    g.set(xscale="log")
    g.savefig("read_count_occurence_frequencies.png")
    plt.close(g.figure)

    for strand in ["OT", "OB"]:
        for column, label in [("nMethylated", "methylated"), ("nUnmethylated", "unmethylated")]:
            g = sns.FacetGrid(
                mbias_df[mbias_df["Strand"] == strand],
                col="replicate", height=4, aspect=2, col_wrap=2, hue="Read", palette="Set1", sharey=False,
            )
            g.map(sns.lineplot, "Position", column)
            g.set(ylim=(0, None))
            g.add_legend()
            g.savefig(f"read_bias_in_{label}_{strand}_strand.png")
            plt.close(g.figure)


if __name__ == "__main__":
    main()
