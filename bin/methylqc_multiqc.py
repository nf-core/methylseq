#!/usr/bin/env python

"""
Convert EvolvDx methyl_qc outputs into MultiQC custom content (*_mqc.json).

Every input is optional; a section is written only for the inputs that are given.
All sections are grouped under one "EvolvDx methylation QC" heading in the report.
"""

import argparse
import json

import numpy as np
import pandas as pd

__version__ = "0.0.1"

PARENT = {
    "parent_id": "evolvdx_methylqc",
    "parent_name": "EvolvDx methylation QC",
    "parent_description": "Cohort-level methylation QC from the EvolvDx methyl_qc subworkflow (<code>--run_methylqc</code>).",
}

SPIKE_IN_NAMES = {"phage_lambda": "Lambda", "plasmid_puc19c": "pUC19"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control_stats", help="control_stats.csv from CONTROL_STATS")
    parser.add_argument("--read_metrics", help="metrics_table.csv from BASIC_READ_STATISTICS")
    parser.add_argument("--meth_distribution", help="tab_id_value_counts_df.csv from BASIC_READ_STATISTICS")
    parser.add_argument("--depth_distribution", help="rnf_counts_id_df.csv from BASIC_READ_STATISTICS")
    parser.add_argument("--mbias", help="merged_mbias.csv from BASIC_READ_STATISTICS")
    parser.add_argument("--chr_percent", help="chr_percent_df.csv from METH_HISTOGRAM_BY_CHROMOSOME")
    parser.add_argument("--twist_annotation", help="twist_annotation.csv from MERGE_ANNOTATION_STATS")
    parser.add_argument("--repeat_annotation", help="repeat_annotation.csv from MERGE_ANNOTATION_STATS")
    parser.add_argument("--tiled_beta", help="merge_tiled_stats.csv from MERGE_TILED_STATS")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    general_stats = {}
    general_headers = {}

    if args.control_stats:
        df = pd.read_csv(args.control_stats, index_col=0)
        for column in df.columns:
            label = SPIKE_IN_NAMES.get(column, column)
            key = f"{column}_meth"
            add_general_stat(general_stats, df[column], key)
            general_headers[key] = {
                "title": f"{label} meth",
                "description": f"Mean CpG methylation of the {label} spike-in control (%)",
                "min": 0, "max": 100, "suffix": "%", "scale": "OrRd", "format": "{:,.1f}",
            }

    if args.read_metrics:
        df = pd.read_csv(args.read_metrics, index_col=0)
        add_general_stat(general_stats, df["mean CpG methylation (%)"], "mean_cpg_meth")
        general_headers["mean_cpg_meth"] = {
            "title": "CpG meth",
            "description": "Mean methylation across covered CpGs (%)",
            "min": 0, "max": 100, "suffix": "%", "scale": "RdYlBu-rev", "format": "{:,.1f}",
        }
        add_general_stat(general_stats, df["median read #"], "median_cpg_depth")
        general_headers["median_cpg_depth"] = {
            "title": "CpG depth",
            "description": "Median read depth across covered CpGs",
            "min": 0, "suffix": "X", "scale": "BuPu", "format": "{:,.0f}",
        }
        write_section("methylqc_read_metrics", {
            "id": "methylqc_read_metrics",
            "section_name": "CpG read metrics",
            "description": "Per-sample CpG coverage and methylation summary from the MethylDackel bedGraphs (primary chromosomes only).",
            "plot_type": "table",
            "pconfig": {"id": "methylqc_read_metrics_table", "title": "methyl_qc: CpG read metrics", "namespace": "methyl_qc"},
            "headers": {
                "CpGs covered": {"format": "{:,.0f}", "scale": "Greens"},
                "median read #": {"title": "Median depth", "format": "{:,.1f}", "scale": "BuPu"},
                "mean read #": {"title": "Mean depth", "format": "{:,.1f}", "scale": "BuPu"},
                "% CpGs with < 10 reads": {"min": 0, "max": 100, "suffix": "%", "scale": "OrRd"},
                "mean CpG methylation (%)": {"min": 0, "max": 100, "suffix": "%", "scale": "RdYlBu-rev"},
                "% CpGs partially methylated": {"min": 0, "max": 100, "suffix": "%", "scale": "Purples"},
            },
            "data": table_data(df),
        })

    if args.repeat_annotation:
        df = pd.read_csv(args.repeat_annotation, index_col=0)
        add_general_stat(general_stats, df.iloc[:, 0], "repeat_cpg_pct")
        general_headers["repeat_cpg_pct"] = {
            "title": "% CpG in repeats",
            "description": f"Percent of covered CpGs overlapping {df.columns[0].replace('Repeats from ', '')}",
            "min": 0, "max": 100, "suffix": "%", "scale": "Blues", "format": "{:,.1f}",
        }

    if general_stats:
        write_section("methylqc_general_stats", {
            "id": "methylqc_general_stats",
            "plot_type": "generalstats",
            "pconfig": [{key: header} for key, header in general_headers.items()],
            "data": general_stats,
        })

    if args.meth_distribution:
        df = pd.read_csv(args.meth_distribution)
        data = {}
        for sample, d in df.groupby("replicate", sort=False):
            total = d["count"].sum()
            data[str(sample)] = [[int(x), round(100 * c / total, 4)] for x, c in zip(d["percent_methylation"], d["count"])]
        write_section("methylqc_meth_distribution", {
            "id": "methylqc_meth_distribution",
            "section_name": "CpG methylation distribution",
            "description": "Distribution of per-CpG methylation levels. Healthy cfDNA typically shows a large peak near 100% and a smaller one near 0%.",
            "plot_type": "linegraph",
            "pconfig": {
                "id": "methylqc_meth_distribution_plot",
                "title": "methyl_qc: CpG methylation distribution",
                "xlab": "CpG methylation (%)", "ylab": "% of covered CpGs",
                "xmin": 0, "xmax": 100, "ymin": 0,
            },
            "data": data,
        })

    if args.depth_distribution:
        df = pd.read_csv(args.depth_distribution)
        data = {}
        for sample, d in df.groupby("replicate", sort=False):
            data[str(sample)] = [[int(x), int(c)] for x, c in zip(d["read count"], d["number of read count occurrences"])]
        write_section("methylqc_depth_distribution", {
            "id": "methylqc_depth_distribution",
            "section_name": "CpG read depth distribution",
            "description": "Number of CpGs at each read depth (methylated + unmethylated reads).",
            "plot_type": "linegraph",
            "pconfig": {
                "id": "methylqc_depth_distribution_plot",
                "title": "methyl_qc: CpG read depth",
                "xlab": "Read depth", "ylab": "Number of CpGs",
                "xlog": True, "ymin": 0,
            },
            "data": data,
        })

    if args.mbias:
        df = pd.read_csv(args.mbias)
        df["pct"] = 100 * df["nMethylated"] / (df["nMethylated"] + df["nUnmethylated"]).replace(0, np.nan)
        datasets, labels = [], []
        for strand in ["OT", "OB"]:
            for read in sorted(df["Read"].unique()):
                d = df[(df["Strand"] == strand) & (df["Read"] == read)].sort_values("Position")
                if d.empty:
                    continue
                datasets.append({
                    str(sample): [[int(p), round(float(v), 3)] for p, v in zip(g["Position"], g["pct"]) if pd.notna(v)]
                    for sample, g in d.groupby("replicate", sort=False)
                })
                labels.append({"name": f"{strand} read {read}", "ylab": "% CpG methylated"})
        write_section("methylqc_mbias", {
            "id": "methylqc_mbias",
            "section_name": "M-bias",
            "description": "CpG methylation by position in the read (MethylDackel mbias). Lines should be flat; slopes at the read ends suggest trimming.",
            "plot_type": "linegraph",
            "pconfig": {
                "id": "methylqc_mbias_plot",
                "title": "methyl_qc: M-bias",
                "xlab": "Position in read (bp)", "ylab": "% CpG methylated",
                "ymin": 0, "ymax": 100,
                "data_labels": labels,
            },
            "data": datasets,
        })

    if args.chr_percent:
        df = pd.read_csv(args.chr_percent, index_col=0)
        pivot = df.pivot_table(index="Sample", columns="Chromosome", values="Percent CpGs Methylated", sort=False)
        cpgs = df.pivot_table(index="Sample", columns="Chromosome", values="CpG Count", aggfunc="sum", sort=False)
        # Drop chromosomes with no CpGs in any sample (e.g. absent spike-ins)
        pivot = pivot.loc[:, cpgs.sum(axis=0) > 0]
        write_section("methylqc_chr_methylation", {
            "id": "methylqc_chr_methylation",
            "section_name": "Methylation by chromosome",
            "description": "Percent of covered CpGs with any methylation, per chromosome and spike-in control.",
            "plot_type": "heatmap",
            "pconfig": {
                "id": "methylqc_chr_methylation_heatmap",
                "title": "methyl_qc: % CpGs methylated by chromosome",
                "xlab": "Chromosome", "ylab": "Sample",
                "min": 0, "max": 100, "square": False,
            },
            "xcats": [str(c) for c in pivot.columns],
            "ycats": [str(s) for s in pivot.index],
            "data": [[round(float(v), 2) if pd.notna(v) else None for v in row] for row in pivot.to_numpy()],
        })

    if args.twist_annotation:
        df = pd.read_csv(args.twist_annotation, index_col=0).T
        write_section("methylqc_twist_annotation", {
            "id": "methylqc_twist_annotation",
            "section_name": "Twist target annotation",
            "description": "Percent of covered CpGs overlapping each Twist methylome target annotation.",
            "plot_type": "table",
            "pconfig": {"id": "methylqc_twist_annotation_table", "title": "methyl_qc: Twist annotation", "namespace": "methyl_qc"},
            "headers": {str(c): {"min": 0, "max": 100, "suffix": "%", "scale": "Blues", "format": "{:,.1f}"} for c in df.columns},
            "data": table_data(df),
        })

    if args.tiled_beta:
        df = pd.read_csv(args.tiled_beta, index_col=0)
        bins = np.linspace(0, 1, 21)
        centres = (bins[:-1] + bins[1:]) / 2
        data = {}
        for sample in df.columns:
            values = df[sample].dropna().to_numpy()
            if values.size == 0:
                continue
            counts, _ = np.histogram(values, bins=bins)
            data[str(sample)] = [[round(float(x), 3), round(100 * c / values.size, 3)] for x, c in zip(centres, counts)]
        write_section("methylqc_tiled_beta", {
            "id": "methylqc_tiled_beta",
            "section_name": "Tiled beta values",
            "description": "Distribution of methylation beta values across 1 kb tiles.",
            "plot_type": "linegraph",
            "pconfig": {
                "id": "methylqc_tiled_beta_plot",
                "title": "methyl_qc: tile beta distribution",
                "xlab": "Beta value", "ylab": "% of tiles",
                "xmin": 0, "xmax": 1, "ymin": 0,
            },
            "data": data,
        })


def add_general_stat(general_stats: dict, series: pd.Series, key: str) -> None:
    for sample, value in series.items():
        if pd.notna(value):
            general_stats.setdefault(str(sample), {})[key] = float(value)


def table_data(df: pd.DataFrame) -> dict:
    return {
        str(sample): {str(k): (float(v) if pd.notna(v) else None) for k, v in row.items()}
        for sample, row in df.iterrows()
    }


# Line graph data is written as [x, y] pairs: MultiQC sorts dict keys as strings (1, 10, 100, 2, ...)
def write_section(name: str, content: dict) -> None:
    with open(f"{name}_mqc.json", "w") as fh:
        json.dump({**PARENT, **content}, fh, indent=2)


if __name__ == "__main__":
    main()
