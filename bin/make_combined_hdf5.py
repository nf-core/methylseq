#!/usr/bin/env python

"""
Gather per-sample MethylDackel bedGraph files into a single HDF5 file.

Each sample is stored as a group named by its sample id with two datasets:
    <id>/chr   chromosome name of each CpG
    <id>/data  integer matrix: start, end, percent methylation, methylated reads, unmethylated reads
Only the primary human chromosomes (chr1-22, X, Y, M) are kept.
"""

import argparse
import re

import h5py
import pandas as pd

KEEP_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def main():
    parser = argparse.ArgumentParser(description="Gather bedGraph files into a single HDF5 file")
    parser.add_argument("--input_files", required=True, help="Space-separated bedGraph files, in the same order as --meta")
    parser.add_argument("--meta", required=True, help="Sample metadata, e.g. [[id:healthy_1, group:healthy], [id:disease_1, group:disease]]")
    parser.add_argument("--id_file", required=True, help="Output CSV mapping sample id to bedGraph file")
    parser.add_argument("--hdf5", required=True, help="Output HDF5 file")
    args = parser.parse_args()

    input_files = args.input_files.split()
    sample_ids = re.findall(r"id:([^,\]]+)", args.meta)
    if len(sample_ids) != len(input_files):
        raise SystemExit(f"Got {len(input_files)} bedGraph files but {len(sample_ids)} sample ids")

    pd.DataFrame({"id": sample_ids, "file": input_files}).to_csv(args.id_file, index=False)

    with h5py.File(args.hdf5, "w") as hdf5:
        for sample, bedgraph in zip(sample_ids, input_files):
            insert_bedgraph_into_hdf5(bedgraph, sample, hdf5)


def insert_bedgraph_into_hdf5(bedgraph_file: str, sample: str, hdf5: h5py.File) -> None:
    columns = ["chr", "start", "end", "percent_methylation", "meth_reads", "unmeth_reads"]
    try:
        tab = pd.read_table(bedgraph_file, header=None, skiprows=1, names=columns)
    except pd.errors.EmptyDataError:
        tab = pd.DataFrame(columns=columns)
    tab = tab[tab["chr"].isin(KEEP_CHROMOSOMES)]

    hdf5.create_group(sample)
    hdf5.create_dataset(
        f"{sample}/chr",
        data=tab["chr"].to_numpy().astype(object),
        dtype=h5py.string_dtype(encoding="utf-8"),
        compression="gzip",
    )
    hdf5.create_dataset(f"{sample}/data", data=tab.iloc[:, 1:].to_numpy(dtype="i"), compression="gzip")


if __name__ == "__main__":
    main()
