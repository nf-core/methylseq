#!/usr/bin/env python

"""
Merge tiled methylation counts (methylKit methylRawList from tile_methyl_counts.R) into one
table of beta values: one row per tile (chr_start_end), one column per sample.

Tiles outside chr1-22, X, Y, M are dropped. Tiles not covered in a sample are left empty.
"""

import argparse

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

KEEP_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

importr("methylKit")

TILES_TO_DATAFRAME = ro.r("""
function(obj) {
    ids <- getSampleID(obj)
    do.call(rbind, lapply(seq_along(obj), function(i) {
        d <- getData(obj[[i]])
        data.frame(
            sample = rep(ids[i], nrow(d)),
            chr    = as.character(d$chr),
            tile   = paste(d$chr, d$start, d$end, sep = "_"),
            numCs  = d$numCs,
            numTs  = d$numTs,
            stringsAsFactors = FALSE
        )
    }))
}
""")


def main():
    parser = argparse.ArgumentParser(description="Merge tiled beta values across samples")
    parser.add_argument("--tiled_methobj", required=True, help="Tiled methylRawList RDS file")
    parser.add_argument("--meta", required=False, help="Unused; sample ids are read from the RDS file")
    parser.add_argument("--output", required=True, help="Output CSV of tile beta values")
    args = parser.parse_args()

    tiled_methobj = ro.r["readRDS"](args.tiled_methobj)
    with localconverter(ro.default_converter + pandas2ri.converter):
        tiles = ro.conversion.rpy2py(TILES_TO_DATAFRAME(tiled_methobj))

    tile_beta_df = tile_beta_table(tiles)
    tile_beta_df.to_csv(args.output)


def tile_beta_table(tiles: pd.DataFrame) -> pd.DataFrame:
    tiles = tiles[tiles["chr"].isin(KEEP_CHROMOSOMES)].drop_duplicates(subset=["sample", "tile"])
    depth = tiles["numCs"] + tiles["numTs"]
    tiles = tiles.assign(beta=tiles["numCs"] / depth.where(depth > 0))
    sample_order = list(dict.fromkeys(tiles["sample"]))
    table = tiles.pivot(index="tile", columns="sample", values="beta")
    return table.reindex(columns=sample_order)


if __name__ == "__main__":
    main()
