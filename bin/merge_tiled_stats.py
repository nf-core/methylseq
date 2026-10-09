#!/usr/bin/env python

import argparse
import rpy2
from rpy2.robjects.packages import importr
import pandas as pd
import numpy as np
import numpy as np
import re
import rpy2.rinterface as rinterface

as_char = rinterface.baseenv['as.character'] 	
rutils = importr('utils')



def main():
	parser = argparse.ArgumentParser(description='Calculate and merge tiled beta methylation scores across samples')
	parser.add_argument('--tiled_methobj', required=True, help='Path to tiled methobj RDS file')
	parser.add_argument("--meta", required=True, help="Sample ID and group information metadata (i.e. [[id:healthy_1, group:healthy], [id:healthy_2, group:healthy], [id:disease_1, group:disease], [id:disease_2, group:disease]]")
	parser.add_argument('--output', required=True, help='Path to output merged tiled beta values CSV file')
	args = parser.parse_args()
	
	tiled_methobj = rpy2.robjects.r.readRDS(args.tiled_methobj)
	rep_list = re.findall(r"id:([^,\]]+)", args.meta)

	tile_beta_df = tile_beta_to_csv(rep_list,tiled_methobj)
	tile_beta_df.to_csv(args.output)

def beta(meth_reads:list, unmeth_reads:list, a:float) -> np.array:
    ''' Calculate beta values for DNA methylation. 
        Beta is essentially the fraction methylated in a sample at a specific base.

        Parameters:
        meth_reads: List of methylated read counts.
        unmeth_reads: List of unmethylated read counts.
        a: Stabilization parameter, often set to either 0 or 100 in the literature.

        The variables meth_reads and unmeth_reads must be the same length. Each element of the list corresponds to a coordinate.
    '''

    numerator = np.array(meth_reads)
    denominator = np.array(meth_reads) + np.array(unmeth_reads) + a
 
    try:
        return np.array(numerator/denominator)
    except:
        print('Incorrect format for the lists or the a parameter.')

def tile_beta_to_csv(rep_list:list,
                     methobj_tiles #:rpy2.robjects.vectors.ListVector,
                     ) -> pd.DataFrame:
	'''
	Function taking rpy2 object methobj_tiles and saving the tiles beta values to a csv named tile_beta_df_fn.
	'''
	keep_chr_names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',  'chrY', 'chrM']
    

    # Indices for (tiled) methylation R object
	chromosome_column_index = 0
	start_column_index = 1
	end_column_index = 2
	strand_column_index = 3
	coverage_column_index = 4
	meth_column_index = 5
	unmeth_column_index = 6


	tile_betaL = []
	for i,tiled_counts_obj in enumerate(methobj_tiles):
		# if verbose:
		# 	print('Tiling beta values for '+rep_list[i])
		# We must be careful to index the R object with 1-based indexing, so the range of indices is shifted by 1.
		chr_list = [as_char(tiled_counts_obj[chromosome_column_index].rx[int(i)])[0] for i in np.arange(1, len(tiled_counts_obj[chromosome_column_index])+1)]
		coord_matrix = pd.DataFrame({
 	 	 	'chr': chr_list,
 	 	 	'start': list(map(str,np.array(tiled_counts_obj[start_column_index]))),
 	 	 	'end': list(map(str,np.array(tiled_counts_obj[end_column_index]))),
 	 	 	'strand': list(map(str,np.array(tiled_counts_obj[strand_column_index])))})

		is_in_keep_chr = coord_matrix.loc[:,'chr'].isin(keep_chr_names)
		meth_column = pd.Series(tiled_counts_obj[meth_column_index]).where(is_in_keep_chr).dropna().values
		unmeth_column = pd.Series(tiled_counts_obj[unmeth_column_index]).where(is_in_keep_chr).dropna().values
		coord_matrix = coord_matrix.where(is_in_keep_chr).dropna().values
		tile_index =  list(map(lambda x: '_'.join(x), coord_matrix))

		tile_beta_aux = pd.Series(beta(meth_column, unmeth_column, 0), index=tile_index, name = rep_list[i])

		tile_betaL.append(tile_beta_aux)

	dups = [q.index.where(q.index.duplicated()).dropna().values for q in tile_betaL]
	tile_betaL_filt = [tile_betaL[i].mask(tile_betaL[i].index.isin(q)).dropna() for i,q in enumerate(dups)]
	tile_beta_df = pd.concat(tile_betaL_filt, axis=1).fillna(0)
	
	return tile_beta_df

if __name__ == '__main__':
    main()