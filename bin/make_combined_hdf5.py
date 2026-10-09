#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import h5py
import re





def main():
	parser = argparse.ArgumentParser(description="Gather BedGraph files into single hdf5 file")
	parser.add_argument("--input_files", required=True, help="List of input bedGraph files (i.e. 'healthy_1.bedGraph healthy_2.bedGraph disease_1.bedGraph disease_2.bedGraph")
	parser.add_argument("--meta", required=True, help="Sample ID and group information metadata (i.e. [[id:healthy_1, group:healthy], [id:healthy_2, group:healthy], [id:disease_1, group:disease], [id:disease_2, group:disease]]")
	parser.add_argument("--id_file", required=True, help="Output CSV for id and file names")
	parser.add_argument("--hdf5", required=True, help="Output file for hdf5 output")
	args = parser.parse_args()

	input_files = args.input_files.split()

	rep_list = re.findall(r"id:([^,\]]+)", args.meta)

	file_id_df = pd.DataFrame.from_dict({
		"id": rep_list,
		"file": input_files
	})
	file_id_df.to_csv(args.id_file, index = None )

	with h5py.File(args.hdf5, "w") as tab_filtL_hdf5:
		dt = h5py.string_dtype(encoding="utf-8")
		
		for index, row in file_id_df.iterrows():
			insert_bedGraph_into_HDF5(
				row["file"], 
				row["id"],
				tab_filtL_hdf5,
				dt
			)



def insert_bedGraph_into_HDF5(bedGraph_file:str, rep_name:str, tab_filtL_hdf5, dt):
	all_chr_names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',  'chrY', 'chrM']
	
	tab = pd.read_table(bedGraph_file, header=None, skiprows=1)
	tab.columns = ['chr','coord0','coord1','coverage','meth_reads','unmeth_reads']
	read_num = tab.loc[:,'meth_reads'].add(tab.loc[:,'unmeth_reads']);
	tab_filt = tab.where(tab.loc[:,'chr'].isin(all_chr_names)).dropna();

	shape_param = (tab_filt.shape[0], tab_filt.shape[1]-1)
	dset = tab_filtL_hdf5.create_group(rep_name)
	dset = tab_filtL_hdf5.create_dataset(rep_name+'/'+tab_filt.columns[0], shape=(shape_param[0],), dtype=dt, compression='gzip')
	dset = tab_filtL_hdf5.create_dataset(rep_name+'/data', shape=shape_param, dtype='i', compression='gzip')

	tab_filtL_hdf5[rep_name+'/'+tab_filt.columns[0]][...]=tab_filt.iloc[:,0].values.astype(dt)
	tab_filtL_hdf5[rep_name+'/data'][...]=tab_filt.iloc[:,1:].values

if __name__ == "__main__":
    main()