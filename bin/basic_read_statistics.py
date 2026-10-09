#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import h5py
import re
import matplotlib.pyplot as plt
import seaborn as sns


__version__ = '0.0.1'


def main():
	parser = argparse.ArgumentParser(description="Gather BedGraph files into single hdf5 file")
	parser.add_argument("--bedGraph_hdf5", required=True, help="Single hdf5 file with bedGraph data from all samples")
	parser.add_argument("--mbias_files", required=True, help="List of input bedGraph files (i.e. 'healthy_1.mbias.txt healthy_2.mbias.txt disease_1.mbias.txt disease_2.mbias.txt")
	parser.add_argument("--meta", required=True, help="Sample ID and group information metadata (i.e. [[id:healthy_1, group:healthy], [id:healthy_2, group:healthy], [id:disease_1, group:disease], [id:disease_2, group:disease]]")
	parser.add_argument("--id_file", required=True, help="Output CSV for id and file names")
	parser.add_argument("--merged_mbias", required=True, help="Output merged mbias CSV file")
	parser.add_argument('--version', action='version', version=__version__)
	
	args = parser.parse_args()

	mbias_files = args.mbias_files.split()

	rep_list = re.findall(r"id:([^,\]]+)", args.meta)

	file_id_df = pd.DataFrame.from_dict({
		"id": rep_list,
		"file": mbias_files
	})
	file_id_df.to_csv(args.id_file, index = None )

	
	mbias_dataL = [pd.read_csv(mbias_file, sep = '\t') for mbias_file in mbias_files]
	rep_series_mbiasL = map_rep_id_series(mbias_dataL, file_id_df["id"])
	mbias_data_idL = concat_ids_onto_df_list(mbias_dataL, rep_series_mbiasL)
	mbias_data_id_df = pd.concat(mbias_data_idL, axis=0)
	mbias_data_id_df.to_csv(args.merged_mbias, index=None)

	tab_id_value_countsL = []

	rnf_countsL = []
	rnf_medianL = []
	rnf_meanL = []


	with h5py.File(args.bedGraph_hdf5, "r") as tab_filtL_hdf5:
		for index, row in file_id_df.iterrows():
			tab_id = tab_filtL_hdf5[row["id"] +'/data'][:,2]
			tab_id_value_counts = np.unique(tab_id, return_counts=True)
			replicate_name_array = [row["id"] for _ in tab_id_value_counts[0]]
			tab_id_value_counts_df = pd.DataFrame(data={"coverage":tab_id_value_counts[0], "count":tab_id_value_counts[1], "replicate":replicate_name_array})
			tab_id_value_countsL.append(tab_id_value_counts_df)

			read_num_filt = tab_filtL_hdf5[row["id"]+'/data'][:,3]+tab_filtL_hdf5[row["id"]+'/data'][:,4]
			rnf_value_counts = np.unique(read_num_filt, return_counts=True)
			rnf_value_counts_df = pd.DataFrame(rnf_value_counts, index=['read count','number of read count occurrences']).T
			rnf_countsL.append(rnf_value_counts_df)

			rnf_medianL.append(np.median(read_num_filt))
			rnf_meanL.append(read_num_filt.mean())
	
	tab_id_value_counts_df = pd.concat(tab_id_value_countsL).reset_index()
	tab_id_value_counts_df.to_csv("tab_id_value_counts_df.csv", index = None)

	rep_series_rnfL = map_rep_id_series(rnf_countsL, file_id_df["id"])
	rnf_counts_idL = concat_ids_onto_df_list(rnf_countsL, rep_series_rnfL)
	rnf_counts_id_df = pd.concat(rnf_counts_idL, axis=0)
	rnf_counts_id_df.to_csv("rnf_counts_id_df.csv", index = None)

	report_table_precision = 3
	tab_count_by_rep = tab_id_value_counts_df.loc[:,["coverage", "replicate"]].groupby("replicate").count()
	tab_count_by_rep.to_csv("tab_count_by_rep.csv")
	
	lt100_tab_count_by_rep = tab_id_value_counts_df.where(tab_id_value_counts_df.loc[:,'coverage'].lt(100)).dropna().groupby("replicate").count()
	lt100_tab_count_by_rep.to_csv("lt100_tab_count_by_rep.csv")
	
	percent_partial_read_coverage = lt100_tab_count_by_rep.div(tab_count_by_rep).mul(100).loc[:,'coverage']
	percent_partial_read_coverage = percent_partial_read_coverage.map(lambda x: np.round(x, report_table_precision))
	percent_partial_read_coverage.name = '% of total reads that are partial'
	metrics_table = pd.concat([percent_partial_read_coverage, pd.DataFrame({'median read #':rnf_medianL, 'mean read #':rnf_meanL}, index=file_id_df["id"])], axis=1)
	metrics_table.to_csv("metrics_table.csv")

	plot_basic_read_statistics(rnf_counts_id_df,mbias_data_id_df,tab_id_value_counts_df,metrics_table)




def rep_id_series(tab: pd.DataFrame, rep_name: str, series_name: str = 'replicate') -> pd.Series:
	''' tab: pd.DataFrame \n rep_name: str \n Make a Pandas Series of the replicate name rep_name for a replicate with DataFrame tab that shares its index. '''
	return pd.Series([rep_name for _ in tab.index], index = tab.index, name = series_name)
def map_rep_id_series(tab_filtL: list, rep_list: list) -> list:
    ''' tab_filtL: A list of pd.DataFrame that is the chromosome-filtered BedGraph data \n rep_list: List of all rep names in pd.Series format. \n Create a list containing the replicate ID to all of the pd.DataFrame in tab_filtL. Assumes that the list of pd.DataFrame and the list of replicates are the same length and order.'''
    return list(map(lambda tab, rep_name: rep_id_series(tab, rep_name), tab_filtL, rep_list))
def concat_ids_onto_df_list(tabL: list, rep_seriesL: list) -> list:
    ''' tab_filtL: A list of pd.DataFrame that is the chromosome-filtered BedGraph data \n rep_list: List of all rep names in pd.Series format. \n
        Append a column containing the replicate ID to all of the pd.DataFrame in tab_filtL.
        Assumes that the list of pd.DataFrame and the list of replicate pd.Series are the same
        length and same order.'''
    return list(map(lambda tab, rep_series: pd.concat([tab, rep_series], axis=1), tabL, rep_seriesL))

def plot_basic_read_statistics(rnf_counts_id_df:pd.DataFrame,
                                mbias_data_id_df:pd.DataFrame,
                                tab_id_value_counts_df:pd.DataFrame,
                                metrics_table:pd.DataFrame,
                                rnf_apart:bool = True,
                                rep_designator:str = 'replicate',
                                statistic:str = 'coverage') -> None:

	# Combined Coverage Distribution
	# print_header('Combined Coverage Distributions')
	f,ax = plt.subplots()
	sns.lineplot(data=tab_id_value_counts_df, x=statistic, y='count', hue=rep_designator, palette='Set2');
	plt.yscale('log')
	sns.despine(left=True, bottom=True);
	f.savefig("combined_coverage_distributions.png")


    # Coverage distribution
    # print_header('Individual Coverage Distributions')
    # g = sns.FacetGrid(tab_id_value_counts_df, col=rep_designator, height=2, aspect=2, col_wrap=2)
    # g.map(sns.lineplot, statistic, 'count')
    # plt.yscale('log')

    # Read count occurrence histograms
	# print_header('Read Count Occurrence Frequencies')
	sns.lineplot(data=rnf_counts_id_df, x='read count', y='number of read count occurrences', hue='replicate', palette='Set2')
	sns.scatterplot(data=rnf_counts_id_df, x='read count', y='number of read count occurrences', s=20, color='gray')
	plt.xscale('log')
	if rnf_apart:
		g = sns.FacetGrid(rnf_counts_id_df.reset_index(), col=rep_designator, height=4, aspect=2, col_wrap=2)
		g.map(sns.lineplot, 'read count', 'number of read count occurrences')
		g.map(sns.scatterplot, 'read count', 'number of read count occurrences', s=20, color='gray')
		plt.xscale('log')
	plt.savefig("read_count_occurence_frequencies.png")

    # Read bias plots: Methylated OT Strand
	# print_header('Read Bias in Methylated OT Strand')
	g = sns.FacetGrid(mbias_data_id_df.where(mbias_data_id_df.loc[:,'Strand'].eq('OT')), col=rep_designator, height=4, aspect=2, col_wrap=2, hue='Read', palette='Set1', sharey=False)
	g.map(sns.lineplot, 'Position', 'nMethylated')
	plt.ylim(bottom=0)
	g.add_legend()
	plt.savefig("read_bias_in_methylated_OT_strand.png")

    # Read bias plots: Methylated OB Strand
	# print_header('Read Bias in Methylated OB Strand')
	g = sns.FacetGrid(mbias_data_id_df.where(mbias_data_id_df.loc[:,'Strand'].eq('OB')), col=rep_designator, height=4, aspect=2, col_wrap=2, hue='Read', palette='Set1', sharey=False)
	g.map(sns.lineplot, 'Position', 'nMethylated')
	plt.ylim(bottom=0)
	g.add_legend()
	plt.savefig("read_bias_in_methylated_OB_strand.png")

    # Read bias plots: Unmethylated OT Strand
	# print_header('Read Bias in Unmethylated OT Strand')
	g = sns.FacetGrid(mbias_data_id_df.where(mbias_data_id_df.loc[:,'Strand'].eq('OT')), col=rep_designator, height=4, aspect=2, col_wrap=2, hue='Read', palette='Set1', sharey=False)
	g.map(sns.lineplot, 'Position', 'nUnmethylated')
	plt.ylim(bottom=0)
	g.add_legend()
	plt.savefig("read_bias_in_unmethylated_OT_strand.png")

    # Read bias plots: Unmethylated OB Strand
	# print_header('Read Bias in Unmethylated OB Strand')
	g = sns.FacetGrid(mbias_data_id_df.where(mbias_data_id_df.loc[:,'Strand'].eq('OB')), col=rep_designator, height=4, aspect=2, col_wrap=2, hue='Read', palette='Set1', sharey=False)
	g.map(sns.lineplot, 'Position', 'nUnmethylated')
	plt.ylim(bottom=0)
	g.add_legend()
	plt.savefig("read_bias_in_unmethylated_OB_strand.png")
# Output catch-all read statistics ##########

if __name__ == "__main__":
    main()