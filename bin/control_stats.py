#!/usr/bin/env python

import argparse
import rpy2
from rpy2.robjects.packages import importr
import pandas as pd
import numpy as np

methylKit = importr('methylKit')



def main():
	parser = argparse.ArgumentParser(description='Create stats for methylation controls')
	parser.add_argument('--methobj', required=True, help='Path to methobj file')
	parser.add_argument('--stats_out', required=True, help='Path to output control stats CSV file')
	args = parser.parse_args()
	
	methobj = rpy2.robjects.r.readRDS(args.methobj)
	control_stats_df = create_control_stats_table(methobj)
	control_stats_df.to_csv(args.stats_out)



def create_control_stats_table(methobj, spike_in_controls:list = ['phage_lambda', 'plasmid_puc19c']) -> pd.DataFrame:
	chr_index = 0
	coverage_index = 4
	numCs_index = 5
	pct_methylatedL = []
	rep_list = list(methylKit.getSampleID(methobj))
	for i,this_rep in enumerate(rep_list):
		pct_methylatedL_aux = []
		for this_control in spike_in_controls:
			array_selection = np.array(methobj[i][chr_index]) == this_control
			if np.any(array_selection):
				pct_methylated = 100*np.array(methobj[i][numCs_index])[array_selection] / np.array(methobj[i][coverage_index])[array_selection]
				pct_methylated = np.mean(pct_methylated)
				pct_methylatedL_aux.append(pct_methylated)
			else:
				pct_methylatedL_aux.append(0)
		pct_methylatedL.append(pct_methylatedL_aux)
	return pd.DataFrame(pct_methylatedL, index=rep_list, columns=spike_in_controls)

if __name__ == '__main__':
    main()