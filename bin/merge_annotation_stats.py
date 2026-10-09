#!/usr/bin/env python

import argparse
import rpy2
from rpy2.robjects.packages import importr
import rpy2.rinterface as rinterface
import pandas as pd
import numpy as np

as_char = rinterface.baseenv['as.character']
methylKit = importr('methylKit')
genomation = importr('genomation')



def main():
	parser = argparse.ArgumentParser(description='Create flat csv files from twist and repeat annotation RDS files')
	parser.add_argument('--methobj', required=True, help='Path to methobj file')
	parser.add_argument('--methobj_twist_annotation', required=True, help='Path to methobj  twist annotation file')
	parser.add_argument('--cpg_obj', required=True, help='Path to methobj file')
	parser.add_argument('--methobj_repeat_annotation', required=True, help='Path to output control stats CSV file')
	parser.add_argument('--twist_annotation_csv', required=True, help='Output CSV for twist annotation')
	parser.add_argument('--repeat_annotation_csv', required=True, help='Output CSV for repeat annotation')
	args = parser.parse_args()
	
	methobj = rpy2.robjects.r.readRDS(args.methobj)
	methobj_twist_annotation = rpy2.robjects.r.readRDS(args.methobj_twist_annotation)
	cpg_obj = rpy2.robjects.r.readRDS(args.cpg_obj)
	methobj_repeat_annotation = rpy2.robjects.r.readRDS(args.methobj_repeat_annotation)

	twist_stats_df = output_twist_annot_stats(methobj, methobj_twist_annotation,cpg_obj)
	twist_stats_df.to_csv(args.twist_annotation_csv)
	
	repeat_stats_df = output_repeat_annot_stats(methobj, methobj_repeat_annotation)
	repeat_stats_df.to_csv(args.repeat_annotation_csv)


# Twist Annotation Statistics ##########
def output_twist_annot_stats(methobj,methobj_twist_annotation,cpg_obj) -> None:
	rep_list = list(methylKit.getSampleID(methobj))
	cpg_obj_names = np.array(rpy2.robjects.r.names(cpg_obj))
	cpg_stats = [np.array(genomation.getFeatsWithTargetsStats(q, percentage=True))[0] for q in methobj_twist_annotation]
	cpg_stats_df = pd.DataFrame(np.array(np.array_split(cpg_stats,len(cpg_obj_names))), index=cpg_obj_names, columns=rep_list)
	return cpg_stats_df


# Repeat Annotation Statistics ##########
def output_repeat_annot_stats(methobj,methobj_repeat_annotation) -> None:
	annotation_description:str = 'UCSC RepeatMasker chm13v2.0_rmsk'
	rep_list = list(methylKit.getSampleID(methobj))
	cpg_stats = [np.array(genomation.getFeatsWithTargetsStats(q, percentage=True))[0] for q in methobj_repeat_annotation]
	cpg_stats_df = pd.DataFrame(cpg_stats, index=rep_list, columns=['Repeats from '+annotation_description])
	return cpg_stats_df
	

if __name__ == '__main__':
    main()