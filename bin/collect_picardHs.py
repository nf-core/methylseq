#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import pathlib
import seaborn as sns
import os
import yaml

__version__ = '0.0.1'

# Picard on off target barplot
def generate_on_off_target_tsv(picard_path, sample_name, header):
    
    # Format picard result data
    try:
        picard_sample = pd.read_csv(f'{picard_path}/{sample_name}.CollectHsMetrics.coverage_metrics', sep='\t',
                    skiprows=6, nrows=1)
    except:
        return
    barplot_data = picard_sample[['ON_BAIT_BASES', 'NEAR_BAIT_BASES', 'OFF_BAIT_BASES']]
    barplot_data = barplot_data.melt(var_name='type', value_name='bases').set_index('type')
    barplot_data.index.name = None
    
    
    # Add header
    with open(header, "r") as input:
        with open(f"{sample_name}_on_off_target_mqc.tsv", "w") as output:
            for line in input:
                output.write(line)
            for i, row in barplot_data.iterrows():
                output.write('\n')
                output.write(f'{row.name}\t{row[0]}')

def wrapper_tsv(picard_path, sample_sheet, header):
    samplesheet = pd.read_csv(sample_sheet, sep=',')
    for sample in samplesheet['sample']:
        generate_on_off_target_tsv(picard_path, sample, header)


# Picard uniformity lineplot
def picard_uniformity(picard_path):
    column_header = None
    for file in os.listdir(picard_path):
        if file.endswith(".coverage_metrics"):
            sample_name = file.split('.')[0]
            file_path = os.path.join(picard_path, file)
            df = pd.read_csv(file_path, sep='\t',skiprows=6, nrows=1)

            # Get data for uniformity line plot
            lineplot_col = list(df.columns[df.columns.to_series().str.contains('PCT_TARGET_BASES')])
            lineplot_col_filtered = []
            for col in lineplot_col:
                if int(col.replace('PCT_TARGET_BASES_', '').replace('X', '')) < 1000:
                    lineplot_col_filtered.append(col)
            lineplot_data = df[lineplot_col_filtered]

            # Format data for lineplot
            if column_header is None:
                column_header = [float(x.replace('PCT_TARGET_BASES_', '').replace('X', '')) for x in lineplot_col_filtered]
                lineplot_df = pd.DataFrame(columns = column_header)
            lineplot_data.columns = column_header
            lineplot_data.index = [sample_name]
            lineplot_df = pd.concat([lineplot_df, lineplot_data])

    # Build dictionary for lineplot data:
    lineplot_df = lineplot_df*100
    lineplot_dict = dict()
    for i, row in lineplot_df.iterrows():
        row_dict = dict()
        for indices in row.index:
            row_dict[indices] = row[indices]
        lineplot_dict[row.name] = row_dict
    lineplot_data_dict = {'data': lineplot_dict}


    with open('uniformity.yml', 'w') as outfile:
        yaml.dump(lineplot_data_dict, outfile, default_flow_style=False)

# Picard other metrics in table
def picard_other_metrics(picard_path, table_config):
    table_df = pd.DataFrame()
    probe_set = None
    for file in os.listdir(picard_path):
        if file.endswith(".coverage_metrics"):
            sample_name = file.split('.')[0]
            file_path = os.path.join(picard_path, file)
            df = pd.read_csv(file_path, sep='\t',skiprows=6, nrows=1)

            # Get data for table output
            table_data = df[['FOLD_ENRICHMENT','FOLD_80_BASE_PENALTY','MEAN_TARGET_COVERAGE']]
            table_data.index = [sample_name]
            try:
                table_data['FOLD_80_BASE_PENALTY'] = round(table_data['FOLD_80_BASE_PENALTY'],2)
            except:
                pass
            table_df = pd.concat([table_df, table_data])
            
            # Get probe set name
            if probe_set is None:
                probe_set = str(df['BAIT_SET'].unique()[0])

    # Build dictionary for table data:
    table_dict = dict()
    for i, row in table_df.iterrows():
        row_dict = dict()
        for indices in row.index:
            try:
                row_dict[indices] = row[indices].item()
            except:
                row_dict[indices] = row[indices]
        table_dict[row.name] = row_dict
    
    # Convert table config to a dictionary to be edited:
    with open(table_config, 'r') as stream:
        table = yaml.safe_load(stream)
    
    # Add probe set name:
    table['description'] = f'Bait set used to generate all picard result is: {probe_set}' 
    
    # Add data:
    table['data'] = table_dict
    
    with open('picard_other_metrics_mqc.yml', 'w') as outfile:
        yaml.dump(table, outfile, default_flow_style=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Picard Hs visuals"""
    )
    parser.add_argument("picard_path", type=pathlib.Path)
    parser.add_argument("sample_sheet", type=pathlib.Path)
    parser.add_argument("bargraph_header", type=pathlib.Path)
    parser.add_argument("table_config", type=pathlib.Path)
    parser.add_argument('--version', action='version', version=__version__)

    args = parser.parse_args()
    wrapper_tsv(args.picard_path, args.sample_sheet, args.bargraph_header)
    picard_uniformity(args.picard_path)
    picard_other_metrics(args.picard_path, args.table_config)

