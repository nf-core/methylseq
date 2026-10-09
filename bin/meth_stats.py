#!/usr/bin/env python
import argparse
from collections import defaultdict as dd
import pandas as pd
import pathlib

__version__ = '0.0.1'

def get_controls(file: pd.DataFrame):
    df = pd.read_csv(file, sep='\t', skiprows=1, names=['contig', 'start', 'stop', 'percent','C', 'T'])
    df = df[(df['contig'] == 'phage_lambda') | (df['contig'] == 'plasmid_puc19c')]
    df['cov'] = df['C'] + df['T']
    return df

def sn_sv_metrics(frac_c, con_type, thresh=0.99):
    if con_type == 'plasmid_puc19c':
        if frac_c >= thresh:
            return 'TP'
        return 'FN'
    if con_type == 'phage_lambda':
        if frac_c <= thresh:
            return 'TN'
        return 'FP'
    
def process_controls(methyl_frame, min_cov=10, thresh=0.99):
    this_frame = methyl_frame[(methyl_frame['contig'] == 'phage_lambda') | (methyl_frame['contig'] == 'plasmid_puc19c')]
    this_frame['cov'] = this_frame['C'] + this_frame['T']
    this_frame = this_frame[this_frame['cov'] >= min_cov]
    this_frame['frac_c'] = this_frame['C']/this_frame['cov']
    this_frame['frac_t'] = this_frame['T']/this_frame['cov']
    if this_frame.empty:
        this_frame['metrics'] = None
    else:
        this_frame['metrics'] = this_frame.apply(lambda x: sn_sv_metrics(x.frac_c, x.contig, thresh), axis=1)
    return this_frame


class ReportSN_SV:
    def __init__(self, methylframe):
        self.methyl_frame = methylframe
        
        self.mdict = dd(int)
        
        self.tp, self.tn, self.fn, self.fn = 0, 0, 0, 0

        self.sn = 0
        self.sv = 0
        self.tpn = 0
        
        self.__processed = False
        
    def __repr__(self):
        if self.__processed:
            return (f"Sensitivity: {self.sn}\n"
                    f"Specificity: {self.sv}\n"
                    f"Negative Predictive Value: {self.tpn}\n"
                    f"TP: {self.tp}\n"
                    f"FP: {self.fp}\n"
                    f"TN: {self.tn}\n"
                    f"FN: {self.fn}")
        return "Not Processed Yet"
    
    def get_sn_sv(self):
        mdict = dd(int)
        recs = self.methyl_frame.groupby('metrics')['cov'].count().to_dict()
        for k, v in recs.items():
            mdict[k] = v

        self.tp, self.fp, self.fn, self.tn = mdict['TP'], mdict['FP'], mdict['FN'], mdict['TN']

        try:
            self.sn  = self.tp / (self.tp + self.fn)
        except ZeroDivisionError:
            self.sn = 'div0'
        try:
            self.sv  = self.tp / (self.tp + self.fp)
        except ZeroDivisionError:
            self.sv = 'div0'

        try:
            self.tpn = self.tn / (self.tn + self.fn)
        except ZeroDivisionError:
            self.tpn = 'div0'
        self.__processed = True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Method for finding sensitivity and specificity of contols for any given sample")
    parser.add_argument('bedgraph', type=pathlib.Path)
    parser.add_argument('out', type=str)
    parser.add_argument('--threshold', type=float, default=0.99)
    parser.add_argument('--mincov', type=int, default=10)
    
    args = parser.parse_args()
    
    df = get_controls(args.bedgraph)
    pc_df = process_controls(df, min_cov=args.mincov, thresh=args.threshold)
    snsv = ReportSN_SV(pc_df)
    snsv.get_sn_sv()
    with open(args.out, 'w') as fo:
        fo.write(snsv.__repr__())
