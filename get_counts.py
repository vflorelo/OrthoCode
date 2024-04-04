#!/usr/bin/python3
import pandas as pd
import sys
tsv_file = sys.argv[1]
out_file = sys.argv[2]
tsv_df   = pd.read_csv(tsv_file,sep="\t",dtype=str).fillna("")
sp_list  = list(tsv_df.columns[1:])
for species in sp_list:
    tsv_df[species] = tsv_df[species].apply(lambda x: len(list(filter(None,x.split(",")))))
tsv_df.to_csv(out_file,sep="\t",index=False)