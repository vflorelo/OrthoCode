#!/usr/bin/python3
import pandas as pd
import sys
counts_file      = sys.argv[1]
unassigned_file  = sys.argv[2]
full_counts_file = "Orthogroups.FullGeneCount.tsv"
presence_file    = "Orthogroups.Presence.tsv"
counts_df        = pd.read_csv(counts_file,sep="\t").fillna("")
unassigned_df    = pd.read_csv(unassigned_file,sep="\t",dtype=str).fillna("")
species_list     = list(unassigned_df.columns[1:])
for species in species_list:
    unassigned_df[species] = unassigned_df[species].apply(lambda x: 1 if len(x) > 0 else 0)
join_df = pd.concat([counts_df,unassigned_df],ignore_index=True)
join_df["Total"] = join_df[species_list].sum(axis=1)
join_df.to_csv(full_counts_file,sep="\t",index=False)
for species in species_list:
    join_df[species] = join_df[species].apply(lambda x: 1 if x > 0 else 0)
join_df.drop(columns=["Total"]).to_csv(presence_file,sep="\t",index=False)