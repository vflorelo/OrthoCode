#!/usr/bin/python3
import pandas as pd
count_tsv_file = "Orthogroups.GeneCount.tsv"
excl_tsv_file  = "Orthogroups_UnassignedGenes.tsv"
out_tsv_file   = "Orthogroups.FullGeneCount.tsv"
mat_tsv_file   = "Orthogroups.Presence.tsv"
count_df = pd.read_csv(count_tsv_file,sep="\t").fillna("")
excl_df  = pd.read_csv(excl_tsv_file, sep="\t",dtype=str).fillna("")
species_list = list(excl_df.columns[1:])
for species in species_list:
    excl_df[species] = excl_df[species].apply(lambda x: 1 if len(x) > 0 else 0)
join_df = pd.concat([count_df,excl_df],ignore_index=True)
join_df["Total"] = join_df[species_list].sum(axis=1)
join_df.to_csv(out_tsv_file,sep="\t",index=False)
for species in species_list:
    join_df[species] = join_df[species].apply(lambda x: 1 if x > 0 else 0)
join_df.drop(columns=["Total"]).to_csv(mat_tsv_file,sep="\t",index=False)