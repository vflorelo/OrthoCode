#!/usr/bin/python3
import pandas as pd
acc_go_df = pd.read_csv("All_species.go.tsv",sep="\t",dtype=str)
acc_og_df = pd.read_csv("acc_og_mappings.tsv",sep="\t",dtype=str)
acc_go_df = acc_go_df.join(acc_og_df.set_index("accession"),on="accession")
acc_go_df.to_csv("GO_master_file.tsv",sep="\t",index=False)