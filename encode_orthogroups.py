#!/usr/bin/python3
import pandas as pd
import sys
tsv_in_file  = sys.argv[1] # presence matrix
tsv_out_file = sys.argv[2]
sp_list_file = sys.argv[3]
sp_list_fh   = open(sp_list_file, "r") 
sp_list_data = sp_list_fh.read() 
species_list = sp_list_data.split("\n") 
species_list = list(filter(None, species_list))
sp_list_fh.close() 
presence_df  = pd.read_csv(tsv_in_file,sep="\t").fillna("")
counter=0

for species in species_list:

    placement = 2**counter
    counter  += 1
    presence_df[species] = presence_df[species] * placement
presence_df["Total"] = presence_df[species_list].sum(axis=1)
presence_df.to_csv(tsv_out_file,sep="\t",index=False)