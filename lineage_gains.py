#!/usr/bin/python3
import pandas as pd
import re
enc_og_tsv    = "Orthogroups.Encoded.tsv"
og_codes      = "og_codes.tsv"
sp_list_file  = "species_list"
lineages_file = "lineages.tsv"
sp_list_fh    = open(sp_list_file, "r") 
sp_list_data  = sp_list_fh.read() 
species_list  = sp_list_data.split("\n") 
species_list  = list(filter(None, species_list))
species_count = len(species_list)
sp_list_fh.close()
enc_og_df     = pd.read_csv(enc_og_tsv,   sep="\t",dtype=str)
og_codes_df   = pd.read_csv(og_codes,     sep="\t",dtype=str)
lineages_df   = pd.read_csv(lineages_file,sep="\t",dtype=str)
lineage_list  = lineages_df["lineage"].values.flatten().tolist()
def get_code_list(matrix_str):
    matrix_list = list(matrix_str)
    dot_count  = 0
    one_count  = 0
    zero_count = 0
    for matrix_element in range(len(matrix_list)):
        if (matrix_list[matrix_element] == "."):
            dot_count += 1
        if (matrix_list[matrix_element] == "1"):
            one_count += 1
    low_limit  = dot_count - tolerance
    high_limit = one_count + dot_count
    mat_df              = og_codes_df.copy()
    mat_df["match"]     = mat_df["matrix"].apply(lambda x: re.findall(matrix_str,x)[0] if re.findall(matrix_str,x) else "")
    index_list          = mat_df["match"].apply(lambda x: True if x else False)
    mat_df              = mat_df[index_list]
    mat_df["match_sum"] = mat_df["match"].apply(lambda x: sum(list(map(int,list(x)))))
    mat_df              = mat_df[(mat_df["match_sum"]>=low_limit) & (mat_df["match_sum"]<=high_limit)]
    code_list           = mat_df["code"].values.flatten().tolist()
    return code_list
for lineage in lineage_list:
    fixed_species_list = lineages_df[lineages_df["lineage"]==lineage]["fixed"].to_list()[0].split(",")
    cur_species_list   = lineages_df[lineages_df["lineage"]==lineage]["species"].to_list()[0].split(",")
    tolerance          = lineages_df[lineages_df["lineage"]==lineage]["tolerance"]
    cur_index_list     = []
    gain_base_list     = ["0"] * species_count
    for fixed_species in fixed_species_list:
        cur_index = species_list.index(fixed_species)
        gain_base_list[cur_index] = "1"
    for cur_species in cur_species_list:
        cur_index = species_list.index(cur_species)
        gain_base_list[cur_index] = "."
    gain_matrix = ''.join(list_item for list_item in gain_base_list)
    gain_code_list   = get_code_list(gain_matrix)
    core_gain_index  = enc_og_df["Total"].isin(gain_code_list)