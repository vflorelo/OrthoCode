#!/usr/bin/python3
import pandas as pd
#import sys
import re
#enc_og_tsv    = sys.argv[1]
#og_codes      = sys.argv[2]
#sp_list_file  = sys.argv[3]
#lineage_file  = sys.argv[4]
enc_og_tsv    = "Orthogroups.Encoded.tsv"
og_codes      = "og_codes.txt"
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
lineage_list  = lineages_df["lineage_name"].values.flatten().tolist()
def get_matrix_code(matrix_str,matrix_type):
    matrix_list = list(matrix_str)
    dot_count  = 0
    one_count  = 0
    zero_count = 0
    for matrix_element in range(len(matrix_list)):
        if (matrix_list[matrix_element] == "."):
            dot_count += 1
        if (matrix_list[matrix_element] == "0"):
            zero_count += 1
        if (matrix_list[matrix_element] == "1"):
            one_count += 1
    if (matrix_type == "core_gain"):    # dot tolerance is len(species_list) - 1
        low_limit  = dot_count - 1
        high_limit = dot_count + 1
    elif(matrix_type == "half_gain"):   # dot tolerance is len(species_list) / 2
        low_limit  = dot_count / 2
        high_limit = dot_count - 1
    elif(matrix_type == "core_loss"):   # one tolerance is len(species_list) - 1
        low_limit  = one_count - 1
        high_limit = one_count + 1
    elif(matrix_type == "half_loss"):   # dot tolerance is len(species_list) / 2
        low_limit  = dot_count / 2
        high_limit = dot_count + 1
    mat_df              = og_codes_df.copy()
    mat_df["match"]     = mat_df["matrix"].apply(lambda x: re.findall(matrix_str,x)[0] if re.findall(matrix_str,x) else "")
    index_list          = mat_df["match"].apply(lambda x: True if x else False)
    mat_df              = mat_df[index_list]
    mat_df["match_sum"] = mat_df["match"].apply(lambda x: sum(list(map(int,list(x)))))
    mat_df              = mat_df[[mat_df["match_sum"]>=low_limit and mat_df["match_sum"]<high_limit]]
    code_list           = mat_df["code"].values.flatten().tolist()
    return code_list
for lineage in lineage_list:
    cur_species_list    = lineages_df[lineages_df["lineage_name"]==lineage]["species_list"].to_list()[0].split(",")
    cur_index_list      = []
    excl_gain_base_list = ["0"] * species_count
    half_gain_base_list = ["0"] * species_count
    core_gain_base_list = ["0"] * species_count
    excl_loss_base_list = ["1"] * species_count
    half_loss_base_list = ["1"] * species_count
    core_loss_base_list = ["."] * species_count
    for cur_species in cur_species_list:
        cur_index = species_list.index(cur_species)
        excl_gain_base_list[cur_index] = "1"
        half_gain_base_list[cur_index] = "."
        core_gain_base_list[cur_index] = "."
        excl_loss_base_list[cur_index] = "0"
        half_loss_base_list[cur_index] = "."
        core_loss_base_list[cur_index] = "0"
    core_gain_matrix = ''.join(list_item for list_item in core_gain_base_list)
    half_gain_matrix = ''.join(list_item for list_item in half_gain_base_list)
    excl_gain_matrix = ''.join(list_item for list_item in excl_gain_base_list)
    excl_gain_matrix = excl_gain_matrix[::-1]
    core_loss_matrix = ''.join(list_item for list_item in core_loss_base_list)
    half_loss_matrix = ''.join(list_item for list_item in half_loss_base_list)
    excl_loss_matrix = ''.join(list_item for list_item in excl_loss_base_list)
    excl_loss_matrix = excl_loss_matrix[::-1]
    core_gain_code   = get_matrix_code(core_gain_matrix,"core_gain")
    half_gain_code   = get_matrix_code(half_gain_matrix,"half_gain")
    excl_gain_code   = int(excl_gain_matrix,base=2)
    core_loss_code   = get_matrix_code(core_loss_matrix,"core_loss")
    half_loss_code   = get_matrix_code(half_loss_matrix,"half_loss")
    excl_loss_code   = int(excl_loss_matrix,base=2)
    core_gain_index  = enc_og_df["Total"].isin(core_gain_code)
    half_gain_index  = enc_og_df["Total"].isin(half_gain_code)
    excl_gain_index  = enc_og_df["Total"] == str(excl_gain_code)
    core_loss_index  = enc_og_df["Total"].isin(core_loss_code)
    half_loss_index  = enc_og_df["Total"].isin(half_loss_code)
    excl_loss_index  = enc_og_df["Total"] == str(excl_loss_code)
    core_gain_list   = enc_og_df[core_gain_index]["Orthogroup"].values.flatten().tolist()
    half_gain_list   = enc_og_df[half_gain_index]["Orthogroup"].values.flatten().tolist()
    excl_gain_list   = enc_og_df[excl_gain_index]["Orthogroup"].values.flatten().tolist()
    core_loss_list   = enc_og_df[core_loss_index]["Orthogroup"].values.flatten().tolist()
    half_loss_list   = enc_og_df[half_loss_index]["Orthogroup"].values.flatten().tolist()
    excl_loss_list   = enc_og_df[excl_loss_index]["Orthogroup"].values.flatten().tolist()