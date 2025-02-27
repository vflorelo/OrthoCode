#!/usr/bin/python3
import itertools
import numpy   as np
import pandas  as pd
import re
import seaborn as sns
import sys
from   tqdm    import tqdm
from   pyDOE2  import ff2n
lineages_file     = sys.argv[1]
species_list_file = sys.argv[2]
og_tsv_file       = sys.argv[3]
og_enc_tsv_file   = sys.argv[4]
og_codes_file     = sys.argv[5]
species_list_fh   = open(species_list_file, "r") 
species_list_data = species_list_fh.read() 
species_list      = species_list_data.split("\n") 
species_list      = list(filter(None, species_list))
species_list_fh.close()
species_count     = len(species_list)
encoded_og_df     = pd.read_csv(og_enc_tsv_file,sep="\t",dtype=str).fillna(0)
full_og_df        = pd.read_csv(og_tsv_file,    sep="\t",dtype=str).fillna('')
og_codes_df       = pd.read_csv(og_codes_file,  sep="\t",dtype=str).fillna('')
og_matrix_list    = og_codes_df["matrix"].values.flatten().tolist()
lineages_df       = pd.read_csv(lineages_file,  sep="\t",dtype=str).fillna('')
lineage_list      = lineages_df["lineage"].values.flatten().tolist()
def filter_matrix_list(cur_og_matrix_list,sp_list,list_type,tolerance):
    dot_base_str = ["."] * species_count
    target_values   = []
    target_str_list = []
    for sp_index,sp_name in enumerate(species_list):
        if sp_name in sp_list:
            target_values.append(sp_index)
    target_len = len(sp_list)
    sub_mats = [list(i) for i in itertools.product([0, 1], repeat=target_len)]
    mats = []
    for sub_mat in sub_mats:
        if(sum(sub_mat) >= target_len - tolerance):
            mats.append(sub_mat)
    del(sub_mats)
    re_list = []
    for mat in mats:
        dot_str = dot_base_str.copy()
        sp_mat = list(set(map(lambda x, y: x * y, mat, target_values)))
        if (0 in sp_mat):
            sp_mat.remove(0)
        dot_arr = np.array(dot_str)
        if (list_type == "fixed" or list_type == "include"):
            dot_arr[sp_mat] = 1
        elif (list_type == "exclude"):
            dot_arr[sp_mat] = 0
        dot_list = dot_arr.tolist()
        re_str = "".join(dot_list)
        re_list.append(re_str)
        del(dot_str)
        del(dot_arr)
        del(re_str)
    del(mats)
    re_count  = len(re_list)
    if(re_count >= 100):
        re_mod    = re_count % 100
        bin_count = int(((re_count - re_mod) / 100) + 1)
    else:
        bin_count = 1
    for bin_num in tqdm(range(bin_count)):
        start_pos   = bin_num * 100
        re_sub_list = re_list[start_pos:][0:100]
        re_obj_str  = "|".join(list(re_sub_list))
        re_obj_str  = "\\b("+re_obj_str+")\\b"
        re_obj      = re.compile(re_obj_str)
        strings     = list(filter(lambda x: re_obj.match(x), cur_og_matrix_list))
        target_str_list += strings
        target_str_list = list(set(target_str_list))
        del(re_obj)
        del(strings)
    target_str_list = list(set(target_str_list))
    return target_str_list
for lineage in lineage_list:
    print("Processing "+lineage+" lineage\n")
    lin_df         = lineages_df.copy()
    lin_df         = lin_df[lin_df["lineage"]==lineage]
    fixed_list     = lin_df["fixed"].values.flatten().tolist()[0].split(",")
    include_list   = lin_df["include"].values.flatten().tolist()[0].split(",")
    exclude_list   = lin_df["exclude"].values.flatten().tolist()[0].split(",")
    fixed_tol      = int(lin_df["fix_tol"].values[0])
    include_tol    = int(lin_df["incl_tol"].values[0])
    exclude_tol    = int(lin_df["excl_tol"].values[0])
    del(lin_df)
    print("Selecting matrices targetting fixed species list")
    matrix_list    = filter_matrix_list(og_matrix_list,fixed_list,"fixed",fixed_tol)
    print(str(len(matrix_list))+" matrices kept")
    print("\n")
    print("Selecting matrices targetting included species list")
    matrix_list    = filter_matrix_list(matrix_list,include_list,"include",include_tol)
    print(str(len(matrix_list))+" matrices kept")
    print("\n")
    print("Removing matrices targetting exluded species list")
    matrix_list    = filter_matrix_list(matrix_list,exclude_list,"exclude",exclude_tol)
    print(str(len(matrix_list))+" matrices kept")
    print("\n")
    lin_og_code_df = og_codes_df.copy()
    lin_enc_og_df  = encoded_og_df.copy()
    lin_og_code_df = lin_og_code_df[lin_og_code_df["matrix"].isin(matrix_list)]
    lin_code_list  = lin_og_code_df["code"].values.flatten().tolist()
    lin_enc_og_df  = lin_enc_og_df[lin_enc_og_df["Total"].isin(lin_code_list)]
    lin_og_list    = lin_enc_og_df["Orthogroup"].values.flatten().tolist()
    del(lin_og_code_df)
    del(lin_enc_og_df)
    for species in fixed_list:
        species_og_df = full_og_df.copy()
        species_og_df = species_og_df[species_og_df["Orthogroup"].isin(lin_og_list)]
        sp_tmp_list   = species_og_df[species].values.flatten().tolist()
        sp_gene_list  = ",".join(sp_tmp_list).split(",")
        list_file = species + "_" + lineage + ".idlist"
        with open(list_file,'w') as f:
            for sp_gene in sp_gene_list:
                f.write(f"{sp_gene}\n")
        del(species_og_df)