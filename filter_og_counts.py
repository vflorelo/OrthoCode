#!/usr/bin/python3
import pandas as pd
import sys
in_file      = sys.argv[1]
out_file     = sys.argv[2]
max_fam_size = int(sys.argv[3])
max_missing  = int(sys.argv[4])
def get_family_type(orthogroup):
    count_list = df[df["Orthogroup"]==orthogroup][species_list].values.flatten().tolist()
    count_len  = len(count_list)
    zero_count = count_list.count(0)
    max_count  = int(max(count_list))
    min_count  = int(min(count_list))
    fam_range  = (max_count - min_count) + 1
    if(zero_count >= max_missing):
        fam_type = "exclusive"
    else:
        fam_type = "extended"
    return fam_type, fam_range
df = pd.read_csv(in_file,sep="\t")
df = df.drop(columns=["Total"])
species_list = list(df.columns)[1:]
big_fam_list = df.index[(df[species_list] > max_fam_size).any(axis=1)].tolist()
df = df.drop(index=big_fam_list)
df[["fam_type","fam_range"]] = df.apply(lambda x: get_family_type(x["Orthogroup"]),axis=1,result_type='expand')
exc_fam_list = df.index[df["fam_type"]=="exclusive"].tolist()
df = df.drop(index=exc_fam_list)
wild_fam_list = df.index[df["fam_range"]>=15].tolist()
df = df.drop(index=wild_fam_list)
df = df.drop(columns=["fam_type","fam_range"])
df.to_csv(out_file,sep="\t",index=False)
