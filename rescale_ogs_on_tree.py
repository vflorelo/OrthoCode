#!/usr/bin/python3
from   ete3    import Tree, TreeStyle, NodeStyle, faces
from   pathlib import Path
import pandas  as pd
import sys
import re
sp_tree_file  = "species_tree.ultra.nwk"
enc_og_tsv    = "Orthogroups.Encoded.tsv"
og_codes      = "og_codes.tsv"
sp_list_file  = "species_list"
go_terms_file = "GO_master_file.tsv"
sp_list_fh    = open(sp_list_file, "r") 
sp_list_data  = sp_list_fh.read() 
species_list  = sp_list_data.split("\n") 
species_list  = list(filter(None, species_list))
species_count = len(species_list)
root_species  = "Cryptosporidium_parvum"
sp_list_fh.close()
species_tree  = Tree(sp_tree_file)
species_tree.set_outgroup(root_species)
species_tree.convert_to_ultrametric()
species_tree.ladderize()
enc_og_df   = pd.read_csv(enc_og_tsv,sep="\t",dtype=str)
og_codes_df = pd.read_csv(og_codes,  sep="\t",dtype=str)
go_terms_df = pd.read_csv(go_terms_file,sep="\t",dtype=str)
def get_matrix_code(matrix_str,matrix_type):
    matrix_list = list(matrix_str)
    dot_count = 0
    one_count = 0
    for matrix_element in range(len(matrix_list)):
        if matrix_list[matrix_element] == ".":
            dot_count += 1
        if matrix_list[matrix_element] == "1":
            one_count += 1
    loss_tolerance = dot_count - 1
    gain_tolerance = one_count - 1
    if (matrix_type == "gain"):
        tolerance = gain_tolerance
    elif(matrix_type == "loss"):
        tolerance = loss_tolerance
    mat_df              = og_codes_df.copy()
    mat_df["match"]     = mat_df["matrix"].apply(lambda x: re.findall(matrix_str,x)[0] if re.findall(matrix_str,x) else "")
    index_list          = mat_df["match"].apply(lambda x: True if x else False)
    mat_df              = mat_df[index_list]
    mat_df["match_sum"] = mat_df["match"].apply(lambda x: sum(list(map(int,list(x)))))
    mat_df              = mat_df[mat_df["match_sum"]>=tolerance]
    code_list           = mat_df["code"].values.flatten().tolist()
    return code_list
counter         = 0
node_list       = []
core_gains      = []
core_gain_count = []
core_losses     = []
core_loss_count = []
excl_gains      = []
excl_gain_count = []
excl_losses     = []
excl_loss_count = []
parent          = []
parent_count    = []
for node in species_tree.traverse():
    counter += 1
    node.add_feature("node_id",counter)
    excl_gain_base_list = ["0"] * species_count
    core_gain_base_list = ["0"] * species_count
    excl_loss_base_list = ["1"] * species_count
    core_loss_base_list = ["."] * species_count
    cur_species_list    = []
    cur_index_list      = []
    if (node.is_leaf() == True):
        cur_species_list = [node.name]
    else:
        for sub_node in node.get_descendants():
            if( sub_node.is_leaf() == True):
                cur_species_list.append(sub_node.name)
    for cur_species in cur_species_list:
        cur_index = species_list.index(cur_species)
        cur_index_list.append(cur_index)
    for cur_index in cur_index_list:
        excl_gain_base_list[cur_index] = "1"
        core_gain_base_list[cur_index] = "."
        excl_loss_base_list[cur_index] = "0"
        core_loss_base_list[cur_index] = "0"
    if (node.is_root()==True):
        parent_id = ""
    else:
        parent_id = node.up.node_id
    print(cur_species_list)
    core_gain_matrix = ''.join(list_item for list_item in core_gain_base_list)
    excl_gain_matrix = ''.join(list_item for list_item in excl_gain_base_list)
    excl_gain_matrix = excl_gain_matrix[::-1]
    core_loss_matrix = ''.join(list_item for list_item in core_loss_base_list)
    excl_loss_matrix = ''.join(list_item for list_item in excl_loss_base_list)
    excl_loss_matrix = excl_loss_matrix[::-1]
    core_gain_code   = get_matrix_code(core_gain_matrix,"gain")
    excl_gain_code   = int(excl_gain_matrix,base=2)
    core_loss_code   = get_matrix_code(core_loss_matrix,"loss")
    excl_loss_code   = int(excl_loss_matrix,base=2)
    core_gain_index  = enc_og_df["Total"].isin(core_gain_code)
    excl_gain_index  = enc_og_df["Total"] == str(excl_gain_code)
    core_loss_index  = enc_og_df["Total"].isin(core_loss_code)
    excl_loss_index  = enc_og_df["Total"] == str(excl_loss_code)
    core_gain_list   = enc_og_df[core_gain_index]["Orthogroup"].values.flatten().tolist()
    excl_gain_list   = enc_og_df[excl_gain_index]["Orthogroup"].values.flatten().tolist()
    core_loss_list   = enc_og_df[core_loss_index]["Orthogroup"].values.flatten().tolist()
    excl_loss_list   = enc_og_df[excl_loss_index]["Orthogroup"].values.flatten().tolist()
    core_gain = len(core_gain_list)
    excl_gain = len(excl_gain_list)
    core_loss = len(core_loss_list)
    excl_loss = len(excl_loss_list)
    node_list.append(counter)
    core_gains.append(core_gain_list)
    core_losses.append(core_loss_list)
    excl_gains.append(excl_gain_list)
    excl_losses.append(excl_loss_list)
    parent.append(parent_id)
    core_gain_count.append(core_gain)
    core_loss_count.append(core_loss)
    excl_gain_count.append(excl_gain)
    excl_loss_count.append(excl_loss)
gains_dict = {"node"           : node_list,
              "core_gains"     : core_gains,
              "core_losses"    : core_losses,
              "excl_gains"     : excl_gains,
              "excl_losses"    : excl_losses,
              "core_gain_count": core_gain_count,
              "excl_gain_count": excl_gain_count,
              "core_loss_count": core_loss_count,
              "excl_loss_count": excl_loss_count,
              "parent":parent}
gains_df = pd.DataFrame(gains_dict)
def get_parent_ogs(node_id):
    parent_id      = gains_df[gains_df["node"]==node_id]["parent"].values.flatten().tolist()[0]
    if (parent_id!=""):
        parent_og_list = gains_df[gains_df["node"]==parent_id]["core_gains"].values.flatten().tolist()[0]
    else:
        parent_og_list = gains_df[gains_df["node"]==node_id]["core_gains"].values.flatten().tolist()[0]
    return parent_og_list
gains_df["parent_og_list"] = gains_df["node"].apply(lambda x: get_parent_ogs(x))
go_types=["BP","CC","MF"]
node_list = gains_df["node"]
# event_types = ["core_gains","core_losses","excl_gains","excl_losses"]
# for node in node_list:
#     node_directory = "node_"+str(node)
#     parent_og_list = gains_df[gains_df["node"]==node]["parent_og_list"].values.flatten().tolist()[0]
#     Path(node_directory).mkdir(parents=True, exist_ok=True) 
#     for go_type in go_types:
#         bg_df = go_terms_df.copy()
#         bg_df = bg_df[(bg_df["orthogroup"].isin(parent_og_list))&(bg_df["go_type"]==go_type)][["orthogroup","go_id"]]
#         bg_df = bg_df.drop_duplicates()
#         bg_file_name = node_directory+"/bg."+str(go_type)+".txt"
#         bg_df.to_csv(bg_file_name,header=False,index=False,sep="=")
#         for event_type in event_types:
#             og_list = gains_df[gains_df["node"]==node][event_type].values.flatten().tolist()[0]
#             tmp_df  = go_terms_df.copy()
#             tmp_df  = tmp_df[(tmp_df["orthogroup"].isin(og_list))&(tmp_df["go_type"]==go_type)][["orthogroup","go_id"]]
#             tmp_file_name = node_directory+"/"+str(event_type)+"."+str(go_type)+".txt"
#             tmp_df.to_csv(tmp_file_name,header=False,index=False,sep="=")