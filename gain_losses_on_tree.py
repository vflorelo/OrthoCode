#!/usr/bin/python3
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, faces
import sys
import re
#tree_file     = sys.argv[1]
#enc_og_tsv    = sys.argv[2]
#og_codes      = sys.argv[3]
#sp_list_file  = sys.argv[4]
tree_file     = "SpeciesTree.ultra.nwk"
svg_file      = "SpeciesTree.ultra.svg"
enc_og_tsv    = "Orthogroups.Encoded.tsv"
og_codes      = "og_codes.txt"
sp_list_file  = "species_list"
sp_list_fh    = open(sp_list_file, "r") 
sp_list_data  = sp_list_fh.read() 
species_list  = sp_list_data.split("\n") 
species_list  = list(filter(None, species_list))
species_count = len(species_list)
sp_list_fh.close()
tree        = Tree(tree_file)
enc_og_df   = pd.read_csv(enc_og_tsv,sep="\t",dtype=str)
og_codes_df = pd.read_csv(og_codes,  sep="\t",dtype=str)
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
counter = 0
for n in tree.traverse():
    counter += 1
    core_gain_file = "lists/node_"+str(counter)+"_core_gain.idlist"
    excl_gain_file = "lists/node_"+str(counter)+"_excl_gain.idlist"
    core_loss_file = "lists/node_"+str(counter)+"_core_loss.idlist"
    excl_loss_file = "lists/node_"+str(counter)+"_excl_loss.idlist"
    excl_gain_base_list = ["0"] * species_count
    core_gain_base_list = ["0"] * species_count
    excl_loss_base_list = ["1"] * species_count
    core_loss_base_list = ["."] * species_count
    cur_species_list    = []
    cur_index_list      = []
    if (n.is_leaf() == True):
        cur_species_list = [n.name]
    else:
        for sub_n in n.get_descendants():
            if( sub_n.is_leaf() == True):
                cur_species_list.append(sub_n.name)
    for cur_species in cur_species_list:
        cur_index = species_list.index(cur_species)
        cur_index_list.append(cur_index)
    for cur_index in cur_index_list:
        excl_gain_base_list[cur_index] = "1"
        core_gain_base_list[cur_index] = "."
        excl_loss_base_list[cur_index] = "0"
        core_loss_base_list[cur_index] = "0"
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
    with open(core_gain_file, 'w') as f:
        for gain in core_gain_list:
            f.write(f"{gain}\n")
    with open(excl_gain_file, 'w') as f:
        for gain in excl_gain_list:
            f.write(f"{gain}\n")
    with open(core_loss_file, 'w') as f:
        for loss in core_loss_list:
            f.write(f"{loss}\n")
    with open(excl_loss_file, 'w') as f:
        for loss in excl_loss_list:
            f.write(f"{loss}\n")
    core_gain        = len(enc_og_df[core_gain_index])
    excl_gain        = len(enc_og_df[excl_gain_index])
    core_loss        = len(enc_og_df[core_loss_index])
    excl_loss        = len(enc_og_df[excl_loss_index])
    n.add_feature("nid",counter)
    n.add_feature("core_gain",core_gain)
    n.add_feature("excl_gain",excl_gain)
    n.add_feature("core_loss",core_loss)
    n.add_feature("excl_loss",excl_loss)
def node_layout(sub_node):
    core_gain = str(sub_node.core_gain)
    excl_gain = str(sub_node.excl_gain)
    core_loss = str(sub_node.core_loss)
    excl_loss = str(sub_node.excl_loss)
    nid       = str(sub_node.nid)
    node_str  = nid+" (/++"+core_gain+"/+"+excl_gain+"/--"+core_loss+"/-"+excl_loss+")"
    node_text = faces.TextFace(node_str,fsize=20)
    faces.add_face_to_node(node_text,sub_node,0,position="branch-bottom")
ts = TreeStyle()
ts.layout_fn = node_layout
nstyle = NodeStyle()
for n in tree.traverse():
    n.set_style(nstyle)
#tree.show(tree_style=ts)
tree.render(svg_file,tree_style=ts)