#!/usr/bin/python3
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, faces
import sys
import re
tree_file    = sys.argv[1]
og_counts    = sys.argv[2]
og_codes_txt = sys.argv[3]
sp_list_file = sys.argv[4]
sp_list_fh   = open(sp_list_file, "r") 
sp_list_data = sp_list_fh.read() 
species_list = sp_list_data.split("\n") 
sp_list_fh.close() 
#code   matrix
#1      100000000000000000
#2      010000000000000000
#species_list = ["Pfalciparum","Preichenowi","Pgaboni","Pberghei","Pyoelii","Pchabaudi","Hepatocystis","Pknowlesi","Pvivax","Pgallinaceum","Htartakovskyi","Tgondii","Ncaninum","Sneurona","Etenella","Ccayetanensis","Tannulata","Bbovis"]
tree         = Tree(tree_file)
og_counts_df = pd.read_csv(og_counts,sep="\t",dtype=str)
og_codes_df  = pd.read_csv(og_codes_txt,sep="\t",dtype=str)
def get_matrix_code(matrix_str):
    matrix_list = list(matrix_str)
    one_count = 0
    for matrix_element in range(len(matrix_list)):
        if matrix_list[matrix_element] == "1":
            one_count += 1
    eff_matrix_len = len(matrix_list) - one_count
    tmp_df = og_codes_df.copy()
    tmp_df["match"]     = tmp_df["matrix"].apply(lambda x: re.findall(matrix_str,x)[0] if re.findall(matrix_str,x) else "")
    tmp_index_list = tmp_df["match"].apply(lambda x: True if x else False)
    tmp_df = tmp_df[tmp_index_list]
    tmp_df["match_sum"] = tmp_df["match"].apply(lambda x: sum(list(map(int,list(x)))))
    tmp_df["match_frac"] = (tmp_df["match_sum"] - one_count)/(eff_matrix_len)
    tmp_df = tmp_df[tmp_df["match_frac"]>=0.8]
    tmp_code_list = tmp_df["code"].values.flatten().tolist()
    return tmp_code_list
diff_list = []
counter = 0
for n in tree.traverse():
    counter += 1
    gain_base_list = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    loss_base_list = [".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".","."]
    cur_species_list = []
    cur_index_list = []
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
        gain_base_list[cur_index] = 1
        loss_base_list[cur_index] = 0
    gain_matrix_str = ''.join(str(list_item) for list_item in gain_base_list)
    loss_matrix_str = ''.join(str(list_item) for list_item in loss_base_list)
    gain_matrix_str = gain_matrix_str[::-1]
    node_gain_code  = int(gain_matrix_str,base=2)
    node_loss_code  = get_matrix_code(loss_matrix_str)
    node_gain_og_index = og_counts_df["Total"] == str(node_gain_code)
    node_loss_og_index = og_counts_df["Total"].isin(node_loss_code)
    node_gains  = len(og_counts_df[node_gain_og_index])
    node_losses = len(og_counts_df[node_loss_og_index])
    node_diff   = node_gains+node_losses
    n.add_feature("nid",counter)
    n.add_feature("gains",node_gains)
    n.add_feature("losses",node_losses)
    n.add_feature("diff",node_diff)
def pie_layout(mynode):
    gains     = mynode.gains
    losses    = mynode.losses
    diff      = mynode.diff
    nid       = mynode.nid
    gain_frac = (gains/diff)*100
    loss_frac = (losses/diff)*100
    node_data = [gain_frac,loss_frac]
    node_cols = ["blue","red"]
    piechart  = faces.PieChartFace(node_data,colors=node_cols,width=50,height=50)
    faces.add_face_to_node(piechart,mynode,0,position="branch-top")
    diff_text_str = str(nid)+"(+"+str(gains)+"/-"+str(losses)+")"
    diff_text = faces.TextFace(diff_text_str,fsize=20)
    faces.add_face_to_node(diff_text,mynode,0,position="branch-bottom")
ts = TreeStyle()
ts.layout_fn = pie_layout
nstyle = NodeStyle()
for n in tree.traverse():
    n.set_style(nstyle)
tree.render("test_tree.svg",tree_style=ts)