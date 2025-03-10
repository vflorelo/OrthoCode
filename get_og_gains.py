#!/usr/bin/python3
import argparse
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, faces
import re
import numpy as np
import os
parser = argparse.ArgumentParser(prog="get_og_gains.py",
                                 description="A program for calculating the number of ancestral-, core- and exclusive-gains of orthogroups guided by the topology of a species tree",
                                 epilog="""LICENSE:
                                 Distributed under the GNU General Public License (GPLv3). See License.md
                                 https://github.com/vflorelo/OrthoCode""")
parser.add_argument('-t',
                    '--species_tree',
                    type=argparse.FileType('r'),
                    required=True,
                    help='A plain text file containing a species tree')
parser.add_argument('-r',
                    '--root_species',
                    type=str,
                    required=True,
                    help='The species name to use for rooting the tree')
parser.add_argument('-e',
                    '--encoded_ogs',
                    type=argparse.FileType('r'),
                    required=True,
                    help='A tsv file containing the encoded distribution of orthogroups (see encode_orthogroups.py)' )
parser.add_argument('-c',
                    '--og_codes',
                    type=argparse.FileType('r'),
                    required=True,
                    help='A tsv file containing the presence/absence matrices for each encoded orthogroup type')
parser.add_argument('-l',
                    '--species_list',
                    type=argparse.FileType('r'),
                    required=True,
                    help='A plain-text file with the list of species to analyse, it should be in the same order as the presence/absence matrix')
parser.add_argument('-g',
                    '--orthogroups',
                    type=argparse.FileType('r'),
                    required=True,
                    help="A tsv file containing the orhtogroups from which the gain lists will be built.")


args         = parser.parse_args()
species_tree = args.species_tree
root_species = args.root_species
encoded_ogs  = args.encoded_ogs
og_codes     = args.og_codes
species_list = args.species_list
orthogroups  = args.orthogroups


if (os.path.isdir("gain_lists") == False):
    os.mkdir("gain_lists")
    os.mkdir("gain_lists/ancestral")
    os.mkdir("gain_lists/core")
    os.mkdir("gain_lists/exclusive")
else:
    if (os.path.isdir("gain_lists/ancestral") == False):
        os.mkdir("gain_lists/ancestral")
    if (os.path.isdir("gain_lists/core") == False):
        os.mkdir("gain_lists/core")
    if (os.path.isdir("gain_lists/exclusive") == False):
        os.mkdir("gain_lists/exclusive")

out_svg_file      = "annotated_tree.svg"

species_list_fh   = open(species_list, "r") 
species_list_data = species_list_fh.read() 
species_list      = species_list_data.split("\n") 
species_list      = list(filter(None, species_list))
species_count     = len(species_list)
species_list_fh.close()



tree = Tree(species_tree)
tree.set_outgroup(root_species)
tree.convert_to_ultrametric()
tree.ladderize()


encoded_ogs_df = pd.read_csv(encoded_ogs,sep="\t",dtype=str)
og_codes_df    = pd.read_csv(og_codes,  sep="\t",dtype=str)
orthogroups_df = pd.read_csv(orthogroups,sep="\t",dtype=str)
counter = 0

def get_matrix_codes(matrix_str,matrix_type,a_index,b_index):
    a_len = len(a_index)
    b_len = len(b_index)
    matrix_df = og_codes_df.copy()
    if(matrix_type == "exclusive"):
        matrix_df = matrix_df[matrix_df["matrix"]==matrix_str]
        code_list = matrix_df["code"].values.flatten().tolist()
    else:
        if (matrix_type == "core" ):
            a_min = int(a_len/2) + 1
            b_min = int(b_len/2) + 1
        elif(matrix_type == "ancestral"):
            a_min = 1
            b_min = 1
        matrix_df["match"]     = matrix_df["matrix"].apply(lambda x: list(map(int,re.findall(matrix_str,x)[0])) if re.findall(matrix_str,x) else [])
        index_list          = matrix_df["match"].apply(lambda x: True if x else False)
        matrix_df              = matrix_df[index_list]
        matrix_df["ch_a_sum"]  = matrix_df["match"].apply(lambda x: sum(list(np.array(x)[a_index])) if x else 0)
        matrix_df["ch_b_sum"]  = matrix_df["match"].apply(lambda x: sum(list(np.array(x)[b_index])) if x else 0)
        matrix_df              = matrix_df[(matrix_df["ch_a_sum"] >= a_min ) & (matrix_df["ch_b_sum"] >= b_min)]
        code_list           = matrix_df["code"].values.flatten().tolist()
    return code_list
for tree_node in tree.traverse():
    counter += 1
    exclusive_base_list = ["0"] * species_count
    core_base_list      = ["0"] * species_count
    ancestral_base_list = ["0"] * species_count
    cur_sp_list    = []
    if (tree_node.is_leaf() == True):
        cur_sp_list = [tree_node.name]
    else:
        for sub_node in tree_node.get_descendants():
            if( sub_node.is_leaf() == True):
                cur_sp_list.append(sub_node.name)
    if(len(cur_sp_list)<4):
        for cur_sp in cur_sp_list:
            sp_index = species_list.index(cur_sp)
            exclusive_base_list[sp_index] = "1"
            core_base_list[sp_index] = "1"
            ancestral_base_list[sp_index]  = "1"
        a_index     = []
        b_index     = []
        exclusive_matrix = ''.join(exclusive_base_list)
        core_matrix      = ''.join(core_base_list)
        ancestral_matrix = ''.join(ancestral_base_list)
        exclusive_codes  = get_matrix_codes(exclusive_matrix,"exclusive",a_index,b_index)
        core_codes       = get_matrix_codes(exclusive_matrix,"exclusive",a_index,b_index)
        ancestral_codes  = get_matrix_codes(exclusive_matrix,"exclusive",a_index,b_index)
    else:
        children = tree_node.get_children()
        child_a = children[0]
        child_b = children[1]
        ch_a_sp_list = []
        ch_b_sp_list = []
        a_index      = []
        b_index      = []
        if (child_a.is_leaf() == True):
            ch_a_sp_list.append(child_a.name)
            sp_index = species_list.index(child_a.name)
            a_index.append(sp_index)
        else:
            for sub_ch_node in child_a.get_descendants():
                if(sub_ch_node.is_leaf() == True):
                    ch_a_sp_list.append(sub_ch_node.name)
                    sp_index = species_list.index(sub_ch_node.name)
                    a_index.append(sp_index)
        if (child_b.is_leaf() == True):
            ch_b_sp_list.append(child_b.name)
            sp_index = species_list.index(child_b.name)
            b_index.append(sp_index)
        else:
            for sub_ch_node in child_b.get_descendants():
                if(sub_ch_node.is_leaf() == True):
                    ch_b_sp_list.append(sub_ch_node.name)
                    sp_index = species_list.index(sub_ch_node.name)
                    b_index.append(sp_index)
        if  (len(ch_a_sp_list) == 1):
            outgroup = ch_a_sp_list[0]
        elif(len(ch_b_sp_list) == 1):
            outgroup = ch_b_sp_list[0]
        else:
            outgroup = ""
        if(outgroup == ""):
            for ch_sp in ch_a_sp_list:
                sp_index = species_list.index(ch_sp)
                exclusive_base_list[sp_index] = "1"
                core_base_list[sp_index] = "."
                ancestral_base_list[sp_index]  = "."
            for ch_sp in ch_b_sp_list:
                sp_index = species_list.index(ch_sp)
                exclusive_base_list[sp_index] = "1"
                core_base_list[sp_index] = "."
                ancestral_base_list[sp_index]  = "."
        else:
            outgroup_index = species_list.index(outgroup)
            exclusive_base_list[outgroup_index] = "1"
            core_base_list[outgroup_index] = "1"
            ancestral_base_list[outgroup_index]  = "1"
            if  (outgroup in ch_a_sp_list):
                for ch_sp in ch_b_sp_list:
                    sp_index = species_list.index(ch_sp)
                    exclusive_base_list[sp_index] = "1"
                    core_base_list[sp_index] = "."
                    ancestral_base_list[sp_index]  = "."
            elif(outgroup in ch_b_sp_list):
                for ch_sp in ch_a_sp_list:
                    sp_index = species_list.index(ch_sp)
                    exclusive_base_list[sp_index] = "1"
                    core_base_list[sp_index] = "."
                    ancestral_base_list[sp_index]  = "."
        exclusive_matrix = ''.join(exclusive_base_list)
        core_matrix = ''.join(core_base_list)
        ancestral_matrix = ''.join(ancestral_base_list)
        exclusive_codes  = get_matrix_codes(exclusive_matrix,"exclusive",a_index,b_index)
        core_codes       = get_matrix_codes(core_matrix,"core",a_index,b_index)
        ancestral_codes  = get_matrix_codes(ancestral_matrix,"ancestral",a_index,b_index)
    exclusive_index = encoded_ogs_df["Total"].isin(exclusive_codes)
    core_index      = encoded_ogs_df["Total"].isin(core_codes)
    ancestral_index = encoded_ogs_df["Total"].isin(ancestral_codes)
    exclusive_ogs   = encoded_ogs_df[exclusive_index]["Orthogroup"]
    core_ogs        = encoded_ogs_df[core_index]["Orthogroup"]
    ancestral_ogs   = encoded_ogs_df[ancestral_index]["Orthogroup"]
    exclusive_out_file = "gain_lists/exclusive/node_"+str(counter)+".og_list"
    core_out_file      = "gain_lists/core/node_"+str(counter)+".og_list"
    ancestral_out_file = "gain_lists/ancestral/node_"+str(counter)+".og_list"
    exclusive_ogs.to_csv(exclusive_out_file,sep="\t",index=None,header=None)
    core_ogs.to_csv(core_out_file,sep="\t",index=None,header=None)
    ancestral_ogs.to_csv(ancestral_out_file,sep="\t",index=None,header=None)
    exclusive_gains = len(exclusive_ogs)
    core_gains      = len(core_ogs)
    ancestral_gains = len(ancestral_ogs)
    tree_node.add_feature("nid",counter)
    tree_node.add_feature("exclusive_gains",exclusive_gains)
    tree_node.add_feature("core_gains",core_gains)
    tree_node.add_feature("ancestral_gains", ancestral_gains)
def node_layout(sub_node):
    exclusive_gains = str(sub_node.exclusive_gains)
    core_gains = str(sub_node.core_gains)
    ancestral_gains  = str(sub_node.ancestral_gains)
    nid        = str(sub_node.nid)
    node_str   = nid+" (+++"+ancestral_gains+"/++"+core_gains+"/+"+exclusive_gains+")"
    node_text  = faces.TextFace(node_str,fsize=20)
    faces.add_face_to_node(node_text,sub_node,0,position="branch-bottom")
ts = TreeStyle()
ts.layout_fn = node_layout
ts.force_topology = True
nstyle = NodeStyle()
for tree_node in tree.traverse():
    tree_node.set_style(nstyle)
tree.render(out_svg_file,tree_style=ts)