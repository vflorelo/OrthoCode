#!/usr/bin/python3
from ete3 import Tree
import sys
tree_file = sys.argv[1]
species_1 = sys.argv[2]
species_2 = sys.argv[3]
umt_file  = sys.argv[4]
tree      = Tree(tree_file)
xg_name = tree.get_common_ancestor(species_1,species_2)
tree.set_outgroup(xg_name)
tree.convert_to_ultrametric()
print(tree.write(outfile=umt_file))
