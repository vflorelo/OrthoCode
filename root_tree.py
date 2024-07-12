#!/usr/bin/python3
from ete3 import Tree
import sys
tree_file     = sys.argv[1]
root_species  = sys.argv[2]
out_tree_file = sys.argv[3]
tree          = Tree(tree_file)
tree.set_outgroup(root_species)
print(tree.write(outfile=out_tree_file))
