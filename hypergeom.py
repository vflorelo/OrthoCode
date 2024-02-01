#!/usr/bin/python3
import sys
import scipy.stats as stats
global_prot_count = int(sys.argv[1])
global_comp_count = int(sys.argv[2])
cur_prot_count    = int(sys.argv[3])
cur_comp_count    = int(sys.argv[4])
stats.hypergeom.cdf(cur_comp_count,global_prot_count,global_comp_count,cur_prot_count)
print (str(stats.hypergeom.sf(cur_comp_count-1,global_prot_count,global_comp_count,cur_prot_count)))