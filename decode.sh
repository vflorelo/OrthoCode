#!/bin/bash
code=$1
species_list_file=$2
species_list=$(cat ${species_list_file})
num_taxa=$(cat ${species_list_file} | sort -V | uniq | grep -v ^$ | wc -l)
mat_str=$(echo "obase=2; ibase=10; ${code}" | bc | rev | perl -pe 's/\n//')
mat_len=$(echo "${mat_str}" | awk '{print length($1)}')
if [ "${mat_len}" -lt "${num_taxa}" ]
then
  pad_len=$(echo "${mat_len}" | awk -v num_taxa="${num_taxa}" '{print num_taxa - $1}')
  pad_cmd=$(echo "printf '0%.0s' {1..${pad_len}}")
  pad_str=$(eval ${pad_cmd})
else
  pad_str=""
fi
full_mat_str="${mat_str}${pad_str}"
echo -e "${code}\t${full_mat_str}"