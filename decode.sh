#!/bin/bash
code=$1
species_list_file=$2
run_mode=$3
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
echo "${full_mat_str}"
if [ "${run_mode}" == "full" ]
then
  col_str=$(echo "${full_mat_str}" | perl -pe 's/1/1\n/g;s/0/0\n/g' | grep -n 1 | cut -d\: -f1 | perl -pe 's/\n/\,/g' | perl -pe 's/\,$//')
  echo "${species_list}" | perl -pe 's/\n/\t/g' | cut -f${col_str} | perl -pe 's/\t/\n/g'
fi