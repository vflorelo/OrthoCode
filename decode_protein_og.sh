#!/bin/bash
prot_id=$1
og_tsv_file=$2
og_bin_file=$3
bin_col=$(head -n1 ${og_bin_file} | perl -pe 's/\t/\n/g' | grep -n . | cut -d\: -f1 | tail -n1 )
prot_type=$(grep -w $(grep -w ${prot_id} ${og_tsv_file} | cut -f1) ${og_bin_file} | cut -f${bin_col})
echo -e "${prot_id}\t${prot_type}"