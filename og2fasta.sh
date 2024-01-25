#!/bin/bash
prot_id=$1
folder_name=$2
fasta_file=$3
og_id=$(grep -w ${prot_id} ${folder_name}/Orthogroups.tsv | cut -f1)
prot_list=$(grep -w ${og_id} ${folder_name}/Orthogroups.tsv | cut -f1 --complement | perl -pe 's/\ /\n/g;s/\,/\n/g;s/\t/\n/g' | sort -V | uniq | grep -i ^[a-z])
seqtk subseq ${fasta_file} <(echo "${prot_list}") > ${og_id}.fasta