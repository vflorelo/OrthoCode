#!/bin/bash
og_id=$1
tsv_file=$2
grep -w ^"${og_id}" "${tsv_file}" | cut -f2- | perl -pe 's/\t/\n/g;s/\,/\n/g;s/\ /\n/g' | sort -V | uniq | grep -v ^$ | grep . | perl -pe "s/$/\t${og_id}/"