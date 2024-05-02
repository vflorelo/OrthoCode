#!/bin/bash
tsv_file=$1
head -n 1 "${tsv_file}" | perl -pe 's/\t/\n/g' | grep -n .
