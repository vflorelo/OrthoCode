#!/bin/bash
tree_file=$1
base_name=$(echo ${tree_file} | rev | cut -d\. -f1 --complement | rev)
perl -pe 's/\(/\n\(/g;s/\)/\n\)/g;s/\,/\n\,/g;s/\:/\n\:/g;s/\).*\//\)/g;s/\n//g' ${tree_file} > ${base_name}.boot.nwk
perl -pe 's/\(/\n\(/g;s/\)/\n\)/g;s/\,/\n\,/g;s/\:/\n\:/g;s/\/.*/\)/g;s/\n//g'   ${tree_file} > ${base_name}.alrt.nwk