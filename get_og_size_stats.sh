#!/bin/bash
tsv_file=$1
datablock=$(grep -wv 0 ${tsv_file} | grep -wv Orthogroup | rev | cut -f1 --complement | rev)
og_list=$(echo "${datablock}" | cut -f1 | sort -V | uniq)
for og in ${og_list}
do
    og_size_list=$(echo   "${datablock}" | grep -w ${og} | cut -f2- | perl -pe 's/\t/\n/g' | sort -n)
    og_size_median=$(echo "${og_size_list}" | awk '{a[i++]=$1;}END{x=int((i+1)/2);if(x<(i+1)/2){print(a[x-1]+a[x])/2}else {print a[x-1]}}')
    og_size_mean=$(echo   "${og_size_list}" | awk '{sum+=$1}END{print sum/NR}')
    og_size_mode=$(echo   "${og_size_list}" | sort -n | uniq -c | sort -n | tail -n1 | awk '{print $2}')
    echo -e "${og}\t${og_size_median}\t${og_size_mean}\t${og_size_mode}"
done