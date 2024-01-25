#!/bin/bash
code=$1
mat_str=$(echo "obase=2; ibase=10; ${code}" | bc | rev)
mat_len=$(echo "${mat_str}" | wc -c)
if [ "${mat_len}" -lt 19 ]
then
  pad_len=$(echo "${mat_len}" | awk '{print 19-$1}')
  pad_cmd=$(echo "printf '0%.0s' {1..${pad_len}}")
  pad_str=$(eval ${pad_cmd})
else
  pad_str=""
fi
full_mat_str="${mat_str}${pad_str}"
echo "${full_mat_str}"
