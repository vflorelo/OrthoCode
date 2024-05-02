#!/bin/bash
usage="$(dirname $0)/usage.sh"
source ${usage}
while [ "$1" != "" ]
do
    case $1 in
        --tsv_dir         )
            shift
            tsv_dir=$(realpath $1)
            ;;
        --counts_file     )
            shift
            counts_file=$1
            ;;
        --unassigned_file )
            shift
            unassigned_file=$1
            ;;
        --species_list    )
            shift
            species_list=$1
            ;;
		--threads       )
            shift
            threads=$1
            ;;
		--help          )
            usage
            exit 0
            ;;
	esac
	shift
done
if [ ! -d "${tsv_dir}" ] && [ ! -z "${tsv_dir}" ]
then
    echo "Missing fasta directory. Exiting"
    exit 0
elif [ -z "${tsv_dir}" ]
then
    if [ -z "${counts_file}" ] || [ -z "${unassigned_file}" ]
    then
        echo "Missing required options. Exiting"
        exit 0
    fi
    if [ ! -f "${counts_file}" ] || [ ! -f "${unassigned_file}" ]
    then
        echo "Missing required files. Exiting"
        exit 0
    fi
fi

num_proc=$(nproc)
if [ -z "${threads}" ]
then
    threads=${num_proc}
else
    threads_test=$(echo -e "${threads}\t${num_proc}" | awk '{if((int($1)==$1) && ($1<=$2)){print $1}}')
    if [ -z "${threads_test}" ]
    then
        echo "Invalid threads value. Exiting"
        exit
    fi
fi


uuid=$(uuidgen | cut -d\- -f5)
cur_date=$(date +%Y-%m-%d)
log_file="OrthoCode.${cur_date}.log"
err_file="OrthoCode.${cur_date}.err"

if [ -z "${species_list}" ] || [ ! -f "${species_list}" ]
then
    echo "Warning: Missing or unspecified species list file"
    echo "Generating species list from counts file:"
    echo "species_list.${uuid}.txt"
    head -n1 "${counts_file}" | cut -f2- | perl -pe 's/\t/\n/g' | head -n-1 > species_list.${uuid}.txt
    species_list="species_list.${uuid}.txt"
fi

echo "Creating full counts file"
orthocounts.py "${counts_file}" "${unassigned_file}" > "${log_file}" 2> "${err_file}"
exit_code=$?
if [ "${exit_code}" -gt 0 ]
then
    echo "Something went wrong with orthocounts.py"
    echo "check ${log_file} and ${err_file} for more details"
    exit 0
else
    echo "Full gene counts file created (Orthogroups.FullGeneCount.tsv)"
    echo "Presence/absence matrix created (Orthogroups.Presence.tsv)"
    full_counts_file="Orthogroups.FullGeneCount.tsv"
    presence_file="Orthogroups.Presence.tsv"
fi

echo "Encoding Orthogroups table"
encode_orthogroups.py "${presence_file}" "${species_list}" >> ${log_file} 2>> ${err_file}
exit_code=$?
if [ "${exit_code}" -gt 0 ]
then
    echo "Something went wrong with encode_orthogroups.py"
    echo "check ${log_file} and ${err_file} for more details"
    exit 0
else
    echo "Encoded orthogroups table created (Orthogroups.Encoded.tsv)"
    encoded_file="Orthogroups.Encoded.tsv"
fi

echo "Creating orthogroup code index"
code_list=$(tail -n+2 "${encoded_file}" | sort -n | uniq )
echo "${code_list}" | awk -v species_list="${species_list}" '{print "decode.sh",$1,species_list}' | parallel -j ${threads} > codes.tmp.tsv 2>> ${err_file}
exit_code=$?
if [ "${exit_code}" -gt 0 ]
then
    echo "Orthogroup code indexing failed, see ${err_file} for details"
    exit 0
else
    echo -e "code\tmatrix" > og_codes.tsv
    sort -n codes.tmp.tsv >> og_codes.tsv
    echo "Codes and matrices created in og_codes.tsv"
fi
rm codes.tmp.tsv