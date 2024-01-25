#!/bin/bash
prot_code=$1
species_list=$(echo -e "Plasmodium_falciparum_3D7\nPlasmodium_reichenowi_CDC\nPlasmodium_gaboni_SY75\nPlasmodium_berghei_ANKA\nPlasmodium_yoelii_yoelii_17X\nPlasmodium_chabaudi_chabaudi\nHepatocystis_2019\nPlasmodium_knowlesi_H\nPlasmodium_vivax_P01\nPlasmodium_gallinaceum_8A\nHaemoproteus_tartakovskyi_SISKIN1\nToxoplasma_gondii_ME49\nNeospora_caninum_Liverpool2019\nSarcocystis_neurona_SN3\nEimeria_tenella_Houghton2021\nCyclospora_cayetanensis_CHN_HEN01\nTheileria_annulata_Ankara\nBabesia_bovis_T2Bo")
species_list_len=$(echo "${species_list}" | wc -l)
mat_str=$(echo "obase=2; ibase=10; ${prot_code}" | bc | rev)
mat_len=$(echo "${mat_str}" | wc -c)
pad_len=$(echo -e "${species_list_len}\t${mat_len}" | awk '{print ($1-$2)+1}')
pad_cmd=$(echo "printf '0%.0s' {1..${pad_len}}")
pad_str=$(eval ${pad_cmd})
full_mat_str="${mat_str}${pad_str}"
col_list=$(echo "${full_mat_str}" | perl -pe 's/0/0\n/g;s/1/1\n/g' | grep -wn 1 | cut -d\: -f1 | perl -pe 's/\n/\,/' | perl -pe 's/\,$//')
species_list=$(echo "${species_list}" | perl -pe 's/\n/\t/g'| cut -f${col_list} | perl -pe 's/\t/\n/g')
echo "${species_list}"
