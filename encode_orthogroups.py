#!/usr/bin/python3
import pandas as pd
import sys
mat_tsv_file = sys.argv[1]
out_tsv_file = sys.argv[2]
mat_df       = pd.read_csv(mat_tsv_file,sep="\t").fillna("")
species_list = ['Plasmodium_falciparum_3D7',
                'Plasmodium_reichenowi_CDC',
                'Plasmodium_gaboni_SY75',
                'Plasmodium_berghei_ANKA',
                'Plasmodium_yoelii_yoelii_17X',
                'Plasmodium_chabaudi_chabaudi',
                'Hepatocystis_2019',
                'Plasmodium_knowlesi_H',
                'Plasmodium_vivax_P01',
                'Plasmodium_gallinaceum_8A',
                'Haemoproteus_tartakovskyi_SISKIN1',
                'Toxoplasma_gondii_ME49',
                'Neospora_caninum_Liverpool2019',
                'Sarcocystis_neurona_SN3',
                'Eimeria_tenella_Houghton2021',
                'Cyclospora_cayetanensis_CHN_HEN01',
                'Theileria_annulata_Ankara',
                'Babesia_bovis_T2Bo']
counter=0
for species in species_list:
    placement = 2**counter
    counter  += 1
    mat_df[species] = mat_df[species] * placement
mat_df["Total"] = mat_df[species_list].sum(axis=1)
mat_df.to_csv(out_tsv_file,sep="\t",index=False)