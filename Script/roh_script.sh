#!/bin/bash

#### Step 1 - RUN ROH using plink ####
file="Tiwi_data"
out="tiwi_roh_out"
plink --bfile ${file} \
      --homozyg \
      --homozyg-density 50 \
      --homozyg-gap 1000 \
      --homozyg-kb 1500 \
      --homozyg-snp 50 \
      --homozyg-window-het 1 \
      --homozyg-window-missing 5 \
      --homozyg-window-snp 50 \
      --keep-allele-order \
      --out ${out}

#### Step 2 - Extract a individual sample ROH from plink output ####
Rscript ROH_plink_to_data.R ### this will create matrix for plink file and create a mmap to convert that into plink (map and bed files)

#### Step 3 - convert into mmap and then plink file ####
module load mmap/2022_11_06_21_46.intel   #### load mmap from hpc modules ####
mmap --write_binary_genotype_file --csv_input_filename /mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant/1_Association_files/Tiwi_data_regional_roh_plink.csv --binary_output_filename /mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant/1_Association_files/Tiwi_data_regional_roh_plink.bin --num_skip_fields 8

#### convert into binary file ####
cd /mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant/1_Association_files/
mmap --transpose_binary_genotype_file --binary_input_filename Tiwi_data_regional_roh_plink.bin --binary_output_filename Tiwi_data_plink_SxM  #### convert the marker by subject into subject by marker ####
mmap --subject_by_marker_mmap2plink --binary_input_filename Tiwi_data_plink_SxM --plink_output_prefix Tiwi_roh_plink  #### convert the subject by marker ####

#### use plink to create a binary file (*.bed, *.fam, *.bim files) ####
module load plink/1.9b_6.21-x86_64
plink --file Tiwi_roh_plink --make-bed --keep-allele-order --out Tiwi_roh_plinkOut
