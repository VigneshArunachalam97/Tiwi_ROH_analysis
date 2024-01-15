#!/bin/bash

######### Step 1 - RUN ROH using plink #######
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

######## Step 2 - Extract a individual sample ROH from plink output #####
