# Regional-ROH-association-with-phenotype
This repository provides the pipeline developed to do ROH analysis for the Tiwi indigenous community and their genome wide ROH association with the phenotypes related to kidney individual

Initially, we have used to plink pipeline to find the ROH region in the entire Tiwi population (N = 455). Here, we have used only the common SNPs (MAF > 0.05 or 5%)

**Step 1:** Run a ROH using the "homozyg" function in the plink v1.9. The following plink parameters are used to run the ROH analysis 
        
        plink --bfile Tiwi_data \
              --homozyg \
              --homozyg-density 50 \
              --homozyg-gap 1000 \
              --homozyg-kb 1500 \
              --homozyg-snp 50 \
              --homozyg-window-het 1 \
              --homozyg-window-missing 5 \
              --homozyg-window-snp 50 \
              --keep-allele-order \
              --out tiwi_data_out

Please refer the plink tutorial or website to understand the description for each parameter used (refer: https://www.cog-genomics.org/plink/1.9/ibd)

Plink produces three different output, they are 
                    1. tiwi_data_out.hom.indiv (Individuals - total region of ROH and number of segments of ROH present)
                    2. tiwi_data_out.hom.summary (Each each SNPs wise ROH information by each individuals)
                    3. tiwi_data_out.hom (Overal ROH region - start and end for each individual, length, percentage of homozygous, percentage of heterozygous)

**Step 2:** Use the tiwi_data_out.hom to extract the ROH regions for each samples. Then use bedtools to get the ROH status for each SNPs (present in the ROH/absent in the ROH region) for each individuals. Create a master data.frame for the all the individuals together for 4.9m SNPs. Its really a time consuming task for all the individuals and 4.9 million SNPs. 

**Step 3:** Now the data contains zeros and ones are being converted into plink binary data (.bed, .fam, .bim)

**Step 4:** Calculate the allele frequency and filter maf < 0.01 (1%) for further association analysis

**Step 5:** Use either GCTA or plink to run the association analysis by adding covariates including age, and gender. (population structure PC's 1 and 2, if available)

Please refer the script folder for corresponding plot and analysis.

# Visualization plot - Manhattan plot for ROH based on the SNP status ###
This below script will create manhattan plot based on the SNP frequency over the population. This is automated script and produce a png, pdf, RDS for image object as a ouput in the folder given

        Rscript manhatton_roh.R "/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant" "tiwi_data_out.hom.summary" "tiwi_data_out.hom.indiv" "tiwi_data"

