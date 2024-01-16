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


# Karyogram plot for genome-wide ROH 
The Karyogram plot was adapted from the **github page - https://github.com/tkn132/roh_karyogram**

# Association between regional-wide ROH and phenotypes
We have run association analysis using two different software including plink and GCTA. Before running the association analysis, we have removed the SNPs status frequency less than 0.01 or 1% from the further association analysis. The following parameters are used to the glm using **plink**.
                
        plink --pfile /../1_ROH_analysis_plink/Common_variant/1_Association_files/2_plink_files/Tiwi_roh_plink2Out
              --covar /../3_IBD/MMAP_files/Tiwi_cov_assoc.tsv
              --covar-name age sex
              --covar-variance-standardize  
              --extract /../1_ROH_analysis_plink/Common_variant/1_Association_files/2_plink_files/maf_filtered_snps.snplist
              --glm
              --out /../Common_variant/1_Association_files/3_Association_results/1_withSTD/2_001MAF/Assoc_ROH_001        
              --pheno /../Tiwi_phenotype_assoc_std.tsv

The **GCTA script** for all the available phenotypes are given as follows 
**1. Set working directory**

            cd /home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant/1_Association_files/2_plink_files/
**2. Run GWAS using GCTA** 

           cd /mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/1_MAF0.05/
           COUNTER=1
            line="weight height waist SBP_2 DBP_2 specific_gravity urine_creatinine_1 urine_albumin_mg_dL creatinine_1 albumin_1 hba1c_1 uric_acid_1 urine_osmolality_1 eGFR_predicted ACR BMI"
            for word in $line; do echo $COUNTER; /home/n11142006/vig_Tools/gcta_1.93.3beta2/gcta64 --bfile Tiwi_roh_001maf \
                  --grm Tiwi_data_relationshp_005maf \
                  --pheno /mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Phenotype_data/Tiwi_data_n455_std.tsv \
                  --covar /mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Phenotype_data/Tiwi_data_n455_sex_info.tsv \
                  --qcovar /mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/1_MAF0.05/Tiwi_data_qcovar_2pc_info.txt \
                  --mlma \
                  --out "/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant/1_Association_files/3_Association_results/2_withSTD/2_001MAF_GCTA/Tiwi_data_roh_$word" \
                  --mpheno $COUNTER;COUNTER=$((COUNTER+1));done
**3. Create a manhattan plot for output**

                 for i in *.mlma; \
                        do echo $i; \
                        var="$(echo $i | sed 's/.mlma//g' | sed 's/Tiwi_data_roh_//g')"; \
                        echo ${var}; \
                        Rscript /home/n11142006/TIwi_data/Script/Functions/Manhattan_qq_combined.R $i Assoc_ROH_${var} "${var} - Regional ROH" "GCTA"; \
                done
**4. Run annotation for the output**

      Rscript /home/n11142006/TIwi_data/Script/Functions/Merge_annotation_for_plink.R "../../3_Association_results/../ ".mlma"     

Finally, Now we have completed the regional-wide ROH association analysis for the given population.
