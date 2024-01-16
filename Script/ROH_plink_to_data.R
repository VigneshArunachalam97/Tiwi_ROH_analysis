#!/bin/bash
#### import required packages ####
library(dplyr)
library(data.table) 
setwd("/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant") ### set working directory ###
#### ROH plink out ######
roh <- as.data.frame(fread("/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant/ROH_common_run0.hom"))
sample_id <- unique(roh$FID)
print(length(sample_id))
#### subset the required column and remove the remaining columns ####
roh <- roh[,c("FID", "CHR", "POS1", "POS2")]
print(dim(roh))
#### import the snp for bedtools intersection ####
snp <- "/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/8_Garlic/tiwi_snps"
SNP <- as.data.frame(fread("/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/8_Garlic/tiwi_snps", header = F))
SNPs <- SNP$V5
gc()
print("SNP data imported and extracted")
##### subset use bedtools to intersect the file ####
out1 <- "/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant/1_Association_files"
samples <- list()
print("Loop started")
for(sample in sample_id){
  tmp <- subset(roh, FID == sample)
	tmp <- tmp[,c("CHR", "POS1", "POS2")]
	tmp$CHR <- paste0("chr", tmp$CHR)
	write.table(tmp, paste0(out1, "/", sample, ".tsv"), row.names = F, quote = F, col.names = F, sep = "\t")
	cmd <- paste0("bedtools intersect -a ", snp, " -b ", paste0(out1, "/", sample, ".tsv"), "> ", paste0(out1, "/", sample, "_roh"))
	system(cmd)
  unlink(paste0(out1, "/", sample, ".tsv"))
	tmp <- as.data.frame(fread(paste0(out1, "/", sample, "_roh"), header = F))
	old <- data.frame(SNP = tmp$V5, ROH_status = 1)
	new <- data.frame(SNP = setdiff(SNPs, tmp$V5), ROH_status = 0)
	old <- rbind(old, new)
  gc()
	rm(new)
	rm(tmp)
	colnames(old) <- c("SNP", sample)
	rownames(old) <- old$SNP
  gc()
	old <- old[SNPs,]
	samples[[sample]] <- old
}
gc()
print("loop ended")
snp_full <- do.call(cbind, samples)
rownames(snp_full) <- snp_full[,1]
snp_full <- snp_full[,-grep("SNP", colnames(snp_full))]
snp_full$SNP <- rownames(snp_full)
gc()
y <- ncol(snp_full)
snp_full <- snp_full[,c(y, 1:y-1)]
print("order changed")
backup_snp_full <- snp_full
gc()
colnames(snp_full) <- gsub("\\..*$", "", colnames(snp_full))
write.table(snp_full, paste0(out1, "/Tiwi_data_plink_roh_regional.tsv"), quote = F, row.names = F)

#### add the sample which has a zero ROH ####
setwd("/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant")
roh_tiwi <- read.table("ROH_common_run0.hom", header = T)
roh_indivi <- read.table("ROH_common_run0.indiv", header = T)
sum <- read.table("ROH_common_run0.hom.summary", header = T)
y <- setdiff(unique(roh_tiwi$FID), roh_indivi$FID)
no_roh <- data.frame(matrix(0, ncol = 19, nrow = length(sum$SNP)))
colnames(no_roh) <- c("SNP", y)
no_roh$SNP <- sum$SNP
write.table(no_roh, row.names = F, quote = F, "1_Association_files/Tiwi_data_plink_roh_regional_noROHsample.tsv")

#### merge in R itself ###
no_roh <- read.table("Tiwi_data_plink_roh_regional_noROHsample.tsv", header = T)
roh <- read.table("Tiwi_data_plink_roh_regional.tsv", header = T)
rm(no_roh)
rm(roh)
new <- merge(roh, no_roh, by.x = "SNP", by.y = "SNP", all = TRUE)
colnames(new) <- gsub("X", "", colnames(new))
colnames(new) <- gsub("\\.", "-", colnames(new))
rownames(new) <- new$SNP
sum <- read.table("../ROH_common_run0.hom.summary", header = T)
new <- new[sum$SNP,]   ### regarrange based on the sum dataset ###
rm(sum)
#### rearragne the columns ###
roh_indivi <- read.table("../ROH_common_run0.hom.indiv", header = T)
new <- new[,c("SNP", roh_indivi$IID)]
write.table(new, "Tiwi_data_plink_roh_reginal.tsv", quote = F, row.names = F)  ##### save the output file ROH format #####
#### make garlic file and convert that into plink file ##
csv_file <- "/mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/1_MAF0.05/3_IBD/MMAP_files/Tiwi_data_mmap_first_sex_columns.csv"
col <- read.csv(csv_file)
print(head(col))
rownames(col) <- col$SNPNAME
print(table(col$SNPNAME == rownames(new))) ### check the snps are in the same order
##### merge the two dataset ####
combined_data <- cbind(col, new)
write.csv(combined_data, row.names = F, quote = F, "/mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Runs_of_Homozygosity/1_ROH_analysis_plink/Common_variant/1_Association_files/Tiwi_data_regional_roh_plink.csv")
