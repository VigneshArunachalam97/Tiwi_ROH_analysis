#!/bin/bash
args = commandArgs(trailingOnly=TRUE)
print(args)
##args[1] - file path
##args[2] - pattern either ".loco.mlma or glm.linear"
##args[3] - "plink" or "gcta"
print("Hi, Starting to Run")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(data.table))
print("imported packages")
x<-args[1]
print(x)
setwd(x)
gc()
print("data importing")
annotate<-as.data.frame(fread("/mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/Annovar_annotation/Tiwi_annotation.hg38_multianno_subset.txt"))
print(dim(annotate))
#annotate<-read.table("/home/n11142006/TIwi_data/Genotypic_data/Annotation_Final.txt", header = FALSE, na.strings = ".")
#colnames(annotate)<-c("Gene_position", "Gene", "Distance", "Chr", "Start", "rsID", "Ref", "Alt", "Chr_id")
if(args[3] == "plink") {
annotate$ID<-paste0(annotate$Chr,":",annotate$Start,":",annotate$Ref,":",annotate$Alt)
} else {
annotate$SNP<-paste0(annotate$Chr,":",annotate$Start,":",annotate$Ref,":",annotate$Alt)
}
annotate <- annotate[,-c(1:5)]
print(dim(annotate))
#annotate$Distance<-gsub("dist=NONE;","",gsub(";dist=NONE","",annotate$Distance))
print("imported annotation files")
head(annotate)
gc()
files<-list.files(pattern=args[2], include.dirs = FALSE, full.names = FALSE)
print(files)
gc()
print("Start importing the files")
for(i in 1:length(files)){
gwas_out<-read.csv(files[i], header = T, sep = "\t")
if(args[3] == "plink"){
         gwas_out <- subset(gwas_out, TEST == "ADD")
         gwas_out <- left_join(gwas_out, annotate, by = "ID")
         write.table(gwas_out, sep = "\t", file = gsub(args[2],"_annotated.tsv",files[i]), row.names=FALSE, quote = F)
         write.table(subset(gwas_out, P<1e-5), sep = "\t", file = gsub(args[2], "_annotated_subset.tsv", files[i]), row.names = FALSE, quote=F)
 } else {
         gwas_out<-left_join(gwas_out, annotate, by = "SNP")
         write.table(gwas_out, sep = "\t", file = gsub(args[2],"_annotated.tsv",files[i]), row.names=FALSE, quote = F)
         write.table(subset(gwas_out, p<1e-5), sep = "\t", file = gsub(args[2], "_annotated_subset.tsv", files[i]), row.names = FALSE, quote=F)
 }
}
#gwas_out$Chr_id<-paste0("chr",gwas_out$Chr,":",gwas_out$bp,":",gwas_out$A2,":",gwas_out$A1)
#gc()
#gwas_out$Gene_region<-left_join(gwas_out, annotate, by = "Chr_id")$Gene_position
#gc()
#gwas_out$Gene<-left_join(gwas_out, annotate, by = "Chr_id")$Gene
#gc()
#gwas_out$rsID_annovar<-left_join(gwas_out, annotate,by="Chr_id")$rsID
#gc()
#gwas_out$Distance<-left_join(gwas_out, annotate, by = "Chr_id")$Distance
gc()
#print(files[i])
print("Files merged and saved")
#head(gwas_out)
#write.table(gwas_out, sep = "\t", file = gsub(args[2],"_annotated.tsv",files[i]), row.names=FALSE, quote = F)
#write.table(subset(gwas_out, P<1e-5), sep = "\t", file = gsub(args[2], "_annotated_subset.tsv", files[i]), row.names = FALSE, quote=F)
gc()
print(paste0("All files saved in ", args[1]))
print("Thank you")
