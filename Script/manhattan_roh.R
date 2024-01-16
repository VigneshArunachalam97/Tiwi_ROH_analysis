#!/bin/bash
args <- commandArgs(trailingOnly = TRUE)
print(args)
##args1 - folder to set working directory
##args2 - ROH summary file for individual SNPs
##args3 - ROH indiv file (plink out) to calculate the percentage of SNPs
##args4 - output prefix with folder 
#### create a manhattan plot for SNP accurance in the ROH frequency ######
#### function for the manhattan plot data creation ###
data_gwas<-function(gwas.data){
		gwas.data$Chr<-as.factor(gwas.data$Chr)
		max_bp<-list()
		for(i in levels(gwas.data$Chr)){
				x<-subset(gwas.data, Chr == i)
				max_bp[[i]]<-max(as.numeric(x$bp))
		}
		max_bp<-as.vector(unlist(max_bp))
		data_cum<-data.frame(Chr = levels(gwas.data$Chr), bp_add = lag(cumsum(as.numeric(max_bp)), default = 0))
		data_cum$Chr<-as.factor(data_cum$Chr)
		gwas_data <- gwas.data %>% inner_join(data_cum, by = "Chr") %>% mutate(bp_cum = bp + bp_add)
		gwas_data <-gwas_data %>% mutate(SNP_2 = ifelse(grepl("\\.", SNP) == TRUE, paste0("Chr",Chr,":",bp), SNP))
		d<-data.frame(SNP = head(gwas_data[order(gwas_data$p),]$SNP_2, 20))
		axis_set<-list()
		for(i in levels(gwas_data$Chr)){
				x<-subset(gwas_data, Chr == i)
				axis_set[[i]]<-mean(x$bp_cum) 
		}
		axis_set<-data.frame(Chr = levels(gwas.data$Chr), center = unlist(axis_set))		
ab<-list(gwas_data, axis_set)
names(ab)<-c("gwas_data", "axis_set")
return(ab)
}
### import required package ###
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(data.table))
####### import the files #####
setwd(args[1]) ### set working directory
data <- as.data.frame(fread(args[2]))
indiv_data <- as.data.frame(fread(args[3]))
print(dim(indiv_data))
indiv_data <- subset(indiv_data, KB > 0) ####remove the samples with zero ROH 
paste0("The samples with non zero ROH dimensions are ", dim(indiv_data)[1])
data$prop <- (data$UNAFF/dim(indiv_data)[1])  ### calculate the proportion of the SNP in ROH regional
colnames(data) <- c("Chr", "SNP", "bp", "AFF", "UNAFF", "p")
data <- subset(data, p != 0)
print("data imported and subsetted")
### set colors for the manhattan plot
colr<-unname(unlist(rep(as.data.frame(matrix(c("#CCCCCC", "#403C3C"), nrow = 2)), each = 11)))
### create a data for manhattan plot ###
manhattan <- data_gwas(data)
manh_plot <- ggplot(data = manhattan[[1]], aes(x = bp_cum, y = p, 
                                  color = as.factor(Chr))) +
  geom_hline(yintercept = 0.15, color = "DarkRed", linetype = "dashed") + 
  geom_hline(yintercept = 0.20, color = "Darkblue", linetype = "dashed") +
  geom_point(size = 0.4) +
  scale_x_continuous(label = manhattan[[2]]$Chr, breaks = manhattan[[2]]$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(manhattan[[1]]$p + 0.03))) +
  scale_color_manual(values = colr) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chr", 
       y = "proportion of SNP in ROH region") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 5, vjust = 1, angle = 90), 
		axis.text.y = element_text(size = 7),
		axis.title.x = element_text(size = 8, hjust = 0.5),
		axis.title.y = element_text(size = 8, hjust = 0.5),
		plot.title = element_text(face = "bold", size = 12, color = "Black", hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 14, color = "darkblue", hjust = 0.6), 
    legend.position = "none") 

saveRDS(manh_plot, paste0(args[4], "_roh_manhattan_obj.RDS"))
print("saved the plot RDS object")
print("plot created and saving")
png(paste0(args[4], "_roh_manhattan_plot.png"), type = "cairo", h = 1000, w= 1800, res = 350, bg = "transparent")
manh_plot
dev.off()

pdf(paste0(args[4], "_roh_manhattan_plot.pdf"), h=5,w=10, bg="white", colormodel="cmyk")
manh_plot
dev.off()
print("plot created and saved in the given output folder")
