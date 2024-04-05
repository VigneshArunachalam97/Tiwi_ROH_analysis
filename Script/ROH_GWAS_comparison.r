#### import library ###
library(data.table)
library(dplyr)
library(ggvenn)
library(ggplot2)
library(cowplot)

### set the working directory ###
setwd("/mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/2_comparison_GWAS_ROH")

## import those two GWAS data ####
roh_gwas <- fread("Tiwi_data_roh_ACR.mlma")
gwas <- fread("Tiwi_data_results_ACR_annotated.tsv")


#### subset the nominal significance form both the files ####
roh_id<- subset(roh_gwas, p <= 1e-5)$SNP
gwas_id <- subset(gwas, p <= 1e-5)$SNP

id <- unique(c(roh_id, gwas_id))

###### subset the data #######
roh_gwas_subset <- as.data.frame(subset(roh_gwas, SNP %in% id))
colnames(roh_gwas_subset)[colnames(roh_gwas_subset) == "b"] <- "b_roh"
colnames(roh_gwas_subset)[colnames(roh_gwas_subset) == "se"] <- "se_roh"
colnames(roh_gwas_subset)[colnames(roh_gwas_subset) == "p"] <- "p_roh"
head(roh_gwas_subset)
gwas_subset <- as.data.frame(subset(gwas, SNP %in% id))

#### combine them both using left join ####
roh_gwas_subset$SNP <- as.character(roh_gwas_subset$SNP)
gwas_subset$SNP <- as.character(gwas_subset$SNP)
gwas_subset$b_roh <- left_join(gwas_subset, roh_gwas_subset, by = "SNP")$b_roh
gwas_subset$se_roh <- left_join(gwas_subset, roh_gwas_subset, by = "SNP")$se_roh
gwas_subset$p_roh <- left_join(gwas_subset, roh_gwas_subset, by = "SNP")$p_roh

write.table(gwas_subset, "Tiwi_data_GWAS_ROH_association_nominal_significance.tsv", row.names = F, quote = F)
####### genome wide ##########
roh_id_geno <- subset(roh_gwas, p <= 1.21e-8)$SNP
print(length(roh_id_geno))
gwas_id_geno <- subset(gwas, p <= 5e-8)$SNP
print(length(gwas_id_geno))
######## check for any intersection #####
inter <- intersect(roh_id_geno, gwas_id_geno)
print(length(inter))

##### venn diagram ###
x <- list(ROH_Associaton = roh_id_geno, 
           Traditional_GWAS = gwas_id_geno)
		   
pl <- ggvenn(x, 
			show_percentage = TRUE, 
			fill_color = c("#0073C2FF", "#EFC000FF"),
			digits = 1, 
			stroke_size = 1, 
			set_name_size = 4, text_size = 4)		   
saveRDS(pl, "venn_diagram_gene.RDS")
pl
png("venn.png", type = "cairo", w = 2000, h = 2000, units = "px", res = 300)
pl
dev.off()


geno_id <- c(roh_id_geno, gwas_id_geno)

###### subset the data #######
roh_gwas_subset <- as.data.frame(subset(roh_gwas, SNP %in% geno_id))
colnames(roh_gwas_subset)[colnames(roh_gwas_subset) == "b"] <- "b_roh"
colnames(roh_gwas_subset)[colnames(roh_gwas_subset) == "se"] <- "se_roh"
colnames(roh_gwas_subset)[colnames(roh_gwas_subset) == "p"] <- "p_roh"
head(roh_gwas_subset)
gwas_subset <- as.data.frame(subset(gwas, SNP %in% geno_id))

gwas_subset$b_roh <- left_join(gwas_subset, roh_gwas_subset, by = "SNP")$b_roh
gwas_subset$se_roh <- left_join(gwas_subset, roh_gwas_subset, by = "SNP")$se_roh
gwas_subset$p_roh <- left_join(gwas_subset, roh_gwas_subset, by = "SNP")$p_roh

write.table(gwas_subset, "Tiwi_data_GWAS_ROH_association_genome_significance.tsv", row.names = F, quote = F)

print(summary(gwas_subset$b_roh))
print(summary(gwas_subset$b))

############ scatter plot ##########
scat <- ggplot(gwas_subset, aes(x = b, y = b_roh)) + 
        geom_point(alpha = 0.6, size = 0.9) +
		xlab("Effect size - GWAS") + 
		ylab("Effect size - ROH ") +
		geom_hline(yintercept = 0, linetype = "dotted") +
		geom_vline(xintercept = 0, linetype = "dotted") + 
		theme(axis.title.x = element_text(size = 10, hjust = 0.5, color = "black", face = "bold"), 
		      axis.title.y = element_text(size = 10, hjust = 0.5, color = "black", face = "bold"),
			  axis.text = element_text(size = 8, color = "black"))

saveRDS(scat, "scatter_diagram_gene.RDS")
		
png("scatter.png", type = "cairo", w = 1080, h = 1080, units = "px", res = 300)
scat
dev.off()
### combine the plot and save as RDS object ###
whole <- plot_grid(pl, scat, nrow = 1, ncol = 2, rel_widths = c(1, 0.7), labels = c("a", "b"))
saveRDS(whole, "Whole_plot.RDS")

png("venn_scatter_plot.png", type = "cairo-png", w = 3080, h = 1080, units = "px", res = 300)
whole
dev.off()

tiff("venn_scatter_plot.tiff", type = "cairo", w = 3080, h = 1080, units = "px", res = 300)
whole
dev.off()


pdf("venn_scatter_plot.pdf", w = 8, h = 4)
whole
dev.off()


##### import the RDS object ####
Whole <- readRDS("Whole_plot.RDS")
png("venn_scatter_plot.png", type = "cairo", w = 3080, h = 1080, units = "px", res = 300)
Whole
dev.off()

pdf("venn_scatter_plot.pdf", w = 8, h = 4)
Whole
dev.off()



