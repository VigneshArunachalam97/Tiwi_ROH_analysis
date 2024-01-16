#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
#### import packages ###
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(qqman))
suppressPackageStartupMessages(library(cowplot))
print("Packages imported")
print("import the different function script")
source("/home/n11142006/TIwi_data/Script/Functions/Functions.R")
ggplotColors <- function(n = 6, h = c(0, 360) + 15){
                 if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
                hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
#library(RColorBrewer)
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
#colr<-sample(col_vector, 22)
#colr<-c("#CBD5E8","#FDDAEC","#FB8072","#66A61E","#6A3D9A","#FFED6F","#E31A1C","#A6761D","#FFF2AE","#984EA3","#B15928","#FED9A6","#A6CEE3","#66C2A5","#DECBE4","#F0027F","#CCCCCC","#666666","#FDBF6F","#FFD92F","#F4CAE4","#B3B3B3")
colr<-unname(unlist(rep(as.data.frame(matrix(c("#CCCCCC", "#403C3C"), nrow = 2)), each = 11)))
#### args[1]  ---- File with path
#### args[2] --- outcome prefix file name with path
#### args[3] -- title of the plot ##3
#### args[4] --- plink or GCTA
if(args[4] == "plink") {
     results<-as.data.frame(fread(args[1]))
     print(head(results))
     results <- subset(results, TEST == "ADD")
     print(dim(results))
     results <- results[,c("#CHROM", "POS", "ID", "REF", "ALT", "TEST", "BETA", "SE", "T_STAT", "P")]
     print(colnames(results))
     colnames(results) <- c("Chr", "bp", "SNP", "Ref", "Alt", "TEST", "b", "se", "test_stat", "p")
} else {
     results <- as.data.frame(fread(args[1]))
}
print("Data imported")
######### calculated genomic inflation factor (lambda) ####
GI<-round(lamba(results), 3)
print(paste0("The genomic inflation factor (lambda) is ", GI))

############## plotting started ####
head(results)
results$Chr<-as.factor(results$Chr)
names(colr)<-levels(results$Chr)

max_bp<-list()
for(i in levels(results$Chr)){
 x<-subset(results, Chr == i)
 max_bp[[i]]<-max(as.numeric(x$bp))
 }
max_bp<-as.vector(unlist(max_bp))
#bp_add<-lag(cumsum(as.numeric(max_bp)), default = 0)
data_cum<-data.frame(Chr = levels(results$Chr), bp_add = lag(cumsum(as.numeric(max_bp)), default = 0))
print("Created data_cum")
# data_cum <- results %>% group_by(Chr)  %>% summarize(max_bp = max(as.numeric(bp))) %>% mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0))
# data_cum
# data_cum<-data_cum[,c(Chr, bp_add)]
print("subset the data_cum")
data_cum$Chr<-as.factor(data_cum$Chr)
gwas_data <- results %>% inner_join(data_cum, by = "Chr") %>% mutate(bp_cum = bp + bp_add)
gwas_data <-gwas_data %>% mutate(SNP_2 = ifelse(grepl("\\.", SNP) == TRUE, paste0("Chr",Chr,":",bp), SNP))
d<-data.frame(SNP = head(gwas_data[order(gwas_data$p),]$SNP_2, 20))
axis_set<-list()
for(i in levels(gwas_data$Chr)){
 x<-subset(gwas_data, Chr == i)
 axis_set[[i]]<-mean(x$bp_cum)
}
axis_set<-data.frame(Chr = levels(results$Chr), center = unlist(axis_set))
#axis_set <- gwas_data %>% group_by(Chr) %>% summarize(center = mean(bp_cum))
ylim <- abs(floor(log10(min(gwas_data$p)))) + 2
print("Set Genome wide significacne as 1e-5 and thershold as 1e-3")
print(dim(gwas_data))
# png(paste0(args[2], "_mannhattan_plot_", format(Sys.time(), "%b_%d_%y"), ".png"), w =1800 , h = 1000, res = 350, type = "cairo", bg ="transparent")
x<-ggplot(gwas_data, aes(x = bp_cum, y = -log10(p),
                                  color = as.factor(Chr))) +
  #ggtitle("Manhattan Plot") +
  geom_hline(yintercept = -log10(1e-8), color = "DarkRed", linetype = "dashed") +
  geom_hline(yintercept = -log10(1e-5), color = "Darkblue", linetype = "dashed") +
  geom_point(size = 0.2) +
  scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = colr) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome",
       y = quote(-log10(p))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 4.5, vjust = 1),
                axis.text.y = element_text(size = 5),
                axis.title.x = element_text(size = 7, hjust = 0.5),
                axis.title.y = element_text(size = 7, hjust = 0.5),
                plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 14, color = "darkblue", hjust = 0.6),
    legend.position = "none")
# x
# dev.off()
png(paste0(args[2], "_mannhattan_plot_", format(Sys.time(), "%b_%d_%y"), ".png"), w =1800 , h = 1000, res = 350, type = "cairo", bg ="transparent")
x
dev.off()
print("GWAS plot created and saved")
#y<-qqplot(results$p, "Q-Q Plot")
print("QQ plot created")
#combine<-plot_grid(x, y, ncol =2, labels = c("a)", "b)"), rel_widths = c(5.5,3), align = "h")  ### combine the plot using plot_grid function
#### create a separate plot for title of the variable #######
#z<-ggplot() + geom_blank() + ggtitle("Albumin to Creatinine Ratio") + theme(plot.title = element_text(hjust = 0.5, size = 16, color = "Darkblue"), panel.background=element_blank())
#z<-ggplot() + geom_blank() + ggtitle(args[3]) + theme(plot.title = element_text(hjust = 0.5, size = 16, color = "Darkblue", face = "bold"), panel.background=element_blank())
#if(args[4] == "GCTA") {
     #tiff(paste0(args[2], "_manhattan_qq_combined.tiff"), res=310, type = "cairo", width=2100, height = 800, units="px")
     #plot_grid(z, combine, nrow = 2, ncol = 1, rel_heights = c(1.2, 8), align = "v")
     #dev.off()
     #print("tiff file saved")
     #svg(paste0(args[2], "_manhattan_qq_combined.svg"), width =7, height = 5)
     #plot_grid(z, combine, nrow = 2, ncol = 1, rel_heights = c(1.2, 8), align = "v")
     #dev.off()
#} else {
     #print("plot saved")
#}
