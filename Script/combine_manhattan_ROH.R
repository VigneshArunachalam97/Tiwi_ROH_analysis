args <- commandArgs(trailingOnly = TRUE)
print(args)
### args1 - folder contains the RDS object for other population
#!/bin/bash
#### import required packages ####
library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)

setwd(args[1]) #set working directory #
files <- list.files(pattern = ".RDS")
names <- gsub("_roh_manhattan_obj.RDS", "", 
              gsub("ukb23352_", "",
		   gsub("_data", "", files)))       
names <- str_to_sentence(names)
names <- gsub("Ni", "Norfolk Islander", 
              gsub("Pacific", "Pacific Islander", names))
			  
plot <- list()
for(i in 1:length(files)){
    pl <- readRDS(files[i])
	pl <- pl + ggtitle(names[i]) + theme(plot.title = element_text(size = 14, face = "bold", color = "Darkblue"), 
	                                     axis.title.y = element_text(size = 10, face = "bold", color = "black"))
	plot[[i]] <- pl
}

final_plot <- plot_grid(plotlist = plot, ncol = 2, nrow = (length(files)/5), rel_widths = 1, rel_heights = 1)

png("ROH_genomeWide_all_Population.png", w = 3600, h = 1000*(length(files)/5), res = 300, type = "cairo")
final_plot
dev.off()
print("plot saved in the given folder")
