# ################################################################################################################################################################################################################
# GO.step83.replicate_variability.create_violin_plots.R
#
# CC-BY Lee Edsall
#	email: le49@duke.edu
#	Twitter: @LeeEdsall
#
# This script was used to analyze data for Edsall et al. 2019 which compared DNase-seq data from 5 primates (human, chimpanzee, gorilla, orangutan, macaque)

# Create violin plots of the coefficient of variability to assess replicate variability
# 	The input scores are not normalized; so normalize them first
# ###################################################################################################################################################################################################################


# ===================================================================================================================================================================================================================
# load libraries             
# ===================================================================================================================================================================================================================
library("ggplot2")
library("raster")
library("vioplot")

# ===================================================================================================================================================================================================================
# open the connection to the log file             
# ===================================================================================================================================================================================================================

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Open the log file; write date
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sink(file="LOG.step83.replicate_variability.create_violin_plots.R", append=FALSE, type="output", split=TRUE)
cat("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
cat("BEGIN:", date(), "\n")
cat("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n")

# ===================================================================================================================================================================================================================
# Load the score table
# The replicates are in the columns and the locations are in the rows
# The scores have not been normalized
# ===================================================================================================================================================================================================================

cat("Reading score file (all_DHS_sites.passed_coverage_filter.with_non_normalized_scores.zero_filtered.txt.with_PLoS_overlap_information) ... ")

scores_with_locations <- as.data.frame(read.table("all_DHS_sites.passed_coverage_filter.with_non_normalized_scores.zero_filtered.txt.with_PLoS_overlap_information", sep="\t", header=FALSE))

colnames(scores_with_locations) <- c("chrom","start","end", "h1", "h2", "h3", "c1", "c2", "c3", "g1", "g2", "g3", "o1", "o2", "o3", "m1", "m2", "m3", "PLoS_overlap")

rownames(scores_with_locations) <- paste(scores_with_locations$chrom, scores_with_locations$start, scores_with_locations$end, sep=":")

cat("done\n\n")


# ===================================================================================================================================================================================================================
# Calculate the coefficient of variation
# ===================================================================================================================================================================================================================

cat("Calculating the coefficient of variation ... ")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create the data frames
#	The first column will be the cv, so the initial values don't matter
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Human
human <- as.data.frame(scores_with_locations$h1)
human$species <- "human"
colnames(human) <- c("cv", "species")

# Chimp
chimp <- as.data.frame(scores_with_locations$h1)
chimp$species <- "chimp"
colnames(chimp) <- c("cv", "species")

# Gorilla
gorilla <- as.data.frame(scores_with_locations$h1)
gorilla$species <- "gorilla"
colnames(gorilla) <- c("cv", "species")

# Orangutan
orangutan <- as.data.frame(scores_with_locations$h1)
orangutan$species <- "orangutan"
colnames(orangutan) <- c("cv", "species")

# Macaque
macaque <- as.data.frame(scores_with_locations$h1)
macaque$species <- "macaque"
colnames(macaque) <- c("cv", "species")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate the cv
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for (lcv in 1:nrow(scores_with_locations))
{

	human[lcv,1] <- cv(c(scores_with_locations[lcv,4], scores_with_locations[lcv,5], scores_with_locations[lcv,6]))
	chimp[lcv,1] <- cv(c(scores_with_locations[lcv,7], scores_with_locations[lcv,8], scores_with_locations[lcv,9]))
	gorilla[lcv,1] <- cv(c(scores_with_locations[lcv,10], scores_with_locations[lcv,11], scores_with_locations[lcv,12]))
	orangutan[lcv,1] <- cv(c(scores_with_locations[lcv,13], scores_with_locations[lcv,14], scores_with_locations[lcv,15]))
	macaque[lcv,1] <- cv(c(scores_with_locations[lcv,16], scores_with_locations[lcv,17], scores_with_locations[lcv,18]))
}


cat("done\n\n")

# ===================================================================================================================================================================================================================
# Generate the plots
# ===================================================================================================================================================================================================================

cat("Generating violin plots ... ")

pdf("coefficient_of_variability_plots.pdf", paper="USr", height=7.5, width=10)

# Set up the plot
violin_plot <- ggplot() + theme_bw() + xlab(NULL) + ylab("Coefficient of Variation (percent)") + ylim(0,200) + scale_x_discrete(limits=c("human", "chimp", "gorilla", "orangutan", "macaque")) + 
	ggtitle("Variation between species replicates across the union set of DHS sites (89,744 sites)") + theme(plot.title = element_text(hjust = 0.5))

# Human
violin_plot <- violin_plot + geom_violin(data=human, aes(species,cv), fill="lightyellow", color="gray50", trim=FALSE)
violin_plot <- violin_plot + geom_boxplot(data=human, aes(species,cv), outlier.shape=NA, width=.25)

# Chimp
violin_plot <- violin_plot + geom_violin(data=chimp, aes(species,cv), fill="lightyellow", color="gray50", trim=FALSE)
violin_plot <- violin_plot + geom_boxplot(data=chimp, aes(species,cv), outlier.shape=NA, width=.25)

# Gorilla
violin_plot <- violin_plot + geom_violin(data=gorilla, aes(species,cv), fill="lightyellow", color="gray50", trim=FALSE)
violin_plot <- violin_plot + geom_boxplot(data=gorilla, aes(species,cv), outlier.shape=NA, width=.25)

# Orangutan
violin_plot <- violin_plot + geom_violin(data=orangutan, aes(species,cv), fill="lightyellow", color="gray50", trim=FALSE)
violin_plot <- violin_plot + geom_boxplot(data=orangutan, aes(species,cv), outlier.shape=NA, width=.25)

# Macaque
violin_plot <- violin_plot + geom_violin(data=macaque, aes(species,cv), fill="lightyellow", color="gray50", trim=FALSE)
violin_plot <- violin_plot + geom_boxplot(data=macaque, aes(species,cv), outlier.shape=NA, width=.25)

print(violin_plot)

dev.off()

# ===================================================================================================================================================================================================================
# Finish up
# ===================================================================================================================================================================================================================
cat("done\n\n") 
cat("Plots are located in coefficient_of_variability_plots.all_DHS_sites.pdf\n")

cat("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
cat("END:", date(), "\n")
cat("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

sink()

