# ################################################################################################################################################################################################################
# GO.step82.replicate_variability.create_scatterplots.R
#
# CC-BY Lee Edsall
#	email: le49@duke.edu
#	Twitter: @LeeEdsall
#
# This script was used to analyze data for Edsall et al. 2019 which compared DNase-seq data from 5 primates (human, chimpanzee, gorilla, orangutan, macaque)

# Create scatterplots of the cut counts to assess replicate variability
#	Plot each species separately; save species plots in a png file; one png file per species (saving all of them in one pdf makes the pdf file huge)
#	The input scores are not normalized; so normalize them first
# ###################################################################################################################################################################################################################

# ===================================================================================================================================================================================================================
# load libraries
# ===================================================================================================================================================================================================================
library("ggplot2")
library("grid")
library("gridExtra")
library("gtable")
library("scales")

# ===================================================================================================================================================================================================================
# set options
# ===================================================================================================================================================================================================================
options(scipen=0)

# ===================================================================================================================================================================================================================
# open connection to log file
# ===================================================================================================================================================================================================================

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Open the log file; write date
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sink(file="LOG.step82.replicate_variability.create_scatterplots.R", append=FALSE, type="output", split=TRUE)
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
# Calculate library sizes
# ===================================================================================================================================================================================================================

cat("Calculating library sizes ... ")

h1_total <- sum(scores_with_locations$h1)
h2_total <- sum(scores_with_locations$h2)
h3_total <- sum(scores_with_locations$h3)

c1_total <- sum(scores_with_locations$c1)
c2_total <- sum(scores_with_locations$c2)
c3_total <- sum(scores_with_locations$c3)

g1_total <- sum(scores_with_locations$g1)
g2_total <- sum(scores_with_locations$g2)
g3_total <- sum(scores_with_locations$g3)

o1_total <- sum(scores_with_locations$o1)
o2_total <- sum(scores_with_locations$o2)
o3_total <- sum(scores_with_locations$o3)

m1_total <- sum(scores_with_locations$m1)
m2_total <- sum(scores_with_locations$m2)
m3_total <- sum(scores_with_locations$m3)

cat("done.\n\n")

# ===================================================================================================================================================================================================================
# Normalize scores by library size
# ===================================================================================================================================================================================================================

cat("Normalizing by library size ... ")

normalized_scores <- scores_with_locations[,4:18]

normalized_scores$h1 <- scores_with_locations$h1 / h1_total
normalized_scores$h2 <- scores_with_locations$h2 / h2_total
normalized_scores$h3 <- scores_with_locations$h3 / h3_total

normalized_scores$c1 <- scores_with_locations$c1 / c1_total
normalized_scores$c2 <- scores_with_locations$c2 / c2_total
normalized_scores$c3 <- scores_with_locations$c3 / c3_total

normalized_scores$g1 <- scores_with_locations$g1 / g1_total
normalized_scores$g2 <- scores_with_locations$g2 / g2_total
normalized_scores$g3 <- scores_with_locations$g3 / g3_total

normalized_scores$o1 <- scores_with_locations$o1 / o1_total
normalized_scores$o2 <- scores_with_locations$o2 / o2_total
normalized_scores$o3 <- scores_with_locations$o3 / o3_total

normalized_scores$m1 <- scores_with_locations$m1 / m1_total
normalized_scores$m2 <- scores_with_locations$m2 / m2_total
normalized_scores$m3 <- scores_with_locations$m3 / m3_total

cat("done.\n\n")

# ###################################################################################################################################################################################################################
# Generate scatterplots
# ###################################################################################################################################################################################################################

# ===================================================================================================================================================================================================================
# Human
# ===================================================================================================================================================================================================================

cat("Generating scatterplots for human replicates ... ")

png("replicate_scatterplots.human.png", height=2250, width=2100, res=300)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 1 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

h1_h1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=h1, y=h1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

h1_h2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=h1, y=h2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) +
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

h1_h3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=h1, y=h3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) +
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 2 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

h2_h1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=h2, y=h1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

h2_h2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=h2, y=h2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

h2_h3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=h2, y=h3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 3 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

h3_h1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=h3, y=h1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

h3_h2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=h3, y=h2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

h3_h3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=h3, y=h3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create final plot
# 	The layout for the gtable is 5 columns x 8 rows (7.25 inches high x 6.5 inches wide)
# 		The 5 columns are:
#			1. plot 1 for the replicate (2 inches)
#			2. spacer (.25 inches)
#			3. plot 2 for the replicate (2 inches)
#			4. spacer (.25 inches)
#			5. plot 3 for the replicate (2 inches)
#	 	The 8 rows are:
#			1: replicate 1 title (.25 inches)
#			2: replicate 1 plots (2 inches)
#			3: spacer (.25 inches)
#			4: replicate 2 title (.25 inches)
#			5: replicate 2 plots (2 inches)
#			6: spacer (.25 inches)
#			7: replicate 3 title (.25 inches)
#			8: replicate 3 plots (2 inches)
#	The parameters are:
#		t --> top
#		l --> left
#		r --> right
#		b- --> bottom
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create the table
human_gtable <- gtable(widths=unit(c(2, .25, 2, .25, 2), "in"), heights=unit(c(.2, 2, .25, .2, 2, .25, .2, 2), "in"))

# Add a white background to all of the cells
human_gtable <- gtable_add_grob(human_gtable, rectGrob(gp = gpar(fill = "white", col="white")), t=1, l=1, b=8, r=5)

# Add the headings for each of the replicate rows
human_gtable <- gtable_add_grob(human_gtable, textGrob("Replicate 1", gp=gpar(col="blue", fontsize=12)), t=1, l=1, r=5)
human_gtable <- gtable_add_grob(human_gtable, textGrob("Replicate 2", gp=gpar(col="blue", fontsize=12)), t=4, l=1, r=5)
human_gtable <- gtable_add_grob(human_gtable, textGrob("Replicate 3", gp=gpar(col="blue", fontsize=12)), t=7, l=1, r=5)

# Add the replicate 1 plots
human_gtable <- gtable_add_grob(human_gtable, h1_h1.plot, t=2, l=1)
human_gtable <- gtable_add_grob(human_gtable, h1_h2.plot, t=2, l=3)
human_gtable <- gtable_add_grob(human_gtable, h1_h3.plot, t=2, l=5)

# Add the replicate 2 plots
human_gtable <- gtable_add_grob(human_gtable, h2_h1.plot, t=5, l=1)
human_gtable <- gtable_add_grob(human_gtable, h2_h2.plot, t=5, l=3)
human_gtable <- gtable_add_grob(human_gtable, h2_h3.plot, t=5, l=5)

# Add the replicate 3 plots
human_gtable <- gtable_add_grob(human_gtable, h3_h1.plot, t=8, l=1)
human_gtable <- gtable_add_grob(human_gtable, h3_h2.plot, t=8, l=3)
human_gtable <- gtable_add_grob(human_gtable, h3_h3.plot, t=8, l=5)

# Plot the elements
grid.arrange(human_gtable, top=textGrob("Human cut counts for the union set of DHS sites (89,744 sites)", gp=gpar(col="black", fontface="bold")), ncol=1, nrow=1)

dev.off()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("done\n\n") 


# ===================================================================================================================================================================================================================
# Chimp
# ===================================================================================================================================================================================================================

cat("Generating scatterplots for chimp replicates ... ")

png("replicate_scatterplots.chimp.png", height=2250, width=2100, res=300)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 1 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

c1_c1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=c1, y=c1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

c1_c2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=c1, y=c2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

c1_c3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=c1, y=c3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 2 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

c2_c1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=c2, y=c1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

c2_c2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=c2, y=c2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

c2_c3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=c2, y=c3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 3 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

c3_c1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=c3, y=c1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

c3_c2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=c3, y=c2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

c3_c3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=c3, y=c3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,3e-04)) + scale_y_continuous(labels = scientific, limits=c(0,3e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create final plot
# 	The layout for the gtable is 5 columns x 8 rows (7.25 inches high x 6.5 inches wide)
# 		The 5 columns are:
#			1. plot 1 for the replicate (2 inches)
#			2. spacer (.25 inches)
#			3. plot 2 for the replicate (2 inches)
#			4. spacer (.25 inches)
#			5. plot 3 for the replicate (2 inches)
#	 	The 8 rows are:
#			1: replicate 1 title (.25 inches)
#			2: replicate 1 plots (2 inches)
#			3: spacer (.25 inches)
#			4: replicate 2 title (.25 inches)
#			5: replicate 2 plots (2 inches)
#			6: spacer (.25 inches)
#			7: replicate 3 title (.25 inches)
#			8: replicate 3 plots (2 inches)
#	The parameters are:
#		t --> top
#		l --> left
#		r --> right
#		b- --> bottom
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create the table
chimp_gtable <- gtable(widths=unit(c(2, .25, 2, .25, 2), "in"), heights=unit(c(.2, 2, .25, .2, 2, .25, .2, 2), "in"))

# Add a white background to all of the cells
chimp_gtable <- gtable_add_grob(chimp_gtable, rectGrob(gp = gpar(fill = "white", col="white")), t=1, l=1, b=8, r=5)

# Add the headings for each of the replicate rows
chimp_gtable <- gtable_add_grob(chimp_gtable, textGrob("Replicate 1", gp=gpar(col="blue", fontsize=12)), t=1, l=1, r=5)
chimp_gtable <- gtable_add_grob(chimp_gtable, textGrob("Replicate 2", gp=gpar(col="blue", fontsize=12)), t=4, l=1, r=5)
chimp_gtable <- gtable_add_grob(chimp_gtable, textGrob("Replicate 3", gp=gpar(col="blue", fontsize=12)), t=7, l=1, r=5)

# Add the replicate 1 plots
chimp_gtable <- gtable_add_grob(chimp_gtable, c1_c1.plot, t=2, l=1)
chimp_gtable <- gtable_add_grob(chimp_gtable, c1_c2.plot, t=2, l=3)
chimp_gtable <- gtable_add_grob(chimp_gtable, c1_c3.plot, t=2, l=5)

# Add the replicate 2 plots
chimp_gtable <- gtable_add_grob(chimp_gtable, c2_c1.plot, t=5, l=1)
chimp_gtable <- gtable_add_grob(chimp_gtable, c2_c2.plot, t=5, l=3)
chimp_gtable <- gtable_add_grob(chimp_gtable, c2_c3.plot, t=5, l=5)

# Add the replicate 3 plots
chimp_gtable <- gtable_add_grob(chimp_gtable, c3_c1.plot, t=8, l=1)
chimp_gtable <- gtable_add_grob(chimp_gtable, c3_c2.plot, t=8, l=3)
chimp_gtable <- gtable_add_grob(chimp_gtable, c3_c3.plot, t=8, l=5)

# Plot the elements
grid.arrange(chimp_gtable, top=textGrob("Chimp cut counts for the union set of DHS sites (89,744 sites)", gp=gpar(col="black", fontface="bold")), ncol=1, nrow=1)

dev.off()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("done\n\n") 

# ===================================================================================================================================================================================================================
# Gorilla
# ===================================================================================================================================================================================================================

cat("Generating scatterplots for gorilla replicates ... ")

png("replicate_scatterplots.gorilla.png", height=2250, width=2100, res=300)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 1 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

g1_g1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=g1, y=g1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

g1_g2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=g1, y=g2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

g1_g3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=g1, y=g3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 2 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

g2_g1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=g2, y=g1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

g2_g2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=g2, y=g2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

g2_g3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=g2, y=g3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 3 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

g3_g1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=g3, y=g1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

g3_g2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=g3, y=g2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

g3_g3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=g3, y=g3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create final plot
# 	The layout for the gtable is 5 columns x 8 rows (7.25 inches high x 6.5 inches wide)
# 		The 5 columns are:
#			1. plot 1 for the replicate (2 inches)
#			2. spacer (.25 inches)
#			3. plot 2 for the replicate (2 inches)
#			4. spacer (.25 inches)
#			5. plot 3 for the replicate (2 inches)
#	 	The 8 rows are:
#			1: replicate 1 title (.25 inches)
#			2: replicate 1 plots (2 inches)
#			3: spacer (.25 inches)
#			4: replicate 2 title (.25 inches)
#			5: replicate 2 plots (2 inches)
#			6: spacer (.25 inches)
#			7: replicate 3 title (.25 inches)
#			8: replicate 3 plots (2 inches)
#	The parameters are:
#		t --> top
#		l --> left
#		r --> right
#		b- --> bottom
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create the table
gorilla_gtable <- gtable(widths=unit(c(2, .25, 2, .25, 2), "in"), heights=unit(c(.2, 2, .25, .2, 2, .25, .2, 2), "in"))

# Add a white background to all of the cells
gorilla_gtable <- gtable_add_grob(gorilla_gtable, rectGrob(gp = gpar(fill = "white", col="white")), t=1, l=1, b=8, r=5)

# Add the headings for each of the replicate rows
gorilla_gtable <- gtable_add_grob(gorilla_gtable, textGrob("Replicate 1", gp=gpar(col="blue", fontsize=12)), t=1, l=1, r=5)
gorilla_gtable <- gtable_add_grob(gorilla_gtable, textGrob("Replicate 2", gp=gpar(col="blue", fontsize=12)), t=4, l=1, r=5)
gorilla_gtable <- gtable_add_grob(gorilla_gtable, textGrob("Replicate 3", gp=gpar(col="blue", fontsize=12)), t=7, l=1, r=5)

# Add the replicate 1 plots
gorilla_gtable <- gtable_add_grob(gorilla_gtable, g1_g1.plot, t=2, l=1)
gorilla_gtable <- gtable_add_grob(gorilla_gtable, g1_g2.plot, t=2, l=3)
gorilla_gtable <- gtable_add_grob(gorilla_gtable, g1_g3.plot, t=2, l=5)

# Add the replicate 2 plots
gorilla_gtable <- gtable_add_grob(gorilla_gtable, g2_g1.plot, t=5, l=1)
gorilla_gtable <- gtable_add_grob(gorilla_gtable, g2_g2.plot, t=5, l=3)
gorilla_gtable <- gtable_add_grob(gorilla_gtable, g2_g3.plot, t=5, l=5)

# Add the replicate 3 plots
gorilla_gtable <- gtable_add_grob(gorilla_gtable, g3_g1.plot, t=8, l=1)
gorilla_gtable <- gtable_add_grob(gorilla_gtable, g3_g2.plot, t=8, l=3)
gorilla_gtable <- gtable_add_grob(gorilla_gtable, g3_g3.plot, t=8, l=5)

# Plot the elements
grid.arrange(gorilla_gtable, top=textGrob("Gorilla cut counts for the union set of DHS sites (89,744 sites)", gp=gpar(col="black", fontface="bold")), ncol=1, nrow=1)

dev.off()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("done\n\n") 

# ===================================================================================================================================================================================================================
# Orangutan
# ===================================================================================================================================================================================================================

cat("Generating scatterplots for orangutan replicates ... ")

png("replicate_scatterplots.orangutan.png", height=2250, width=2100, res=300)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 1 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

o1_o1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=o1, y=o1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

o1_o2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=o1, y=o2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

o1_o3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=o1, y=o3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 2 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

o2_o1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=o2, y=o1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

o2_o2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=o2, y=o2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

o2_o3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=o2, y=o3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 3 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

o3_o1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=o3, y=o1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

o3_o2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=o3, y=o2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

o3_o3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=o3, y=o3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,2e-04)) + scale_y_continuous(labels = scientific, limits=c(0,2e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create final plot
# 	The layout for the gtable is 5 columns x 8 rows (7.25 inches high x 6.5 inches wide)
# 		The 5 columns are:
#			1. plot 1 for the replicate (2 inches)
#			2. spacer (.25 inches)
#			3. plot 2 for the replicate (2 inches)
#			4. spacer (.25 inches)
#			5. plot 3 for the replicate (2 inches)
#	 	The 8 rows are:
#			1: replicate 1 title (.25 inches)
#			2: replicate 1 plots (2 inches)
#			3: spacer (.25 inches)
#			4: replicate 2 title (.25 inches)
#			5: replicate 2 plots (2 inches)
#			6: spacer (.25 inches)
#			7: replicate 3 title (.25 inches)
#			8: replicate 3 plots (2 inches)
#	The parameters are:
#		t --> top
#		l --> left
#		r --> right
#		b- --> bottom
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create the table
orangutan_gtable <- gtable(widths=unit(c(2, .25, 2, .25, 2), "in"), heights=unit(c(.2, 2, .25, .2, 2, .25, .2, 2), "in"))

# Add a white background to all of the cells
orangutan_gtable <- gtable_add_grob(orangutan_gtable, rectGrob(gp = gpar(fill = "white", col="white")), t=1, l=1, b=8, r=5)

# Add the headings for each of the replicate rows
orangutan_gtable <- gtable_add_grob(orangutan_gtable, textGrob("Replicate 1", gp=gpar(col="blue", fontsize=12)), t=1, l=1, r=5)
orangutan_gtable <- gtable_add_grob(orangutan_gtable, textGrob("Replicate 2", gp=gpar(col="blue", fontsize=12)), t=4, l=1, r=5)
orangutan_gtable <- gtable_add_grob(orangutan_gtable, textGrob("Replicate 3", gp=gpar(col="blue", fontsize=12)), t=7, l=1, r=5)

# Add the replicate 1 plots
orangutan_gtable <- gtable_add_grob(orangutan_gtable, o1_o1.plot, t=2, l=1)
orangutan_gtable <- gtable_add_grob(orangutan_gtable, o1_o2.plot, t=2, l=3)
orangutan_gtable <- gtable_add_grob(orangutan_gtable, o1_o3.plot, t=2, l=5)

# Add the replicate 2 plots
orangutan_gtable <- gtable_add_grob(orangutan_gtable, o2_o1.plot, t=5, l=1)
orangutan_gtable <- gtable_add_grob(orangutan_gtable, o2_o2.plot, t=5, l=3)
orangutan_gtable <- gtable_add_grob(orangutan_gtable, o2_o3.plot, t=5, l=5)

# Add the replicate 3 plots
orangutan_gtable <- gtable_add_grob(orangutan_gtable, o3_o1.plot, t=8, l=1)
orangutan_gtable <- gtable_add_grob(orangutan_gtable, o3_o2.plot, t=8, l=3)
orangutan_gtable <- gtable_add_grob(orangutan_gtable, o3_o3.plot, t=8, l=5)

# Plot the elements
grid.arrange(orangutan_gtable, top=textGrob("Orangutan cut counts for the union set of DHS sites (89,744 sites)", gp=gpar(col="black", fontface="bold")), ncol=1, nrow=1)

dev.off()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("done\n\n") 


# ===================================================================================================================================================================================================================
# Macaque
# ===================================================================================================================================================================================================================

cat("Generating scatterplots for macaque replicates ... ")

png("replicate_scatterplots.macaque.png", height=2250, width=2100, res=300)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 1 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

m1_m1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=m1, y=m1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,5e-04)) + scale_y_continuous(labels = scientific, limits=c(0,5e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

m1_m2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=m1, y=m2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,5e-04)) + scale_y_continuous(labels = scientific, limits=c(0,5e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

m1_m3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=m1, y=m3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,5e-04)) + scale_y_continuous(labels = scientific, limits=c(0,5e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 1") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 2 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

m2_m1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=m2, y=m1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,5e-04)) + scale_y_continuous(labels = scientific, limits=c(0,5e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

m2_m2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=m2, y=m2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,5e-04)) + scale_y_continuous(labels = scientific, limits=c(0,5e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

m2_m3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=m2, y=m3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,5e-04)) + scale_y_continuous(labels = scientific, limits=c(0,5e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 2") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Replicate 3 plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

m3_m1.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=m3, y=m1)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,5e-04)) + scale_y_continuous(labels = scientific, limits=c(0,5e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 1") + 
	geom_abline(intercept=0, slope=1, color="red"))

m3_m2.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=m3, y=m2)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,5e-04)) + scale_y_continuous(labels = scientific, limits=c(0,5e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 2") + 
	geom_abline(intercept=0, slope=1, color="red"))

m3_m3.plot <- ggplotGrob(ggplot(normalized_scores, aes(x=m3, y=m3)) + geom_point(size=.5) + 
	scale_x_continuous(labels = scientific, limits=c(0,5e-04)) + scale_y_continuous(labels = scientific, limits=c(0,5e-04)) + 
	theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, hjust = 1)) + theme(axis.text.y=element_text(size=10)) + xlab("Replicate 3") + ylab("Replicate 3") + 
	geom_abline(intercept=0, slope=1, color="red"))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create final plot
# 	The layout for the gtable is 5 columns x 8 rows (7.25 inches high x 6.5 inches wide)
# 		The 5 columns are:
#			1. plot 1 for the replicate (2 inches)
#			2. spacer (.25 inches)
#			3. plot 2 for the replicate (2 inches)
#			4. spacer (.25 inches)
#			5. plot 3 for the replicate (2 inches)
#	 	The 8 rows are:
#			1: replicate 1 title (.25 inches)
#			2: replicate 1 plots (2 inches)
#			3: spacer (.25 inches)
#			4: replicate 2 title (.25 inches)
#			5: replicate 2 plots (2 inches)
#			6: spacer (.25 inches)
#			7: replicate 3 title (.25 inches)
#			8: replicate 3 plots (2 inches)
#	The parameters are:
#		t --> top
#		l --> left
#		r --> right
#		b- --> bottom
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create the table
macaque_gtable <- gtable(widths=unit(c(2, .25, 2, .25, 2), "in"), heights=unit(c(.2, 2, .25, .2, 2, .25, .2, 2), "in"))

# Add a white background to all of the cells
macaque_gtable <- gtable_add_grob(macaque_gtable, rectGrob(gp = gpar(fill = "white", col="white")), t=1, l=1, b=8, r=5)

# Add the headings for each of the replicate rows
macaque_gtable <- gtable_add_grob(macaque_gtable, textGrob("Replicate 1", gp=gpar(col="blue", fontsize=12)), t=1, l=1, r=5)
macaque_gtable <- gtable_add_grob(macaque_gtable, textGrob("Replicate 2", gp=gpar(col="blue", fontsize=12)), t=4, l=1, r=5)
macaque_gtable <- gtable_add_grob(macaque_gtable, textGrob("Replicate 3", gp=gpar(col="blue", fontsize=12)), t=7, l=1, r=5)

# Add the replicate 1 plots
macaque_gtable <- gtable_add_grob(macaque_gtable, m1_m1.plot, t=2, l=1)
macaque_gtable <- gtable_add_grob(macaque_gtable, m1_m2.plot, t=2, l=3)
macaque_gtable <- gtable_add_grob(macaque_gtable, m1_m3.plot, t=2, l=5)

# Add the replicate 2 plots
macaque_gtable <- gtable_add_grob(macaque_gtable, m2_m1.plot, t=5, l=1)
macaque_gtable <- gtable_add_grob(macaque_gtable, m2_m2.plot, t=5, l=3)
macaque_gtable <- gtable_add_grob(macaque_gtable, m2_m3.plot, t=5, l=5)

# Add the replicate 3 plots
macaque_gtable <- gtable_add_grob(macaque_gtable, m3_m1.plot, t=8, l=1)
macaque_gtable <- gtable_add_grob(macaque_gtable, m3_m2.plot, t=8, l=3)
macaque_gtable <- gtable_add_grob(macaque_gtable, m3_m3.plot, t=8, l=5)

# Plot the elements
grid.arrange(macaque_gtable, top=textGrob("Macaque cut counts for the union set of DHS sites (89,744 sites)", gp=gpar(col="black", fontface="bold")), ncol=1, nrow=1)

dev.off()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("done\n\n") 


# ###################################################################################################################################################################################################################
# Finish up
# ###################################################################################################################################################################################################################
cat("Plots are located in:\n")
cat("     replicate_scatterplots.chimp.png\n")
cat("     replicate_scatterplots.gorilla.png\n")
cat("     replicate_scatterplots.human.png\n")
cat("     replicate_scatterplots.macaque.png\n")
cat("     replicate_scatterplots.orangutan.png\n\n")

cat("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
cat("END:", date(), "\n")
cat("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

sink()

