# ################################################################################################################################################################################################################
# GO.step81.create_heatmaps.R
#
# CC-BY Lee Edsall
#	email: le49@duke.edu
#	Twitter: @LeeEdsall
#
# This script was used to analyze data for Edsall et al. 2019 which compared DNase-seq data from 5 primates (human, chimpanzee, gorilla, orangutan, macaque)
#
# Creates heatmaps of the normalized scores for all of the differential categories
#	Generates all of the heatmaps and saves them in a single pdf file (this makes it easy to view and compare the heatmaps)
#	Regenerates the heatmaps and stores them individually in separate svg files (so they will be higher quaility images that can be easily be imported into drawing programs)
#
# The glm results file contains non-normalized scores so they will need to be normalized before generating the heatmaps
# The glm results file also contains a flag used to indicate to which category a DHS site belongs
#
# The relevant columns in the glm results file are
#	44-58: scores of the 15 replicates (3 replicates for each of 5 species)
#		44-46: human replicates (column names are h1, h2, h3)
#		47-49: chimpanzee replicates (column names are c1, c2, c3)
#		50-52: gorilla replicates (column names are g1, g2, g3)
#		53-55: orangutan replicates (column names are o1, o2, o3)
#		56-58: macaque replicates (column names are m1, m2, m3)
#	60: flag_change_type
#		This is a 4 character flag that indicates to which category the site belongs
#		Each position is one of 3 letters and indicates the type of change relative to macaque for each non-macaque species.
#		The possible values are s, G, L
#			s = same as macaque
#			G = Greater than macaque
#			L = Less than macaque
#		The order of the species in the flag is human, chimpanzee, gorilla, orangutan
#		Example 1: a human-specific increase has a flag of Gsss
#		Example 2: a chimpanzee/gorilla decrease has a flag of sLLs
# ################################################################################################################################################################################################################

# =====================================================================================================================================================================================
# Load libraries
# =====================================================================================================================================================================================
library(gplots)
library(RColorBrewer)

# =====================================================================================================================================================================================
# Set variables
# =====================================================================================================================================================================================
my_palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 99)
my_breaks <- 100
my_colsep <- c(3,6,9,12)
my_sepwidth <- c(0.1,0.1)
my_cexCol <- .75

# =====================================================================================================================================================================================
# Open connection to log file; write date
# =====================================================================================================================================================================================
sink(file="LOG.step82.create_heatmaps.R", append=FALSE, type="output", split=TRUE)
cat("---------------------------------------------------------------------------------------\n")
cat("BEGIN:", date(), "\n")
cat("---------------------------------------------------------------------------------------\n\n")

# =====================================================================================================================================================================================
# Load the glm results text file
# =====================================================================================================================================================================================
cat("Reading score file (all_DHS_sites.glm_analysis.all_sites.all_information.txt) ... ")

glm_results.original_scores <- as.data.frame(read.table("all_DHS_sites.glm_analysis.all_sites.all_information.txt", sep="\t", header=TRUE))

rownames(glm_results.original_scores) <- paste(glm_results.original_scores$chrom, glm_results.original_scores$start, glm_results.original_scores$end, sep=":")

cat("done\n\n")

# =====================================================================================================================================================================================
# Calculate library sizes to use as a normalization factor
# =====================================================================================================================================================================================

cat("Calculating library sizes ... ")

h1_total <- sum(glm_results.original_scores$h1)
h2_total <- sum(glm_results.original_scores$h2)
h3_total <- sum(glm_results.original_scores$h3)

c1_total <- sum(glm_results.original_scores$c1)
c2_total <- sum(glm_results.original_scores$c2)
c3_total <- sum(glm_results.original_scores$c3)

g1_total <- sum(glm_results.original_scores$g1)
g2_total <- sum(glm_results.original_scores$g2)
g3_total <- sum(glm_results.original_scores$g3)

o1_total <- sum(glm_results.original_scores$o1)
o2_total <- sum(glm_results.original_scores$o2)
o3_total <- sum(glm_results.original_scores$o3)

m1_total <- sum(glm_results.original_scores$m1)
m2_total <- sum(glm_results.original_scores$m2)
m3_total <- sum(glm_results.original_scores$m3)

cat("done.\n\n")

# =====================================================================================================================================================================================
# Normalize the scores by library size
# =====================================================================================================================================================================================

cat("Normalizing by library size ... ")

glm_results.normalized_scores <- glm_results.original_scores

glm_results.normalized_scores$h1 <- glm_results.original_scores$h1 / h1_total
glm_results.normalized_scores$h2 <- glm_results.original_scores$h2 / h2_total
glm_results.normalized_scores$h3 <- glm_results.original_scores$h3 / h3_total

glm_results.normalized_scores$c1 <- glm_results.original_scores$c1 / c1_total
glm_results.normalized_scores$c2 <- glm_results.original_scores$c2 / c2_total
glm_results.normalized_scores$c3 <- glm_results.original_scores$c3 / c3_total

glm_results.normalized_scores$g1 <- glm_results.original_scores$g1 / g1_total
glm_results.normalized_scores$g2 <- glm_results.original_scores$g2 / g2_total
glm_results.normalized_scores$g3 <- glm_results.original_scores$g3 / g3_total

glm_results.normalized_scores$o1 <- glm_results.original_scores$o1 / o1_total
glm_results.normalized_scores$o2 <- glm_results.original_scores$o2 / o2_total
glm_results.normalized_scores$o3 <- glm_results.original_scores$o3 / o3_total

glm_results.normalized_scores$m1 <- glm_results.original_scores$m1 / m1_total
glm_results.normalized_scores$m2 <- glm_results.original_scores$m2 / m2_total
glm_results.normalized_scores$m3 <- glm_results.original_scores$m3 / m3_total

cat("done.\n\n")

# =====================================================================================================================================================================================
# Generate the heatmaps to be saved in a single pdf file
#	Heatmap options
#		Rowv=TRUE		--> reorder the row dendogram
#		Colv=FALSE		--> don't reorder the column dendogram
#		symkey=F		--> don't make the color key symmetric around 0
#		symbreaks=F		--> don't make the breaks symmetric around 0
#		scale="none"		--> don't scale the values
#	Include a title and replicate labels
# =====================================================================================================================================================================================

cat("Generating the heatmaps that will be saved in a pdf file\n")

pdf("heatmap.all.pdf", width=3.5, height=5)
par(cex.main=.5)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 1 of 30: human increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="Gsss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	 main="Human Increases\n(2,539 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 2 of 30: human decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="Lsss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F, 
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human Decreases\n(209 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 3 of 30: chimpanzee increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sGss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Chimpanzee Increases\n(1,539 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 4 of 30: chimpanzee decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sLss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Chimpanzee Decreases\n(196 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 5 of 30: gorilla increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for gorilla increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="ssGs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Gorilla Increases\n(2,737 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 6 of 30: gorilla decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for gorilla decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="ssLs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Gorilla Decreases\n(254 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 7 of 30: orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for orangutan increases... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sssG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Orangutan Increases\n(3,759 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 8 of 30: orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for orangutan decreases... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sssL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Orangutan Decreases\n(646 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 9 of 30: macaque decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for macaque decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GGGG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Macaque Decreases\n(1,317 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 10 of 30: macaque increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for macaque increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LLLL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Macaque Increases\n(4,992 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 11 of 30: human/chimpanzee increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GGss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Chimpanzee Increases\n(1,735 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 12 of 30: human/chimpanzee decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LLss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Chimpanzee Decreases\n(814 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 13 of 30: human/gorilla increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/gorilla increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GsGs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Gorilla Increases\n(997 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 14 of 30: human/gorilla decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/gorilla decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LsLs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Gorilla Decreases\n(221 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 15 of 30: human/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/orangutan increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GssG")[,44:58])

heatmap.2(temp_matrix,dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Orangutan Increases\n(997 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 16 of 30: human/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/orangutan decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LssL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Orangutan Decreases\n(424 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 17 of 30: chimpanzee/gorilla increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/gorilla increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sGGs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Chimpanzee/Gorilla Increases\n(1,568 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 18 of 30: chimpanzee/gorilla decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/gorilla decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sLLs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Chimpanzee/Gorilla Decreases\n(387 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 19 of 30: chimpanzee/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/orangutan increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sGsG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Chimpanzee/Orangutan Increases\n(657 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 20 of 30: chimpanzee/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/orangutan decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sLsL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Chimpanzee/Orangutan Decreases\n(252 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 21 of 30: gorilla/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for gorilla/orangutan increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="ssGG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Gorilla/Orangutan Increases\n(2,120 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 22 of 30: gorilla/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for gorilla/orangutan decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="ssLL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Gorilla/Orangutan Decreases\n(517 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 23 of 30: human/chimpanzee/gorilla increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee/gorilla increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GGGs")[,44:58])

heatmap.2(temp_matrix,dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Chimpanzee/Gorilla Increases\n(1,736 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 24 of 30: human/chimpanzee/gorilla decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee/gorilla decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LLLs")[,44:58])

heatmap.2(temp_matrix,dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Chimpanzee/Gorilla Decreases\n(1,453 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 25 of 30: human/chimpanzee/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee/orangutan increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GGsG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Chimpanzee/Orangutan Increases\n(794 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 26 of 30: human/chimpanzee/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee/orangutan decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LLsL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Chimpanzee/Orangutan Decreases\n(906 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 27 of 30: human/gorilla/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/gorilla/orangutan increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GsGG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Gorilla/Orangutan Increases\n(691 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 28 of 30: human/gorilla/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/gorilla/orangutan decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LsLL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Human/Gorilla/Orangutan Decreases\n(554 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 29 of 30: chimpanzee/gorilla/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/gorilla/orangutan increases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sGGG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Chimpanzee/Gorilla/Orangutan Increases\n(939 sites)")

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pdf plot 30 of 30: chimpanzee/gorilla/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/gorilla/orangutan decreases ... ")

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sLLL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Chimpanzee/Gorilla/Orangutan Decreases\n(686 sites)")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up plots
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev.off()

cat("done.\n\n\n")

# =====================================================================================================================================================================================
# Generate the heatmaps to be saved in the svg files (one file per category)
#	Heatmap options
#		Rowv=TRUE		--> reorder the row dendogram
#		Colv=FALSE		--> don't reorder the column dendogram
#		symkey=F		--> don't make the color key symmetric around 0
#		symbreaks=F		--> don't make the breaks symmetric around 0
#		scale="none"		--> don't scale the values
#	Do not include a title and replicate labels (this way the svg file will contain only the heatmap)
# =====================================================================================================================================================================================

cat("Generating the heatmaps that will be saved in svg files\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 1 of 30: human increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human increases ... ")

svg(file="heatmap.Gsss.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="Gsss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	 main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 2 of 30: human decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human decreases ... ")

svg(file="heatmap.Lsss.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="Lsss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F, 
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 3 of 30: chimpanzee increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee increases ... ")

svg(file="heatmap.sGss.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sGss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 4 of 30: chimpanzee decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee decreases ... ")

svg(file="heatmap.sLss.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sLss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 5 of 30: gorilla increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for gorilla increases ... ")

svg(file="heatmap.ssGs.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="ssGs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 6 of 30: gorilla decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for gorilla decreases ... ")

svg(file="heatmap.ssLs.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="ssLs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 7 of 30: orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for orangutan increases... ")

svg(file="heatmap.sssG.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sssG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 8 of 30: orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for orangutan decreases... ")

svg(file="heatmap.sssL.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sssL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 9 of 30: macaque decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for macaque decreases ... ")

svg(file="heatmap.GGGG.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GGGG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Macaque Decreases\n(1,317 sites)")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 10 of 30: macaque increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for macaque increases ... ")

svg(file="heatmap.LLLL.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LLLL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="Macaque Increases\n(4,992 sites)")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 11 of 30: human/chimpanzee increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee increases ... ")

svg(file="heatmap.GGss.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GGss")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 12 of 30: human/chimpanzee decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee decreases ... ")

svg(file="heatmap.LLss.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LLss")[,44:58])
heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 13 of 30: human/gorilla increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/gorilla increases ... ")

svg(file="heatmap.GsGs.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GsGs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 14 of 30: human/gorilla decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/gorilla decreases ... ")

svg(file="heatmap.LsLs.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LsLs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 15 of 30: human/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/orangutan increases ... ")

svg(file="heatmap.GssG.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GssG")[,44:58])

heatmap.2(temp_matrix,dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 16 of 30: human/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/orangutan decreases ... ")

svg(file="heatmap.LssL.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LssL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 17 of 30: chimpanzee/gorilla increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/gorilla increases ... ")

svg(file="heatmap.sGGs.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sGGs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 18 of 30: chimpanzee/gorilla decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/gorilla decreases ... ")

svg(file="heatmap.sLLs.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sLLs")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 19 of 30: chimpanzee/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/orangutan increases ... ")

svg(file="heatmap.sGsG.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sGsG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 20 of 30: chimpanzee/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/orangutan decreases ... ")

svg(file="heatmap.sLsL.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sLsL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 21 of 30: gorilla/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for gorilla/orangutan increases ... ")

svg(file="heatmap.ssGG.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="ssGG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 22 of 30: gorilla/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for gorilla/orangutan decreases ... ")

svg(file="heatmap.ssLL.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="ssLL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 23 of 30: human/chimpanzee/gorilla increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee/gorilla increases ... ")

svg(file="heatmap.GGGs.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GGGs")[,44:58])

heatmap.2(temp_matrix,dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 24 of 30: human/chimpanzee/gorilla decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee/gorilla decreases ... ")

svg(file="heatmap.LLLs.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LLLs")[,44:58])

heatmap.2(temp_matrix,dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 25 of 30: human/chimpanzee/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee/orangutan increases ... ")

svg(file="heatmap.GGsG.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GGsG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 26 of 30: human/chimpanzee/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/chimpanzee/orangutan decreases ... ")

svg(file="heatmap.LLsL.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LLsL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 27 of 30: human/gorilla/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/gorilla/orangutan increases ... ")

svg(file="heatmap.GsGG.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="GsGG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 28 of 30: human/gorilla/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for human/gorilla/orangutan decreases ... ")

svg(file="heatmap.LLsL.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="LsLL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 29 of 30: chimpanzee/gorilla/orangutan increases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/gorilla/orangutan increases ... ")

svg(file="heatmap.sGGG.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sGGG")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# svg plot 30 of 30: chimpanzee/gorilla/orangutan decreases
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating heatmap for chimpanzee/gorilla/orangutan decreases ... ")

svg(file="heatmap.sLLL.svg", bg="transparent", width=5, height=5)
par(mar=c(0,0,0,0))

temp_matrix <- as.matrix(subset(glm_results.normalized_scores, flag_change_type=="sLLL")[,44:58])

heatmap.2(temp_matrix, dendrogram="none", Rowv=TRUE, Colv=FALSE, labRow="", labCol="", breaks=my_breaks, trace="none", symkey=F, symbreaks=F,
	scale="none", col=my_palette, cexCol=my_cexCol, key=FALSE, lwid=c(0.25,3.25), colsep=my_colsep, sepwidth=my_sepwidth,
	main="")

dev.off()

cat("done.\n\n\n")

# =====================================================================================================================================================================================
# Finish up script
# =====================================================================================================================================================================================

cat("Plots are located in heatmap.all.pdf and heatmap.XXXX.svg where XXXX is the change type flag\n\n")

cat("---------------------------------------------------------------------------------------\n")
cat("END:", date(), "\n")
cat("---------------------------------------------------------------------------------------\n")

sink()
