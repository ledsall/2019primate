# ################################################################################################################################################################################################################
# GO.step10.run_glm.R
#
# CC-BY Lee Edsall
#	email: le49@duke.edu
#	Twitter: @LeeEdsall
#
# This script was used to analyze data for Edsall et al. 2019 which compared DNase-seq data from 5 primates (human, chimpanzee, gorilla, orangutan, macaque)
#
# Perform PCA analysis on all DHS sites using prcomp
# ################################################################################################################################################################################################################


# ==================================================================================================================
# Set variables and open connection to log file
# ==================================================================================================================

plot_symbol <- 21

species_colors <-c("Orange", "Orange", "Orange",
		   "Black", "Black", "Black",
		   "Red", "Red", "Red",
		   "Blue", "Blue", "Blue",
		   "Green", "Green", "Green")

replicate_labels <- c("h1", "h2", "h3",
		      "c1", "c2", "c3",
		      "g1", "g2", "g3",
		      "o1", "o2", "o3",
		      "m1", "m2", "m3")

# ------------------------------------------------------------------------------------------------------------------
# Open the log file; write date
# ------------------------------------------------------------------------------------------------------------------
sink(file="LOG.step91.run_pca.all_sites.R", append=FALSE, type="output", split=TRUE)
cat("---------------------------------------------------------------------------------------\n")
cat("BEGIN:", date(), "\n")
cat("---------------------------------------------------------------------------------------\n\n")


# ==================================================================================================================
# Load the table and add column names and row names
# The replicates are in the columns and the locations are in the rows
# The scores have not been normalized
# ==================================================================================================================
cat("Reading score file (all_DHS_sites.passed_coverage_filter.with_non_normalized_scores.zero_filtered.txt.with_PLoS_overlap_information) ... ")

scores_with_locations <- as.data.frame(read.table("all_DHS_sites.passed_coverage_filter.with_non_normalized_scores.zero_filtered.txt.with_PLoS_overlap_information", sep="\t", header=FALSE))

colnames(scores_with_locations) <- c("chrom","start","end", "h1", "h2", "h3", "c1", "c2", "c3", "g1", "g2", "g3", "o1", "o2", "o3", "m1", "m2", "m3", "PLoS_overlap")

rownames(scores_with_locations) <- paste(scores_with_locations$chrom, scores_with_locations$start, scores_with_locations$end, sep=":")

cat("done\n\n")

# ==================================================================================================================
# Calculate library sizes
# ==================================================================================================================

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

# ==================================================================================================================
# Normalize the scores by library size
# ==================================================================================================================

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

cat("done.\n")

# ==================================================================================================================
# Transpose the matrix so that the locations are in the columns and the samples are in the rows
# Run prcomp
# ==================================================================================================================

cat("Running prcomp ... ")

scores <- t(normalized_scores)

prcomp_results <- prcomp(scores, center=TRUE, scale=TRUE)

cat("done\n\n")

cat("---------------------------------------------------------------------------------------\n")
cat("Summary of the principal component analysis\n")
cat("---------------------------------------------------------------------------------------\n")
print(summary(prcomp_results))
cat("\n")


cat("---------------------------------------------------------------------------------------\n")
cat("Values\n")
cat("---------------------------------------------------------------------------------------\n")
print(prcomp_results$x)
cat("\n")

# ------------------------------------------------------------------------------------------------------------------
# Get the variances for each of the principal components so it can be printed on the axes
# ------------------------------------------------------------------------------------------------------------------
PC1  <- capture.output(cat("PC1 (", sprintf("%.0f",summary(prcomp_results)$importance[2,1]*100), "%)", sep=""))
PC2  <- capture.output(cat("PC2 (", sprintf("%.0f",summary(prcomp_results)$importance[2,2]*100), "%)", sep=""))
PC3  <- capture.output(cat("PC3 (", sprintf("%.0f",summary(prcomp_results)$importance[2,3]*100), "%)", sep=""))
PC4  <- capture.output(cat("PC4 (", sprintf("%.0f",summary(prcomp_results)$importance[2,4]*100), "%)", sep=""))
PC5  <- capture.output(cat("PC5 (", sprintf("%.0f",summary(prcomp_results)$importance[2,5]*100), "%)", sep=""))
PC6  <- capture.output(cat("PC6 (", sprintf("%.0f",summary(prcomp_results)$importance[2,6]*100), "%)", sep=""))
PC7  <- capture.output(cat("PC7 (", sprintf("%.0f",summary(prcomp_results)$importance[2,7]*100), "%)", sep=""))
PC8  <- capture.output(cat("PC8 (", sprintf("%.0f",summary(prcomp_results)$importance[2,8]*100), "%)", sep=""))
PC9  <- capture.output(cat("PC9 (", sprintf("%.0f",summary(prcomp_results)$importance[2,9]*100), "%)", sep=""))
PC10 <- capture.output(cat("PC10 (", sprintf("%.0f",summary(prcomp_results)$importance[2,10]*100), "%)", sep=""))
PC11 <- capture.output(cat("PC11 (", sprintf("%.0f",summary(prcomp_results)$importance[2,11]*100), "%)", sep=""))
PC12 <- capture.output(cat("PC12 (", sprintf("%.0f",summary(prcomp_results)$importance[2,12]*100), "%)", sep=""))
PC13 <- capture.output(cat("PC13 (", sprintf("%.0f",summary(prcomp_results)$importance[2,13]*100), "%)", sep=""))
PC14 <- capture.output(cat("PC14 (", sprintf("%.0f",summary(prcomp_results)$importance[2,14]*100), "%)", sep=""))
PC15 <- capture.output(cat("PC15 (", sprintf("%.0f",summary(prcomp_results)$importance[2,15]*100), "%)", sep=""))

# ==================================================================================================================
# Generate plots
# The species have different colors; the plot symbols are all the same
# Do not use equal scales - it masks the differences in the less important principal components
# ==================================================================================================================

cat("Generating PCA plots ... ")

pdf("pca.all_DHS_sites.pdf")

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC2
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,2], xlab=PC1, ylab=PC2, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC2", 
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,2], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC3
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,3], xlab=PC1, ylab=PC3, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC3",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75) 

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,3], replicate_labels, cex=.5, pos=1)


# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC4
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,4], xlab=PC1, ylab=PC4, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC4",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,4], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC5
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,5], xlab=PC1, ylab=PC5, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC5",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,5], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC6
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,6], xlab=PC1, ylab=PC6, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC6",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,6], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC7
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,7], xlab=PC1, ylab=PC7, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC7",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,7], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC8
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,8], xlab=PC1, ylab=PC8, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC8",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,8], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC9
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,9], xlab=PC1, ylab=PC9, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC9",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,9], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC10
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,10], xlab=PC1, ylab=PC10, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC10",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,10], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC11
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,11], xlab=PC1, ylab=PC11, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC11",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,11], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC12
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,12], xlab=PC1, ylab=PC12, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC12",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,12], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC13
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,13], xlab=PC1, ylab=PC13, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC13",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,13], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC14
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,14], xlab=PC1, ylab=PC14, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC14",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,14], replicate_labels, cex=.5, pos=1)

# ------------------------------------------------------------------------------------------------------------------
# PC1 vs PC15
# ------------------------------------------------------------------------------------------------------------------
plot(prcomp_results$x[,1], prcomp_results$x[,15], xlab=PC1, ylab=PC15, main="Read counts in all DHS sites (89,744 sites)\nPC1 vs PC15",
	bg=species_colors, col=species_colors, type="p", pch=plot_symbol, cex.main=1, cex.axis=.75, cex.lab=.75)

legend("top",c("human", "chimpanzee", "gorilla", "orangutan", "macaque"), pch=plot_symbol, cex=.75, pt.bg=c("orange", "black", "red", "blue", "green"), col=c("orange", "black", "red", "blue", "green"))

text(prcomp_results$x[,1], prcomp_results$x[,15], replicate_labels, cex=.5, pos=1)


# ------------------------------------------------------------------------------------------------------------------
# Finish up plots
# ------------------------------------------------------------------------------------------------------------------
cat("done\n\n") 
cat("Plots are located in pca.all_DHS_sites.pdf\n\n")

dev.off()

# ==================================================================================================================
# Finish up script
# ==================================================================================================================

cat("---------------------------------------------------------------------------------------\n")
cat("END:", date(), "\n")
cat("---------------------------------------------------------------------------------------\n")

sink()

