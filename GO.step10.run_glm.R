# ################################################################################################################################################################################################################
# GO.step10.run_glm.R
#
# CC-BY Lee Edsall
#	email: le49@duke.edu
#	Twitter: @LeeEdsall
#
# This script was used to analyze data for Edsall et al. 2019 which compared DNase-seq data from 5 primates (human, chimpanzee, gorilla, orangutan, macaque)
#
# There are 3 main analysis sections
#	1. Run DSS on all of the DHS sites to get the dispersion and normalization parameters
#	2. Iterate through the DHS sites to fit the glm model and perform the gateway test
#	3. Iterate through the differential DHS sites to perform the constraint tests and determine the type of change
#		Note: it is necessary to re-fit the glm in this step in order to extract the variance-covariance matrix
#		Sections 2 and 3 could be combined, but conceptually they make more sense when split into two sections
#
# There are 10 sections
# 	SECTION 1 OF 10: SETUP
# 	SECTION 2 OF 10: READ THE SCORE FILE AND PREPARE IT FOR THE GLM FUNCTION
# 	SECTION 3 OF 10: RUN DSS TO GET THE DISPERSION AND NORMALIZATION PARAMETERS
# 	SECTION 4 OF 10: GENERATE THE GLM AND PERFORM THE GATEWAY TEST
# 	SECTION 5 OF 10: PERFORM CONSTRAINT TESTS AND DETERMINE TYPE OF CHANGE
# 	SECTION 6 OF 10: UPDATE AND FINALIZE THE "results_all_sites" DATA FRAME
# 	SECTION 7 OF 10: CALCULATE TOTALS AND CREATE TABLES CONTAINING SUBSETS OF THE RESULTS
# 	SECTION 8 OF 10: WRITE RESULTS TO FILES
# 	SECTION 9 OF 10: CALCULATE CHECKSUM AND PRINT STATS
# 	SECTION 10 OF 10: FINISH UP
#
# Input file information
#	1. The replicates are in the columns and the locations are in the rows
#	2. The scores are integers and have not been normalized
#
# ################################################################################################################################################################################################################

# ################################################################################################################################################################################################################
# SECTION 1 OF 10: SETUP
#	1. load libraries
#	2. set global options
#	3. set global variables
#	4. open connection to log file
# ################################################################################################################################################################################################################

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DSS		Dispersion shrinkage for sequencing data
#		Citation:  Wu H, Wang C, Wu Z. 2013. A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data. Biostatistics 14(2):232-243.
#
# MASS		Contains the negative.binomial function needed for the GLM
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(DSS)
library(MASS)

options(scipen=0)
pvalue_threshold = .01

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Open the log file; write date and options
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sink(file="glm_analysis.log", append=FALSE, type="output", split=TRUE)
cat("-------------------------------------------------------------------------------------------------\n")
cat("BEGIN:", date(), "\n")
cat("-------------------------------------------------------------------------------------------------\n\n")

cat("Using the following settings\n")
cat("\tLog file: glm_analysis.log\n")
cat("\tp-value threshold:",  pvalue_threshold, "\n\n")

# ################################################################################################################################################################################################################
# END SECTION 1 OF 10: SETUP
# ################################################################################################################################################################################################################

# ################################################################################################################################################################################################################
# SECTION 2 OF 10: READ THE SCORE FILE AND PREPARE IT FOR THE GLM FUNCTION
#	Locations are in rows and scores are in columns
#	Column information
#		Column  1: chrom		Chromosome
#		Column  2: start		Start
#		Column  3: end			End
#		Column  4: h1			Human replicate 1 read counts
#		Column  5: h2			Human replicate 2 read counts
#		Column  6: h3			Human replicate 3 read counts
#		Column  7: c1			Chimpanzee replicate 1 read counts
#		Column  8: c2			Chimpanzee replicate 2 read counts
#		Column  9: c3			Chimpanzee replicate 3 read counts
#		Column 10: g1			Gorilla replicate 1 read counts
#		Column 11: g2			Gorilla replicate 2 read counts
#		Column 12: g3			Gorilla replicate 3 read counts
#		Column 13: o1			Orangutan replicate 1 read counts
#		Column 14: o2			Orangutan replicate 2 read counts
#		Column 15: o3			Orangutan replicate 3 read counts
#		Column 16: m1			Macaque replicate 1 read counts
#		Column 17: m2			Macaque replicate 2 read counts
#		Column 18: m3			Macaque replicate 3 read counts
#		Column 19: PLoS_overlap		Information about overlapping with the regions from Shibata et al. PLoS Genetics 2012
#
#	Create a data frame to use in the GLM.
#		1. Remove the location information
#		2. Transpose the data frame so the scores are in the rows and the locations are in the columns
# ################################################################################################################################################################################################################

# ================================================================================================================================================================================================================
# Read the score file
# ================================================================================================================================================================================================================

cat("Reading score file (glm_analysis.input_file.txt) ... \n")

scores_with_locations <- as.data.frame(read.table("glm_analysis.input_file.txt", sep="\t", header=FALSE))

colnames(scores_with_locations) <- c("chrom","start","end", "h1", "h2", "h3", "c1", "c2", "c3", "g1", "g2", "g3", "o1", "o2", "o3", "m1", "m2", "m3", "PLoS_overlap")

rownames(scores_with_locations) <- paste(scores_with_locations$chrom, scores_with_locations$start, scores_with_locations$end, scores_with_locations$PLoS_overlap, sep=":")

count_all_sites = nrow(scores_with_locations)

# ================================================================================================================================================================================================================
# Create a version of the score table to use in the GLM formula
#	Remove 4 columns containing non-score information: chrom, start, end, PLoS_overlap
#	Transpose the table so the locations are in the columns and the replicates are in the rows
#	Change the row names (which used to be column names) to use the Y nomenclature
# ================================================================================================================================================================================================================

scores <- t(scores_with_locations[,4:18])
rownames(scores) <- c("Yh1", "Yh2", "Yh3", "Yc1", "Yc2", "Yc3", "Yg1", "Yg2", "Yg3", "Yo1", "Yo2", "Yo3", "Ym1", "Ym2", "Ym3")

# ================================================================================================================================================================================================================
# Finish up
# ================================================================================================================================================================================================================

cat("done.\n\n")
cat("There are", count_all_sites, "sites.\n\n")

# ################################################################################################################################################################################################################
# END SECTION 2 OF 10: READ THE SCORE FILE AND PREPARE IT FOR THE GLM FUNCTION
# ################################################################################################################################################################################################################


# ################################################################################################################################################################################################################
# SECTION 3 OF 10: RUN DSS TO GET THE DISPERSION AND NORMALIZATION PARAMETERS
# ################################################################################################################################################################################################################

cat("Running DSS to get dispersion and normalization parameters ... \n")

# ================================================================================================================================================================================================================
# Set up the read counts matrix
#	Remove 4 columns containing non-score information: chrom, start, end, PLoS_overlap
#	Convert from a data.frame to a matrix
# ================================================================================================================================================================================================================

dss_scores <- as.matrix(scores_with_locations[,4:18])

# ================================================================================================================================================================================================================
# Set up the design matrix used by DSS
#	The species order is human, chimpanzee, gorilla, orangutan, macaque
#	The design matrix gets automatically sorted alphabetically, so to maintain the order we want, add numerals before the species name
# ================================================================================================================================================================================================================

dss_design=data.frame(species=c(rep("1human", 3), rep("2chimpanzee", 3), rep("3gorilla", 3), rep("4orangutan", 3), rep("0macaque", 3)))
rownames(dss_design) <- c("h1", "h2", "h3", "c1", "c2", "c3", "g1", "g2", "g3", "o1", "o2", "o3", "m1", "m2", "m3")
dss_X <- model.matrix(~species, data=dss_design)

# ================================================================================================================================================================================================================
# Run DSS
#       1. Create the seqData object
#       2. Calculate the normalization offsets (stored in library_size_offset_vector)
#       3. Estimate the dispersion parameter (stored in dispersion_parameter_vector)
# ================================================================================================================================================================================================================

dss_seqData=newSeqCountSet(dss_scores, as.data.frame(dss_X))
dss_seqData=estNormFactors(dss_seqData, method="total")
dss_seqData=estDispersion(dss_seqData)

library_size_offset_vector <- normalizationFactor(dss_seqData)
dispersion_parameter_vector <- dispersion(dss_seqData) 

cat("done.\n\n")

# ################################################################################################################################################################################################################
# END SECTION 3 OF 10: RUN DSS TO GET THE DISPERSION AND NORMALIZATION PARAMETERS
# ################################################################################################################################################################################################################


# ################################################################################################################################################################################################################
# SECTION 4 OF 10: GENERATE THE GLM AND PERFORM THE GATEWAY TEST
# Loop through the scores table and for each location, generate the GLM, perform the Chi-squared test, and store the results in a matrix called "results_all_sites"
#
# matrix results_all_sites
#	Locations are in rows and results are in columns
#	Column information
#		Column  1: width			DHS site width
#		Column  2: dispersion			Dispersion parameter (comes from DSS)
#		Column  3: offset_h1			Normalization factor for human replicate 1 (comes from DSS)
#		Column  4: offset_h2			Normalization factor for human replicate 2 (comes from DSS)
#		Column  5: offset_h3			Normalization factor for human replicate 3 (comes from DSS)
#		Column  6: offset_c1			Normalization factor for chimpanzee replicate 1 (comes from DSS)
#		Column  7: offset_c2			Normalization factor for chimpanzee replicate 2 (comes from DSS)
#		Column  8: offset_c3			Normalization factor for chimpanzee replicate 3 (comes from DSS)
#		Column  9: offset_g1			Normalization factor for gorilla replicate 1 (comes from DSS)
#		Column 10: offset_g2			Normalization factor for gorilla replicate 2 (comes from DSS)
#		Column 11: offset_g3			Normalization factor for gorilla replicate 3 (comes from DSS)
#		Column 12: offset_o1			Normalization factor for orangutan replicate 1 (comes from DSS)
#		Column 13: offset_o2			Normalization factor for orangutan replicate 2 (comes from DSS)
#		Column 14: offset_o3			Normalization factor for orangutan replicate 3 (comes from DSS)
#		Column 15: offset_m1			Normalization factor for macaque replicate 1 (comes from DSS)
#		Column 16: offset_m2			Normalization factor for macaque replicate 2 (comes from DSS)
#		Column 17: offset_m3			Normalization factor for macaque replicate 3 (comes from DSS)
#		Column 18: Bh				Human beta value (comes from the GLM)
#		Column 19: Bh_se			Human beta standard error (comes from the GLM)
#		Column 20: Bh_p				Human beta p-value (comes from the GLM)
#		Column 21: Bc				himpanzee beta value (comes from the GLM)
#		Column 22: Bc_se			Chimpanzee beta standard error (comes from the GLM)
#		Column 23: Bc_p				Chimpanzee beta p-value (comes from the GLM)
#		Column 24: Bg				Gorilla beta value (comes from the GLM)
#		Column 25: Bg_se			Gorilla beta standard error (comes from the GLM)
#		Column 26: Bg_p				Gorilla beta p-value (comes from the GLM)
#		Column 27: Bo				Orangutan beta value (comes from the GLM)
#		Column 28: Bo_se			Orangutan beta standard error (comes from the GLM)
#		Column 29: Bo_p				Orangutan beta p-value (comes from the GLM)
#		Column 30: Bm				Macaque beta value (comes from the GLM)
#		Column 31: Bm_se			Macaque beta standard error (comes from the GLM)
#		Column 32: Bm_p				Macaque beta p-value (comes from the GLM)
#		Column 33: gateway_original_pvalue	Gateway original p-value (the p-value for the GLM, calculated using a Chi-squared test)
#		Column 34: gateway_adjusted_pvalue	Gateway adjusted p-value (adjusted using the Benjamini-Hochberg correction)
#		Column 35: human_effect_size		Human effect size (beta value divided by standard error)
#		Column 36: chimpanzee_effect_size	Chimpanzee effect size (beta value divided by standard error)
#		Column 37: gorilla_effect_size		Gorilla effect size (beta value divided by standard error)
#		Column 38: orangutan_effect_size	Orangutan effect size (beta value divided by standard error)
#		Column 39: macaque_effect_size		Macaque effect size (beta value divided by standard error)
#		Column 40: h1				Human replicate 1 read counts (from input file)
#		Column 41: h2				Human replicate 2 read counts (from input file)
#		Column 42: h3				Human replicate 3 read counts (from input file)
#		Column 43: c1				Chimpanzee replicate 1 read counts (from input file)
#		Column 44: c2				Chimpanzee replicate 2 read counts (from input file)
#		Column 45: c3				Chimpanzee replicate 3 read counts (from input file)
#		Column 46: g1				Gorilla replicate 1 read counts (from input file)
#		Column 47: g2				Gorilla replicate 2 read counts (from input file)
#		Column 48: g3				Gorilla replicate 3 read counts (from input file)
#		Column 49: o1				Orangutan replicate 1 read counts (from input file)
#		Column 50: o2				Orangutan replicate 2 read counts (from input file)
#		Column 51: o3				Orangutan replicate 3 read counts (from input file)
#		Column 52: m1				Macaque replicate 1 read counts (from input file)
#		Column 53: m2				Macaque replicate 2 read counts (from input file)
#		Column 54: m3				Macaque replicate 3 read counts (from input file)
# ################################################################################################################################################################################################################

# ================================================================================================================================================================================================================
# Create the X vectors for the GLM function
#	The X vectors associate the samples with the species
#	Create one X vector for each non-macaque species
#	Because we're using macaque as the outgroup ("response element"), we don't need to create a species vector for it
# ================================================================================================================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Xh <- as.data.frame(c(1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0), nrow=15, ncol=1, byrow=TRUE)
rownames(Xh) <- c("h1", "h2", "h3", "c1", "c2", "c3", "g1", "g2", "g3", "o1", "o2", "o3", "m1", "m2", "m3")
colnames(Xh) <- "Xh"

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Xc <- as.data.frame(c(0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0), nrow=15, ncol=1, byrow=TRUE)
rownames(Xc) <- c("h1", "h2", "h3", "c1", "c2", "c3", "g1", "g2", "g3", "o1", "o2", "o3", "m1", "m2", "m3")
colnames(Xc) <- "Xc"

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Xg <- as.data.frame(c(0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0), nrow=15, ncol=1, byrow=TRUE)
rownames(Xg) <- c("h1", "h2", "h3", "c1", "c2", "c3", "g1", "g2", "g3", "o1", "o2", "o3", "m1", "m2", "m3")
colnames(Xg) <- "Xg"

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Xo <- as.data.frame(c(0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0), nrow=15, ncol=1, byrow=TRUE)
rownames(Xo) <- c("h1", "h2", "h3", "c1", "c2", "c3", "g1", "g2", "g3", "o1", "o2", "o3", "m1", "m2", "m3")
colnames(Xo) <- "Xo"

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Combine the individual X vectors into one matrix to use in the GLM
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
X <- as.matrix(cbind(Xh, Xc, Xg, Xo))

# ================================================================================================================================================================================================================
# Create the matrix that will hold the results
# ================================================================================================================================================================================================================

results_all_sites <- matrix(0, nrow=count_all_sites, ncol=54)
colnames(results_all_sites) = c(
	"width", "dispersion", 
	"offset_h1", "offset_h2", "offset_h3", "offset_c1", "offset_c2", "offset_c3", "offset_g1", "offset_g2", "offset_g3", "offset_o1", "offset_o2", "offset_o3", "offset_m1", "offset_m2", "offset_m3",
	"Bh", "Bh_se", "Bh_p", "Bc", "Bc_se", "Bc_p", "Bg", "Bg_se", "Bg_p", "Bo", "Bo_se", "Bo_p", "Bm", "Bm_se", "Bm_p",
	"gateway_original_pvalue", "gateway_adjusted_pvalue",
	"human_effect_size", "chimpanzee_effect_size", "gorilla_effect_size", "orangutan_effect_size", "macaque_effect_size",
	"h1", "h2", "h3", "c1", "c2", "c3", "g1", "g2", "g3", "o1", "o2", "o3", "m1", "m2", "m3")
rownames(results_all_sites) <- colnames(scores)

# ================================================================================================================================================================================================================
# Create a matrix that will hold the flag for the gateway test
#	Populate with "GNS" as the original value
# ================================================================================================================================================================================================================

flag_gateway_test <- matrix("GNS", nrow=count_all_sites, ncol=1)
colnames(flag_gateway_test) = "flag_gateway_test"
rownames(flag_gateway_test) <- colnames(scores)

# ================================================================================================================================================================================================================
# Create a matrix that will hold the flag for the type of change
#	Populate with "ssss" as the original value
# ================================================================================================================================================================================================================

flag_change_type <- matrix("ssss", nrow=count_all_sites, ncol=1)
colnames(flag_change_type) = "flag_change_type"
rownames(flag_change_type) <- colnames(scores)

# ================================================================================================================================================================================================================
# Create a matrix that will hold the p-values for the constraint tests
#	Populate with 1 for the original values
# ================================================================================================================================================================================================================

constraint_test_pvalues <- matrix(1, nrow=count_all_sites, ncol=15)
colnames(constraint_test_pvalues) <- c("ct_H_pvalue", "ct_C_pvalue", "ct_G_pvalue", "ct_O_pvalue", "ct_HCGO_pvalue", 
				       "ct_HC_pvalue", "ct_HG_pvalue", "ct_HO_pvalue", "ct_CG_pvalue", "ct_CO_pvalue", "ct_GO_pvalue",
				       "ct_HCG_pvalue", "ct_HCO_pvalue", "ct_HGO_pvalue", "ct_CGO_pvalue")
rownames(constraint_test_pvalues) <- colnames(scores)

# ================================================================================================================================================================================================================
# Start the loop
# ================================================================================================================================================================================================================

cat("Performing gateway test ... \n")
cat("\tGateway test: analyzing site 1 of", count_all_sites, "\n")

for (lcv in 1:count_all_sites)
{

# ================================================================================================================================================================================================================
# Print periodic status messages
# ================================================================================================================================================================================================================

	if ( (lcv %% 1000) == 0)
	{
		cat("\tGateway test: analyzing site", lcv, "of", count_all_sites, "\n")
	}

# ================================================================================================================================================================================================================
# Run the glm function and calculate the p-value
# ================================================================================================================================================================================================================

	score_data <- as.data.frame(cbind(scores[,lcv], X))
	colnames(score_data) <- c("Y", "Xh", "Xc", "Xg", "Xo")

	fit_glm <- glm(formula = Y ~ Xh + Xc + Xg + Xo + offset(log(library_size_offset_vector)), family = negative.binomial(theta = 1/dispersion_parameter_vector[lcv]), data=score_data)

	gateway_test_statistic <- fit_glm$null.deviance - fit_glm$deviance
	gateway_pvalue <- pchisq(gateway_test_statistic, df=4, lower=FALSE, log.p=FALSE)

# ================================================================================================================================================================================================================
# Update the results_all_sites matrix
# ================================================================================================================================================================================================================

	# DHS site width
	results_all_sites[lcv,1] = scores_with_locations[lcv,3] - scores_with_locations[lcv,2] + 1

	# Dispersion parameter
	results_all_sites[lcv,2] = dispersion_parameter_vector[lcv]

	# Normalization factors for human
	results_all_sites[lcv,3] = library_size_offset_vector[1]
	results_all_sites[lcv,4] = library_size_offset_vector[2]
	results_all_sites[lcv,5] = library_size_offset_vector[3]

	# Normalization factors for chimpanzee
	results_all_sites[lcv,6] = library_size_offset_vector[4]
	results_all_sites[lcv,7] = library_size_offset_vector[5]
	results_all_sites[lcv,8] = library_size_offset_vector[6]

	# Normalization factors for gorilla
	results_all_sites[lcv,9]  = library_size_offset_vector[7]
	results_all_sites[lcv,10] = library_size_offset_vector[8]
	results_all_sites[lcv,11] = library_size_offset_vector[9]

	# Normalization factors for orangutan
	results_all_sites[lcv,12] = library_size_offset_vector[10]
	results_all_sites[lcv,13] = library_size_offset_vector[11]
	results_all_sites[lcv,14] = library_size_offset_vector[12]

	# Normalization factors for macaque
	results_all_sites[lcv,15] = library_size_offset_vector[13]
	results_all_sites[lcv,16] = library_size_offset_vector[14]
	results_all_sites[lcv,17] = library_size_offset_vector[15]

	# Beta value, standard error, and p-value for human
	results_all_sites[lcv,18] = summary(fit_glm)$coefficients[2,1]
	results_all_sites[lcv,19] = summary(fit_glm)$coefficients[2,2]
	results_all_sites[lcv,20] = summary(fit_glm)$coefficients[2,4]

	# Beta value, standard error, and p-value for chimpanzee
	results_all_sites[lcv,21] = summary(fit_glm)$coefficients[3,1]
	results_all_sites[lcv,22] = summary(fit_glm)$coefficients[3,2]
	results_all_sites[lcv,23] = summary(fit_glm)$coefficients[3,4]

	# Beta value, standard error, and p-value for gorilla
	results_all_sites[lcv,24] = summary(fit_glm)$coefficients[4,1]
	results_all_sites[lcv,25] = summary(fit_glm)$coefficients[4,2]
	results_all_sites[lcv,26] = summary(fit_glm)$coefficients[4,4]

	# Beta value, standard error, and p-value for orangutan
	results_all_sites[lcv,27] = summary(fit_glm)$coefficients[5,1]
	results_all_sites[lcv,28] = summary(fit_glm)$coefficients[5,2]
	results_all_sites[lcv,29] = summary(fit_glm)$coefficients[5,4]

	# Beta value, standard error, and p-value for macaque
	results_all_sites[lcv,30] = summary(fit_glm)$coefficients[1,1]
	results_all_sites[lcv,31] = summary(fit_glm)$coefficients[1,2]
	results_all_sites[lcv,32] = summary(fit_glm)$coefficients[1,4]

	# Gateway p-value
	results_all_sites[lcv,33] = gateway_pvalue

	# Effect size for human (the beta value divided by the standard error)
	results_all_sites[lcv,35] = summary(fit_glm)$coefficients[2,1] / summary(fit_glm)$coefficients[2,2]

	# Effect size for chimpanzee (the beta value divided by the standard error)
	results_all_sites[lcv,36] = summary(fit_glm)$coefficients[3,1] / summary(fit_glm)$coefficients[3,2]

	# Effect size for gorilla (the beta value divided by the standard error)
	results_all_sites[lcv,37] = summary(fit_glm)$coefficients[4,1] / summary(fit_glm)$coefficients[4,2]

	# Effect size for orangutan (the beta value divided by the standard error)
	results_all_sites[lcv,38] = summary(fit_glm)$coefficients[5,1] / summary(fit_glm)$coefficients[5,2]

	# Effect size for macaque (the beta value divided by the standard error)
	results_all_sites[lcv,39] = summary(fit_glm)$coefficients[1,1] / summary(fit_glm)$coefficients[1,2]

	# Human scores
	results_all_sites[lcv,40] = scores[1,lcv]
	results_all_sites[lcv,41] = scores[2,lcv]
	results_all_sites[lcv,42] = scores[3,lcv]

	# Chimpanzee scores
	results_all_sites[lcv,43] = scores[4,lcv]
	results_all_sites[lcv,44] = scores[5,lcv]
	results_all_sites[lcv,45] = scores[6,lcv]

	# Gorilla scores
	results_all_sites[lcv,46] = scores[7,lcv]
	results_all_sites[lcv,47] = scores[8,lcv]
	results_all_sites[lcv,48] = scores[9,lcv]

	# Orangutan scores
	results_all_sites[lcv,49] = scores[10,lcv]
	results_all_sites[lcv,50] = scores[11,lcv]
	results_all_sites[lcv,51] = scores[12,lcv]

	# Macaque scores
	results_all_sites[lcv,52] = scores[13,lcv]
	results_all_sites[lcv,53] = scores[14,lcv]
	results_all_sites[lcv,54] = scores[15,lcv]

}  # END: for (lcv in 1:count_all_sites)

# ================================================================================================================================================================================================================
# Calculate the adjusted p-value using the Benjamini-Hochberg correction
# ================================================================================================================================================================================================================

results_all_sites[,34] <- p.adjust(results_all_sites[,33], "BH")

# ================================================================================================================================================================================================================
# Check for significance and update the flag for the gateway test
# ================================================================================================================================================================================================================

for (lcv in 1:count_all_sites)
{
	if (results_all_sites[lcv,34] < pvalue_threshold)
	{
		flag_gateway_test[lcv,1] = "GD"
	}
}

# ================================================================================================================================================================================================================
# Finish up
# ================================================================================================================================================================================================================

cat("done.\n\n")

# ################################################################################################################################################################################################################
# END SECTION 4 OF 10: GENERATE THE GLM AND PERFORM THE GATEWAY TEST
# ################################################################################################################################################################################################################


# ################################################################################################################################################################################################################
# SECTION 5 OF 10: PERFORM CONSTRAINT TESTS AND DETERMINE TYPE OF CHANGE
# ################################################################################################################################################################################################################

# ================================================================================================================================================================================================================
# Create the constraint matrices
# ================================================================================================================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Changes in one species
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human
H_matrix <- matrix(c(1, -1/4, -1/4, -1/4), nrow=1, ncol=4, byrow = TRUE)
colnames(H_matrix) = c("h","c","g","o")

# Chimpanzee
C_matrix <- matrix(c(-1/4, 1, -1/4, -1/4), nrow=1, ncol=4, byrow = TRUE)
colnames(C_matrix) = c("h","c","g","o")

# Gorilla
G_matrix <- matrix(c(-1/4, -1/4, 1, -1/4), nrow=1, ncol=4, byrow = TRUE)
colnames(G_matrix) = c("h","c","g","o")

# Orangutan
O_matrix <- matrix(c(-1/4, -1/4, -1/4, 1), nrow=1, ncol=4, byrow = TRUE)
colnames(O_matrix) = c("h","c","g","o")

# Macaque (also identifies changes in human/chimpanzee/gorilla/orangutan)
HCGO_matrix <- matrix(c(-1/4, -1/4, -1/4, -1/4), nrow=1, ncol=4, byrow = TRUE)
colnames(HCGO_matrix) = c("h","c","g","o")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Changes in two species
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human/chimpanzee
HC_matrix <- matrix(c(1/2, 1/2, -1/3, -1/3), nrow=1, ncol=4, byrow = TRUE)
colnames(HC_matrix) = c("h","c","g","o")

# Human/gorilla
HG_matrix <- matrix(c(1/2, -1/3, 1/2, -1/3), nrow=1, ncol=4, byrow = TRUE)
colnames(HG_matrix) = c("h","c","g","o")

# Human/orangutan
HO_matrix <- matrix(c(1/2, -1/3, -1/3, 1/2), nrow=1, ncol=4, byrow = TRUE)
colnames(HO_matrix) = c("h","c","g","o")

# Chimpanzee/gorilla
CG_matrix <- matrix(c(-1/3, 1/2, 1/2, -1/3), nrow=1, ncol=4, byrow = TRUE)
colnames(CG_matrix) = c("h","c","g","o")

# Chimpanzee/orangutan
CO_matrix <- matrix(c(-1/3, 1/2, -1/3, 1/2), nrow=1, ncol=4, byrow = TRUE)
colnames(CO_matrix) = c("h","c","g","o")

# Gorilla/orangutan
GO_matrix <- matrix(c(-1/3, -1/3, 1/2, 1/2), nrow=1, ncol=4, byrow = TRUE)
colnames(GO_matrix) = c("h","c","g","o")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Changes in three species
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human/chimpanzee/gorilla
HCG_matrix <- matrix(c(1/3, 1/3, 1/3, -1/2), nrow=1, ncol=4, byrow = TRUE)
colnames(HCG_matrix) = c("h","c","g","o")

# Human/chimpanzee/orangutan
HCO_matrix <- matrix(c(1/3, 1/3, -1/2, 1/3), nrow=1, ncol=4, byrow = TRUE)
colnames(HCO_matrix) = c("h","c","g","o")

# Human/gorilla/orangutan
HGO_matrix <- matrix(c(1/3, -1/2, 1/3, 1/3), nrow=1, ncol=4, byrow = TRUE)
colnames(HGO_matrix) = c("h","c","g","o")

# Chimpanzee/gorilla/orangutan
CGO_matrix <- matrix(c(-1/2, 1/3, 1/3, 1/3), nrow=1, ncol=4, byrow = TRUE)
colnames(CGO_matrix) = c("h","c","g","o")

# ================================================================================================================================================================================================================
# Start the loop
# ================================================================================================================================================================================================================

cat("Performing constraint tests to determine type of change ... \n")
cat("\tConstraint tests: analyzing site 1 of", count_all_sites, "\n")

for (lcv in 1:count_all_sites)
{

# ================================================================================================================================================================================================================
# Print periodic status messages
# ================================================================================================================================================================================================================

	if ( (lcv %% 1000) == 0)
	{
		cat("\tConstraint tests: analyzing site", lcv, "of", count_all_sites, "\n")
	}

# ================================================================================================================================================================================================================
# If the site is not significantly differential, skip it
# ================================================================================================================================================================================================================

	if (flag_gateway_test[lcv,1] == "GNS")
	{
		next
	}

# ================================================================================================================================================================================================================
# Re-run the glm to extract the beta values and the variance-covariance matrix
# ================================================================================================================================================================================================================

	score_data <- as.data.frame(cbind(scores[,lcv], X))
	colnames(score_data) <- c("Y", "Xh", "Xc", "Xg", "Xo")

	fit_glm <- glm(formula = Y ~ Xh + Xc + Xg + Xo + offset(log(library_size_offset_vector)), family = negative.binomial(theta = 1/dispersion_parameter_vector[lcv]), data=score_data)

	vcov_matrix <- vcov(fit_glm)[2:5,2:5]
	beta_matrix <- summary(fit_glm)$coefficients[2:5,1]

# ================================================================================================================================================================================================================
# Calculate the constraint values
# ================================================================================================================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Changes in one species
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Human change
	H_pvalue <- pchisq( ( t(H_matrix %*% beta_matrix) %*% solve(H_matrix %*% vcov_matrix %*% t(H_matrix)) %*% (H_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Chimpanzee change
	C_pvalue <- pchisq( ( t(C_matrix %*% beta_matrix) %*% solve(C_matrix %*% vcov_matrix %*% t(C_matrix)) %*% (C_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Gorilla change
	G_pvalue <- pchisq( ( t(G_matrix %*% beta_matrix) %*% solve(G_matrix %*% vcov_matrix %*% t(G_matrix)) %*% (G_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Orangutan change
	O_pvalue <- pchisq( ( t(O_matrix %*% beta_matrix) %*% solve(O_matrix %*% vcov_matrix %*% t(O_matrix)) %*% (O_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Macaque change (also identifies changes in human/chimpanzee/gorilla/orangutan)
	HCGO_pvalue <- pchisq( ( t(HCGO_matrix %*% beta_matrix) %*% solve(HCGO_matrix %*% vcov_matrix %*% t(HCGO_matrix)) %*% (HCGO_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Changes in two species
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Human/chimpanzee change
	HC_pvalue <- pchisq( ( t(HC_matrix %*% beta_matrix) %*% solve(HC_matrix %*% vcov_matrix %*% t(HC_matrix)) %*% (HC_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Human/gorilla change
	HG_pvalue <- pchisq( ( t(HG_matrix %*% beta_matrix) %*% solve(HG_matrix %*% vcov_matrix %*% t(HG_matrix)) %*% (HG_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Human/orangutan change
	HO_pvalue <- pchisq( ( t(HO_matrix %*% beta_matrix) %*% solve(HO_matrix %*% vcov_matrix %*% t(HO_matrix)) %*% (HO_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Chimpanzee/gorilla change
	CG_pvalue <- pchisq( ( t(CG_matrix %*% beta_matrix) %*% solve(CG_matrix %*% vcov_matrix %*% t(CG_matrix)) %*% (CG_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Chimpanzee/orangutan change
	CO_pvalue <- pchisq( ( t(CO_matrix %*% beta_matrix) %*% solve(CO_matrix %*% vcov_matrix %*% t(CO_matrix)) %*% (CO_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Gorilla/orangutan change
	GO_pvalue <- pchisq( ( t(GO_matrix %*% beta_matrix) %*% solve(GO_matrix %*% vcov_matrix %*% t(GO_matrix)) %*% (GO_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Changes in three species
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Human/chimpanzee/gorilla change
	HCG_pvalue <- pchisq( ( t(HCG_matrix %*% beta_matrix) %*% solve(HCG_matrix %*% vcov_matrix %*% t(HCG_matrix)) %*% (HCG_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Human/chimpanzee/orangutan change
	HCO_pvalue <- pchisq( ( t(HCO_matrix %*% beta_matrix) %*% solve(HCO_matrix %*% vcov_matrix %*% t(HCO_matrix)) %*% (HCO_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Human/gorilla/orangutan change
	HGO_pvalue <- pchisq( ( t(HGO_matrix %*% beta_matrix) %*% solve(HGO_matrix %*% vcov_matrix %*% t(HGO_matrix)) %*% (HGO_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

	# Chimpanzee/gorilla/orangutan change
	CGO_pvalue <- pchisq( ( t(CGO_matrix %*% beta_matrix) %*% solve(CGO_matrix %*% vcov_matrix %*% t(CGO_matrix)) %*% (CGO_matrix %*% beta_matrix) ), df=1, lower=FALSE, log.p=FALSE)

# ================================================================================================================================================================================================================
# Update the constraint test p-values matrix
# ================================================================================================================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Changes in one species
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Human change
	constraint_test_pvalues[lcv,1] <- H_pvalue

	# Chimpanzee change
	constraint_test_pvalues[lcv,2] <- C_pvalue

	# Gorilla change
	constraint_test_pvalues[lcv,3] <- G_pvalue

	# Orangutan change
	constraint_test_pvalues[lcv,4] <- O_pvalue

	# Macaque change (also identifies human/chimpanzee/gorilla/orangutan change)
	constraint_test_pvalues[lcv,5] <- HCGO_pvalue

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Changes in two species
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Human/chimpanzee change
	constraint_test_pvalues[lcv,6] <- HC_pvalue

	# Human/gorilla change
	constraint_test_pvalues[lcv,7] <- HG_pvalue

	# Human/orangutan change
	constraint_test_pvalues[lcv,8] <- HO_pvalue

	# Chimpanzee/gorilla change
	constraint_test_pvalues[lcv,9] <- CG_pvalue

	# Chimpanzee/orangutan change
	constraint_test_pvalues[lcv,10] <- CO_pvalue

	# Gorilla/orangutan change
	constraint_test_pvalues[lcv,11] <- GO_pvalue

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Changes in three species
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Human/chimpanzee/gorilla change
	constraint_test_pvalues[lcv,12] <- HCG_pvalue

	# Human/chimpanzee/orangutan change
	constraint_test_pvalues[lcv,13] <- HCO_pvalue

	# Human/gorilla/orangutan change
	constraint_test_pvalues[lcv,14] <- HGO_pvalue

	# Chimpanzee/gorilla/orangutan change
	constraint_test_pvalues[lcv,15] <- CGO_pvalue

# ================================================================================================================================================================================================================
# Identify the smallest p-value and update the flag for the type of change (flag_change_type)
# ================================================================================================================================================================================================================

	minimum_pvalue <- min(H_pvalue, C_pvalue, G_pvalue, O_pvalue, HCGO_pvalue, HC_pvalue, HG_pvalue, HO_pvalue, CG_pvalue, CO_pvalue, GO_pvalue, HCG_pvalue, HCO_pvalue, HGO_pvalue, CGO_pvalue)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 1 of 16: change in one species - human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (H_pvalue < pvalue_threshold) & (H_pvalue == minimum_pvalue) )
	{
		if (results_all_sites[lcv,18] > 0)
		{
			flag_change_type[lcv,1] <- "Gsss"
		} else {
			flag_change_type[lcv,1] <- "Lsss"
		}
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 2 of 16: change in one species - chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (C_pvalue < pvalue_threshold) & (C_pvalue == minimum_pvalue) )
	{
		if (results_all_sites[lcv,21] > 0)
		{
			flag_change_type[lcv,1] <- "sGss"
		} else {
			flag_change_type[lcv,1] <- "sLss"
		}
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 3 of 16: change in one species - gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (G_pvalue < pvalue_threshold) & (G_pvalue == minimum_pvalue) )
	{
		if (results_all_sites[lcv,24] > 0)
		{
			flag_change_type[lcv,1] <- "ssGs"
		} else {
			flag_change_type[lcv,1] <- "ssLs"
		}	
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 4 of 16: change in one species - orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (O_pvalue < pvalue_threshold) & (O_pvalue == minimum_pvalue) )
	{
		if (results_all_sites[lcv,27] > 0)
		{
			flag_change_type[lcv,1] <- "sssG"
		} else {
			flag_change_type[lcv,1] <- "sssL"
		}
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 5 of 16: change in one species - macaque (also identifies changes in human/chimpanzee/gorilla/orangutan)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (HCGO_pvalue < pvalue_threshold) & (HCGO_pvalue == minimum_pvalue) )
	{
		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "LLLL"
		} else {
			flag_change_type[lcv,1] <- "GGGG"
		}
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 6 of 16: change in two species - human/chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (HC_pvalue < pvalue_threshold) & (HC_pvalue == minimum_pvalue) )
	{
		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,21] > 0) )
		{
			flag_change_type[lcv,1] <- "GGss"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,21] < 0) )
		{
			flag_change_type[lcv,1] <- "GLss"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] > 0) )
		{
			flag_change_type[lcv,1] <- "LGss"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] < 0) )
		{
			flag_change_type[lcv,1] <- "LLss"
		}
	}  # END: if ( (HC_pvalue < pvalue_threshold) & (HC_pvalue == minimum_pvalue) )

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 7 of 16: change in two species - human/gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (HG_pvalue < pvalue_threshold) & (HG_pvalue == minimum_pvalue) )
	{
		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,24] > 0) )
		{
			flag_change_type[lcv,1] <- "GsGs"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,24] < 0) )
		{
			flag_change_type[lcv,1] <- "GsLs"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,24] > 0) )
		{
			flag_change_type[lcv,1] <- "LsGs"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,24] < 0) )
		{
			flag_change_type[lcv,1] <- "LsLs"
		}
	}  # END: if ( (HG_pvalue < pvalue_threshold) & (HG_pvalue == minimum_pvalue) )

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 8 of 16: change in two species - human/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (HO_pvalue < pvalue_threshold) & (HO_pvalue == minimum_pvalue) )
	{
		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "GssG"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "GssL"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "LssG"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "LssL"
		}
	}  # END: if ( (HO_pvalue < pvalue_threshold) & (HO_pvalue == minimum_pvalue) )

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 9 of 16: change in two species - chimpanzee/gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (CG_pvalue < pvalue_threshold) & (CG_pvalue == minimum_pvalue) )
	{
		if ( (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,24] > 0) )
		{
			flag_change_type[lcv,1] <- "sGGs"
		}

		if ( (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,24] < 0) )
		{
			flag_change_type[lcv,1] <- "sGLs"
		}

		if ( (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] > 0) )
		{
			flag_change_type[lcv,1] <- "sLGs"
		}

		if ( (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] < 0) )
		{
			flag_change_type[lcv,1] <- "sLLs"
		}
	}  # END: if ( (CG_pvalue < pvalue_threshold) & (CG_pvalue == minimum_pvalue) )

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 10 of 16: change in two species - chimpanzee/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (CO_pvalue < pvalue_threshold) & (CO_pvalue == minimum_pvalue) )
	{
		if ( (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "sGsG"
		}

		if ( (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "sGsL"
		}

		if ( (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "sLsG"
		}

		if ( (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "sLsL"
		}
	}  # END: if ( (CO_pvalue < pvalue_threshold) & (CO_pvalue == minimum_pvalue) )

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 11 of 16: change in two species - gorilla/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (GO_pvalue < pvalue_threshold) & (GO_pvalue == minimum_pvalue) )
	{
		if ( (results_all_sites[lcv,24] > 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "ssGG"
		}

		if ( (results_all_sites[lcv,24] > 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "ssGL"
		}

		if ( (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "ssLG"
		}

		if ( (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "ssLL"
		}
	}  # END: if ( (GO_pvalue < pvalue_threshold) & (GO_pvalue == minimum_pvalue) )

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 12 of 16: change in three species - human/chimpanzee/gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (HCG_pvalue < pvalue_threshold) & (HCG_pvalue == minimum_pvalue) )
	{

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,24] > 0) )
		{
			flag_change_type[lcv,1] <- "GGGs"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,24] < 0) )
		{
			flag_change_type[lcv,1] <- "GGLs"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] > 0) )
		{
			flag_change_type[lcv,1] <- "GLGs"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] < 0) )
		{
			flag_change_type[lcv,1] <- "GLLs"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,24] > 0) )
		{
			flag_change_type[lcv,1] <- "LGGs"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,24] < 0) )
		{
			flag_change_type[lcv,1] <- "LGLs"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] > 0) )
		{
			flag_change_type[lcv,1] <- "LLGs"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] < 0) )
		{
			flag_change_type[lcv,1] <- "LLLs"
		}

	}  # END: if ( (HCG_pvalue < pvalue_threshold) & (HCG_pvalue == minimum_pvalue) )

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 13 of 16: change in three species - human/chimpanzee/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (HCO_pvalue < pvalue_threshold) & (HCO_pvalue == minimum_pvalue) )
	{

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "GGsG"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "GGsL"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "GLsG"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "GLsL"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "LGsG"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "LGsL"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "LLsG"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "LLsL"
		}

	}  # END: if ( (HCO_pvalue < pvalue_threshold) & (HCO_pvalue == minimum_pvalue) )


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 14 of 16: change in three species - human/gorilla/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (HGO_pvalue < pvalue_threshold) & (HGO_pvalue == minimum_pvalue) )
	{

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,24] > 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "GsGG"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,24] > 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "GsGL"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "GsLG"
		}

		if ( (results_all_sites[lcv,18] > 0) & (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "GsLL"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,24] > 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "LsGG"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,24] > 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "LsGL"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "LsLG"
		}

		if ( (results_all_sites[lcv,18] < 0) & (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "LsLL"
		}

	}  # END: if ( (HGO_pvalue < pvalue_threshold) & (HGO_pvalue == minimum_pvalue) )

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 15 of 16: change in three species - chimpanzee/gorilla/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (CGO_pvalue < pvalue_threshold) & (CGO_pvalue == minimum_pvalue) )
	{

		if ( (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,24] > 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "sGGG"
		}

		if ( (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,24] > 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "sGGL"
		}

		if ( (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "sGLG"
		}

		if ( (results_all_sites[lcv,21] > 0) & (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "sGLL"
		}

		if ( (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] > 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "sLGG"
		}

		if ( (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] > 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "sLGL"
		}

		if ( (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] > 0) )
		{
			flag_change_type[lcv,1] <- "sLLG"
		}

		if ( (results_all_sites[lcv,21] < 0) & (results_all_sites[lcv,24] < 0) & (results_all_sites[lcv,27] < 0) )
		{
			flag_change_type[lcv,1] <- "sLLL"
		}

	}  # END: if ( (CGO_pvalue < pvalue_threshold) & (CGO_pvalue == minimum_pvalue) )


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential subset 16 of 16 - other changes
#	If the flag is still "ssss" then none of the other subsets had a p-value below the threshold
#	Change the flag to "OOOO"
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (flag_change_type[lcv,1] == "ssss")
	{
		flag_change_type[lcv,1] <- "OOOO"
	}


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}  # END: for (lcv in 1:count_all_sites)

# ================================================================================================================================================================================================================
# Finish up
# ================================================================================================================================================================================================================

cat("done.\n\n")

# ################################################################################################################################################################################################################
# END SECTION 5 OF 10: PERFORM CONSTRAINT TESTS AND DETERMINE TYPE OF CHANGE
# ################################################################################################################################################################################################################


# ################################################################################################################################################################################################################
# SECTION 6 OF 10: UPDATE AND FINALIZE THE "results_all_sites" DATA FRAME
#	Add the location and PLoS overlap columns
#	Add a column to hold the flag for the gateway test
#	Add a column to hold the flag for the type of change
#	Add the columns to hold the p-values for the constraint tests
#	Column information
#		Column  1: chrom			Chromosome
#		Column  2: start			Start
#		Column  3: end				End
#		Column  4: PLoS_overlap			Information about overlapping with the regions from Shibata et al. PLoS Genetics 2012
#		Column  5: width			DHS site width
#		Column  6: dispersion			Dispersion parameter (comes from DSS)
#		Column  7: offset_h1			Normalization factor for human replicate 1 (comes from DSS)
#		Column  8: offset_h2			Normalization factor for human replicate 2 (comes from DSS)
#		Column  9: offset_h3			Normalization factor for human replicate 3 (comes from DSS)
#		Column 10: offset_c1			Normalization factor for chimpanzee replicate 1 (comes from DSS)
#		Column 11: offset_c2			Normalization factor for chimpanzee replicate 2 (comes from DSS)
#		Column 12: offset_c3			Normalization factor for chimpanzee replicate 3 (comes from DSS)
#		Column 13: offset_g1			Normalization factor for gorilla replicate 1 (comes from DSS)
#		Column 14: offset_g2			Normalization factor for gorilla replicate 2 (comes from DSS)
#		Column 15: offset_g3			Normalization factor for gorilla replicate 3 (comes from DSS)
#		Column 16: offset_o1			Normalization factor for orangutan replicate 1 (comes from DSS)
#		Column 17: offset_o2			Normalization factor for orangutan replicate 2 (comes from DSS)
#		Column 18: offset_o3			Normalization factor for orangutan replicate 3 (comes from DSS)
#		Column 19: offset_m1			Normalization factor for macaque replicate 1 (comes from DSS)
#		Column 20: offset_m2			Normalization factor for macaque replicate 2 (comes from DSS)
#		Column 21: offset_m3			Normalization factor for macaque replicate 3 (comes from DSS)
#		Column 22: Bh				Human beta value (comes from the GLM)
#		Column 23: Bh_se			Human beta standard error (comes from the GLM)
#		Column 24: Bh_p				Human beta p-value (comes from the GLM)
#		Column 25: Bc				himpanzee beta value (comes from the GLM)
#		Column 26: Bc_se			Chimpanzee beta standard error (comes from the GLM)
#		Column 27: Bc_p				Chimpanzee beta p-value (comes from the GLM)
#		Column 28: Bg				Gorilla beta value (comes from the GLM)
#		Column 29: Bg_se			Gorilla beta standard error (comes from the GLM)
#		Column 30: Bg_p				Gorilla beta p-value (comes from the GLM)
#		Column 31: Bo				Orangutan beta value (comes from the GLM)
#		Column 32: Bo_se			Orangutan beta standard error (comes from the GLM)
#		Column 33: Bo_p				Orangutan beta p-value (comes from the GLM)
#		Column 34: Bm				Macaque beta value (comes from the GLM)
#		Column 35: Bm_se			Macaque beta standard error (comes from the GLM)
#		Column 36: Bm_p				Macaque beta p-value (comes from the GLM)
#		Column 37: gateway_original_pvalue	Gateway original p-value (the p-value for the GLM, calculated using a Chi-squared test)
#		Column 38: gateway_adjusted_pvalue	Gateway adjusted p-value (adjusted using the Benjamini-Hochberg correction)
#		Column 39: human_effect_size		Human effect size (beta value divided by standard error)
#		Column 40: chimpanzee_effect_size	Chimpanzee effect size (beta value divided by standard error)
#		Column 41: gorilla_effect_size		Gorilla effect size (beta value divided by standard error)
#		Column 42: orangutan_effect_size	Orangutan effect size (beta value divided by standard error)
#		Column 43: macaque_effect_size		Macaque effect size (beta value divided by standard error)
#		Column 44: h1				Human replicate 1 read counts (from input file)
#		Column 45: h2				Human replicate 2 read counts (from input file)
#		Column 46: h3				Human replicate 3 read counts (from input file)
#		Column 47: c1				Chimpanzee replicate 1 read counts (from input file)
#		Column 48: c2				Chimpanzee replicate 2 read counts (from input file)
#		Column 49: c3				Chimpanzee replicate 3 read counts (from input file)
#		Column 50: g1				Gorilla replicate 1 read counts (from input file)
#		Column 51: g2				Gorilla replicate 2 read counts (from input file)
#		Column 52: g3				Gorilla replicate 3 read counts (from input file)
#		Column 53: o1				Orangutan replicate 1 read counts (from input file)
#		Column 54: o2				Orangutan replicate 2 read counts (from input file)
#		Column 55: o3				Orangutan replicate 3 read counts (from input file)
#		Column 56: m1				Macaque replicate 1 read counts (from input file)
#		Column 57: m2				Macaque replicate 2 read counts (from input file)
#		Column 58: m3				Macaque replicate 3 read counts (from input file)
#		Column 59: flag_gateway_test		Flag for the gateway test (GNS for similar sites; GD for differential sites)
#		Column 60: flag_change_type		Flag for the type of change
#		Column 61: ct_H_pvalue			constraint test p-value for change in human
#		Column 62: ct_C_pvalue			constraint test p-value for change in chimpanzee
#		Column 63: ct_G_pvalue			constraint test p-value for change in gorilla
#		Column 64: ct_O_pvalue			constraint test p-value for change in orangutan
#		Column 65: ct_HCGO_pvalue		constraint test p-value for change in macaque (also identifies changes in human/chimpanzee/gorilla/orangutan)
#		Column 66: ct_HC_pvalue			constraint test p-value for change in human/chimpanzee
#		Column 67: ct_HG_pvalue			constraint test p-value for change in human/gorilla
#		Column 68: ct_HO_pvalue			constraint test p-value for change in human/orangutan
#		Column 69: ct_CG_pvalue			constraint test p-value for change in chimpanzee/gorilla
#		Column 70: ct_CO_pvalue			constraint test p-value for change in chimpanzee/orangutan
#		Column 71: ct_GO_pvalue			constraint test p-value for change in gorilla/orangutan
#		Column 72: ct_HCG_pvalue		constraint test p-value for change in human/chimpanzee/gorilla
#		Column 73: ct_HCO_pvalue		constraint test p-value for change in human/chimpanzee/orangutan
#		Column 74: ct_HGO_pvalue		constraint test p-value for change in human/gorilla/orangutan
#		Column 75: ct_CGO_pvalue		constraint test p-value for change in chimpanzee/gorilla/orangutan
# ################################################################################################################################################################################################################

cat("Finalizing tables ... \n")

# ================================================================================================================================================================================================================
# Extract the location and PLoS overlap information and store in a data frame
# ================================================================================================================================================================================================================

locations_all_sites <- as.data.frame(do.call(rbind,strsplit(rownames(results_all_sites),":")))
colnames(locations_all_sites) <- c("chrom","start","end", "PLoS_overlap")

# ================================================================================================================================================================================================================
# Update the results_all_sites information and store as a data frame
#	cbind: location information; results_all_sites; flag for the gateway test; flag for the type of change; constraint test p-values
# ================================================================================================================================================================================================================

results_all_sites <- cbind.data.frame(locations_all_sites, results_all_sites, flag_gateway_test, flag_change_type, constraint_test_pvalues)

# ================================================================================================================================================================================================================
# Finish up
# ================================================================================================================================================================================================================

cat("done.\n\n")

# ################################################################################################################################################################################################################
# END SECTION 6 OF 10: UPDATE AND FINALIZE THE "results_all_sites" DATA FRAME
# ################################################################################################################################################################################################################


# ################################################################################################################################################################################################################
# SECTION 7 OF 10: CALCULATE TOTALS AND CREATE TABLES CONTAINING SUBSETS OF THE RESULTS
#	For each subset, create 2 tables: one containing all of the information; and one containing only the location information
#	For the table containing all of the information, sort it based on the absolute value of the effect sizes
#	For the table containing only the location information, sort it based on the chromosome and start position
# ################################################################################################################################################################################################################

cat("Calculating totals and creating subsets of the results ...\n")

# ================================================================================================================================================================================================================
# Not significantly differential sites
# ================================================================================================================================================================================================================

results_not_significant_sites <- subset(results_all_sites, flag_gateway_test == "GNS")
results_not_significant_sites <- results_not_significant_sites[with(results_not_significant_sites, order(chrom, as.numeric(as.character(start)))),]

locations_not_significant_sites <- results_not_significant_sites[,1:3]

count_not_significant_sites = nrow(results_not_significant_sites)

# ================================================================================================================================================================================================================
# Differential sites
# ================================================================================================================================================================================================================

results_differential_sites <- subset(results_all_sites, flag_gateway_test == "GD")
results_differential_sites <- results_differential_sites[with(results_differential_sites, order(chrom, as.numeric(as.character(start)))),]

locations_differential_sites <- results_differential_sites[,1:3]

count_differential_sites = nrow(results_differential_sites)

# ================================================================================================================================================================================================================
# Differential subset 1 of 16:  changes in one species - human
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_H_all <- subset(results_all_sites, flag_change_type %in% c("Gsss","Lsss"))
results_H_all <- results_H_all[order(abs(results_H_all$human_effect_size), decreasing=TRUE),]

locations_H_all <- results_H_all[with(results_H_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_H_all = nrow(results_H_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_H_increases <- subset(results_all_sites, flag_change_type == "Gsss")
results_H_increases <- results_H_increases[order(results_H_increases$human_effect_size, decreasing=TRUE),]

locations_H_increases <- results_H_increases[with(results_H_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_H_increases = nrow(results_H_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_H_decreases <- subset(results_all_sites, flag_change_type == "Lsss")
results_H_decreases <- results_H_decreases[order(results_H_decreases$human_effect_size, decreasing=FALSE),]

locations_H_decreases <- results_H_decreases[with(results_H_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_H_decreases = nrow(results_H_decreases)

# ================================================================================================================================================================================================================
# Differential subset 2 of 16:  changes in one species - chimpanzee
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_C_all <- subset(results_all_sites, flag_change_type %in% c("sGss","sLss"))
results_C_all <- results_C_all[order(abs(results_C_all$chimpanzee_effect_size), decreasing=TRUE),]

locations_C_all <- results_C_all[with(results_C_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_C_all = nrow(results_C_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_C_increases <- subset(results_all_sites, flag_change_type == "sGss")
results_C_increases <- results_C_increases[order(results_C_increases$chimpanzee_effect_size, decreasing=TRUE),]

locations_C_increases <- results_C_increases[with(results_C_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_C_increases = nrow(results_C_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_C_decreases <- subset(results_all_sites, flag_change_type == "sLss")
results_C_decreases <- results_C_decreases[order(results_C_decreases$chimpanzee_effect_size, decreasing=FALSE),]

locations_C_decreases <- results_C_decreases[with(results_C_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_C_decreases = nrow(results_C_decreases)

# ================================================================================================================================================================================================================
# Differential subset 3 of 16:  changes in one species - gorilla
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_G_all <- subset(results_all_sites, flag_change_type %in% c("ssGs","ssLs"))
results_G_all <- results_G_all[order(abs(results_G_all$gorilla_effect_size), decreasing=TRUE),]

locations_G_all <- results_G_all[with(results_G_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_G_all = nrow(results_G_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_G_increases <- subset(results_all_sites, flag_change_type == "ssGs")
results_G_increases <- results_G_increases[order(results_G_increases$gorilla_effect_size, decreasing=TRUE),]

locations_G_increases <- results_G_increases[with(results_G_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_G_increases = nrow(results_G_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_G_decreases <- subset(results_all_sites, flag_change_type == "ssLs")
results_G_decreases <- results_G_decreases[order(results_G_decreases$gorilla_effect_size, decreasing=FALSE),]

locations_G_decreases <- results_G_decreases[with(results_G_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_G_decreases = nrow(results_G_decreases)

# ================================================================================================================================================================================================================
# Differential subset 4 of 16:  changes in one species - orangutan
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_O_all <- subset(results_all_sites, flag_change_type %in% c("sssG","sssL"))
results_O_all <- results_O_all[order(abs(results_O_all$orangutan_effect_size), decreasing=TRUE),]

locations_O_all <- results_O_all[with(results_O_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_O_all = nrow(results_O_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_O_increases <- subset(results_all_sites, flag_change_type == "sssG")
results_O_increases <- results_O_increases[order(results_O_increases$orangutan_effect_size, decreasing=TRUE),]

locations_O_increases <- results_O_increases[with(results_O_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_O_increases = nrow(results_O_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_O_decreases <- subset(results_all_sites, flag_change_type == "sssL")
results_O_decreases <- results_O_decreases[order(results_O_decreases$orangutan_effect_size, decreasing=FALSE),]

locations_O_decreases <- results_O_decreases[with(results_O_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_O_decreases = nrow(results_O_decreases)

# ================================================================================================================================================================================================================
# Differential subset 5 of 16:  changes in one species - macaque and human/chimpanzee/gorilla/orangutan
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_M_all <- subset(results_all_sites, flag_change_type %in% c("GGGG","LLLL"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_M_all$temp <- as.numeric(apply(results_M_all, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[40])), abs(as.numeric(x[41])), abs(as.numeric(x[42])) ) }))
results_M_all <- results_M_all[order(results_M_all$temp, decreasing=TRUE),1:75]

locations_M_all <- results_M_all[with(results_M_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_M_all = nrow(results_M_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases (macaque increases are the same as HCGO decreases; the flag is LLLL)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_M_increases <- subset(results_all_sites, flag_change_type == "LLLL")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_M_increases$temp <- as.numeric(apply(results_M_increases, 1, function(x) { max( abs(as.numeric(x[22])), abs(as.numeric(x[25])), abs(as.numeric(x[28])), abs(as.numeric(x[31])) ) }))
results_M_increases <- results_M_increases[order(results_M_increases$temp, decreasing=TRUE),1:75]

locations_M_increases <- results_M_increases[with(results_M_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_M_increases = nrow(results_M_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases (macaque decreases are the same as HCGO increases; the flag is GGGG)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_M_decreases <- subset(results_all_sites, flag_change_type == "GGGG")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_M_decreases$temp <- as.numeric(apply(results_M_decreases, 1, function(x) { max( abs(as.numeric(x[22])), abs(as.numeric(x[25])), abs(as.numeric(x[28])), abs(as.numeric(x[31])) ) }))
results_M_decreases <- results_M_decreases[order(results_M_decreases$temp, decreasing=TRUE),1:75]

locations_M_decreases <- results_M_decreases[with(results_M_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_M_decreases = nrow(results_M_decreases)

# ================================================================================================================================================================================================================
# Differential subset 6 of 16:  changes in two species - human/chimpanzee
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HC_all <- subset(results_all_sites, flag_change_type %in% c("GGss","LLss"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HC_all$temp <- as.numeric(apply(results_HC_all, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[40])) ) }))
results_HC_all <- results_HC_all[order(results_HC_all$temp, decreasing=TRUE),1:75]

locations_HC_all <- results_HC_all[with(results_HC_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_HC_all = nrow(results_HC_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HC_increases <- subset(results_all_sites, flag_change_type == "GGss")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HC_increases$temp <- as.numeric(apply(results_HC_increases, 1, function(x) { max( as.numeric(x[39]), as.numeric(x[40]) ) }))
results_HC_increases <- results_HC_increases[order(results_HC_increases$temp, decreasing=TRUE),1:75]

locations_HC_increases <- results_HC_increases[with(results_HC_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HC_increases = nrow(results_HC_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HC_decreases <- subset(results_all_sites, flag_change_type == "LLss")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HC_decreases$temp <- as.numeric(apply(results_HC_decreases, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[40])) ) }))
results_HC_decreases <- results_HC_decreases[order(results_HC_decreases$temp, decreasing=TRUE),1:75]

locations_HC_decreases <- results_HC_decreases[with(results_HC_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HC_decreases = nrow(results_HC_decreases)

# ================================================================================================================================================================================================================
# Differential subset 7 of 16:  changes in two species - human/gorilla
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HG_all <- subset(results_all_sites, flag_change_type %in% c("GsGs","LsLs"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HG_all$temp <- as.numeric(apply(results_HG_all, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[41])) ) }))
results_HG_all <- results_HG_all[order(results_HG_all$temp, decreasing=TRUE),1:75]

locations_HG_all <- results_HG_all[with(results_HG_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_HG_all = nrow(results_HG_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HG_increases <- subset(results_all_sites, flag_change_type == "GsGs")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HG_increases$temp <- as.numeric(apply(results_HG_increases, 1, function(x) { max( as.numeric(x[39]), as.numeric(x[41]) ) }))
results_HG_increases <- results_HG_increases[order(results_HG_increases$temp, decreasing=TRUE),1:75]

locations_HG_increases <- results_HG_increases[with(results_HG_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HG_increases = nrow(results_HG_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HG_decreases <- subset(results_all_sites, flag_change_type == "LsLs")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HG_decreases$temp <- as.numeric(apply(results_HG_decreases, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[41])) ) }))
results_HG_decreases <- results_HG_decreases[order(results_HG_decreases$temp, decreasing=TRUE),1:75]

locations_HG_decreases <- results_HG_decreases[with(results_HG_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HG_decreases = nrow(results_HG_decreases)

# ================================================================================================================================================================================================================
# Differential subset 8 of 16:  changes in two species - human/orangutan
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HO_all <- subset(results_all_sites, flag_change_type %in% c("GssG","LssL"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HO_all$temp <- as.numeric(apply(results_HO_all, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[42])) ) }))
results_HO_all <- results_HO_all[order(results_HO_all$temp, decreasing=TRUE),1:75]

locations_HO_all <- results_HO_all[with(results_HO_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_HO_all = nrow(results_HO_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HO_increases <- subset(results_all_sites, flag_change_type == "GssG")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HO_increases$temp <- as.numeric(apply(results_HO_increases, 1, function(x) { max( as.numeric(x[39]), as.numeric(x[42]) ) }))
results_HO_increases <- results_HO_increases[order(results_HO_increases$temp, decreasing=TRUE),1:75]

locations_HO_increases <- results_HO_increases[with(results_HO_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HO_increases = nrow(results_HO_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HO_decreases <- subset(results_all_sites, flag_change_type == "LssL")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HO_decreases$temp <- as.numeric(apply(results_HO_decreases, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[42])) ) }))
results_HO_decreases <- results_HO_decreases[order(results_HO_decreases$temp, decreasing=TRUE),1:75]

locations_HO_decreases <- results_HO_decreases[with(results_HO_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HO_decreases = nrow(results_HO_decreases)

# ================================================================================================================================================================================================================
# Differential subset 9 of 16:  changes in two species - chimpanzee/gorilla
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_CG_all <- subset(results_all_sites, flag_change_type %in% c("sGGs","sLLs"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_CG_all$temp <- as.numeric(apply(results_CG_all, 1, function(x) { max( abs(as.numeric(x[40])), abs(as.numeric(x[41])) ) }))
results_CG_all <- results_CG_all[order(results_CG_all$temp, decreasing=TRUE),1:75]

locations_CG_all <- results_CG_all[with(results_CG_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_CG_all = nrow(results_CG_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_CG_increases <- subset(results_all_sites, flag_change_type == "sGGs")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_CG_increases$temp <- as.numeric(apply(results_CG_increases, 1, function(x) { max( as.numeric(x[40]), as.numeric(x[41]) ) }))
results_CG_increases <- results_CG_increases[order(results_CG_increases$temp, decreasing=TRUE),1:75]

locations_CG_increases <- results_CG_increases[with(results_CG_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_CG_increases = nrow(results_CG_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_CG_decreases <- subset(results_all_sites, flag_change_type == "sLLs")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_CG_decreases$temp <- as.numeric(apply(results_CG_decreases, 1, function(x) { max( abs(as.numeric(x[40])), abs(as.numeric(x[41])) ) }))
results_CG_decreases <- results_CG_decreases[order(results_CG_decreases$temp, decreasing=TRUE),1:75]

locations_CG_decreases <- results_CG_decreases[with(results_CG_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_CG_decreases = nrow(results_CG_decreases)

# ================================================================================================================================================================================================================
# Differential subset 10 of 16:  changes in two species - chimpanzee/orangutan
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_CO_all <- subset(results_all_sites, flag_change_type %in% c("sGsG","sLsL"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_CO_all$temp <- as.numeric(apply(results_CO_all, 1, function(x) { max( abs(as.numeric(x[40])), abs(as.numeric(x[42])) ) }))
results_CO_all <- results_CO_all[order(results_CO_all$temp, decreasing=TRUE),1:75]

locations_CO_all <- results_CO_all[with(results_CO_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_CO_all = nrow(results_CO_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_CO_increases <- subset(results_all_sites, flag_change_type == "sGsG")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_CO_increases$temp <- as.numeric(apply(results_CO_increases, 1, function(x) { max( as.numeric(x[40]), as.numeric(x[42]) ) }))
results_CO_increases <- results_CO_increases[order(results_CO_increases$temp, decreasing=TRUE),1:75]

locations_CO_increases <- results_CO_increases[with(results_CO_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_CO_increases = nrow(results_CO_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_CO_decreases <- subset(results_all_sites, flag_change_type == "sLsL")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_CO_decreases$temp <- as.numeric(apply(results_CO_decreases, 1, function(x) { max( abs(as.numeric(x[40])), abs(as.numeric(x[42])) ) }))
results_CO_decreases <- results_CO_decreases[order(results_CO_decreases$temp, decreasing=TRUE),1:75]

locations_CO_decreases <- results_CO_decreases[with(results_CO_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_CO_decreases = nrow(results_CO_decreases)

# ================================================================================================================================================================================================================
# Differential subset 11 of 16:  changes in two species - gorilla/orangutan
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_GO_all <- subset(results_all_sites, flag_change_type %in% c("ssGG","ssLL"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_GO_all$temp <- as.numeric(apply(results_GO_all, 1, function(x) { max( abs(as.numeric(x[41])), abs(as.numeric(x[42])) ) }))
results_GO_all <- results_GO_all[order(results_GO_all$temp, decreasing=TRUE),1:75]

locations_GO_all <- results_GO_all[with(results_GO_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_GO_all = nrow(results_GO_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_GO_increases <- subset(results_all_sites, flag_change_type == "ssGG")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_GO_increases$temp <- as.numeric(apply(results_GO_increases, 1, function(x) { max( as.numeric(x[41]), as.numeric(x[42]) ) }))
results_GO_increases <- results_GO_increases[order(results_GO_increases$temp, decreasing=TRUE),1:75]

locations_GO_increases <- results_GO_increases[with(results_GO_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_GO_increases = nrow(results_GO_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_GO_decreases <- subset(results_all_sites, flag_change_type == "ssLL")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_GO_decreases$temp <- as.numeric(apply(results_GO_decreases, 1, function(x) { max( abs(as.numeric(x[41])), abs(as.numeric(x[42])) ) }))
results_GO_decreases <- results_GO_decreases[order(results_GO_decreases$temp, decreasing=TRUE),1:75]

locations_GO_decreases <- results_GO_decreases[with(results_GO_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_GO_decreases = nrow(results_GO_decreases)

# ================================================================================================================================================================================================================
# Differential subset 12 of 16:  changes in three species - human/chimpanzee/gorilla
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HCG_all <- subset(results_all_sites, flag_change_type %in% c("GGGs","LLLs"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HCG_all$temp <- as.numeric(apply(results_HCG_all, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[40])), abs(as.numeric(x[41])) ) }))
results_HCG_all <- results_HCG_all[order(results_HCG_all$temp, decreasing=TRUE),1:75]

locations_HCG_all <- results_HCG_all[with(results_HCG_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_HCG_all = nrow(results_HCG_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HCG_increases <- subset(results_all_sites, flag_change_type == "GGGs")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HCG_increases$temp <- as.numeric(apply(results_HCG_increases, 1, function(x) { max( as.numeric(x[39]), as.numeric(x[40]), as.numeric(x[41]) ) }))
results_HCG_increases <- results_HCG_increases[order(results_HCG_increases$temp, decreasing=TRUE),1:75]

locations_HCG_increases <- results_HCG_increases[with(results_HCG_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HCG_increases = nrow(results_HCG_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HCG_decreases <- subset(results_all_sites, flag_change_type == "LLLs")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HCG_decreases$temp <- as.numeric(apply(results_HCG_decreases, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[40])), abs(as.numeric(x[41])) ) }))
results_HCG_decreases <- results_HCG_decreases[order(results_HCG_decreases$temp, decreasing=TRUE),1:75]

locations_HCG_decreases <- results_HCG_decreases[with(results_HCG_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HCG_decreases = nrow(results_HCG_decreases)

# ================================================================================================================================================================================================================
# Differential subset 13 of 16:  changes in two species - human/chimpanzee/orangutan
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HCO_all <- subset(results_all_sites, flag_change_type %in% c("GGsG","LLsL"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HCO_all$temp <- as.numeric(apply(results_HCO_all, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[40])), abs(as.numeric(x[42])) ) }))
results_HCO_all <- results_HCO_all[order(results_HCO_all$temp, decreasing=TRUE),1:75]

locations_HCO_all <- results_HCO_all[with(results_HCO_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_HCO_all = nrow(results_HCO_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HCO_increases <- subset(results_all_sites, flag_change_type == "GGsG")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HCO_increases$temp <- as.numeric(apply(results_HCO_increases, 1, function(x) { max( as.numeric(x[39]), as.numeric(x[40]), as.numeric(x[42]) ) }))
results_HCO_increases <- results_HCO_increases[order(results_HCO_increases$temp, decreasing=TRUE),1:75]

locations_HCO_increases <- results_HCO_increases[with(results_HCO_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HCO_increases = nrow(results_HCO_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HCO_decreases <- subset(results_all_sites, flag_change_type == "LLsL")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HCO_decreases$temp <- as.numeric(apply(results_HCO_decreases, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[40])), abs(as.numeric(x[42])) ) }))
results_HCO_decreases <- results_HCO_decreases[order(results_HCO_decreases$temp, decreasing=TRUE),1:75]

locations_HCO_decreases <- results_HCO_decreases[with(results_HCO_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HCO_decreases = nrow(results_HCO_decreases)

# ================================================================================================================================================================================================================
# Differential subset 14 of 16:  changes in three species - human/gorilla/orangutan
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HGO_all <- subset(results_all_sites, flag_change_type %in% c("GsGG","LsLL"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HGO_all$temp <- as.numeric(apply(results_HGO_all, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[41])), abs(as.numeric(x[42])) ) }))
results_HGO_all <- results_HGO_all[order(results_HGO_all$temp, decreasing=TRUE),1:75]

locations_HGO_all <- results_HGO_all[with(results_HGO_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_HGO_all = nrow(results_HGO_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HGO_increases <- subset(results_all_sites, flag_change_type == "GsGG")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HGO_increases$temp <- as.numeric(apply(results_HGO_increases, 1, function(x) { max( as.numeric(x[39]), as.numeric(x[41]), as.numeric(x[42]) ) }))
results_HGO_increases <- results_HGO_increases[order(results_HGO_increases$temp, decreasing=TRUE),1:75]

locations_HGO_increases <- results_HGO_increases[with(results_HGO_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HGO_increases = nrow(results_HGO_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_HGO_decreases <- subset(results_all_sites, flag_change_type == "LsLL")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_HGO_decreases$temp <- as.numeric(apply(results_HGO_decreases, 1, function(x) { max( abs(as.numeric(x[39])), abs(as.numeric(x[41])), abs(as.numeric(x[42])) ) }))
results_HGO_decreases <- results_HGO_decreases[order(results_HGO_decreases$temp, decreasing=TRUE),1:75]

locations_HGO_decreases <- results_HGO_decreases[with(results_HGO_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_HGO_decreases = nrow(results_HGO_decreases)

# ================================================================================================================================================================================================================
# Differential subset 15 of 16:  changes in three species - chimpanzee/gorilla/orangutan
# ================================================================================================================================================================================================================
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All changes
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_CGO_all <- subset(results_all_sites, flag_change_type %in% c("sGGG","sLLL"))

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_CGO_all$temp <- as.numeric(apply(results_CGO_all, 1, function(x) { max( abs(as.numeric(x[40])), abs(as.numeric(x[41])), abs(as.numeric(x[42])) ) }))
results_CGO_all <- results_CGO_all[order(results_CGO_all$temp, decreasing=TRUE),1:75]

locations_CGO_all <- results_CGO_all[with(results_CGO_all, order(chrom, as.numeric(as.character(start)))),1:3]

count_CGO_all = nrow(results_CGO_all)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_CGO_increases <- subset(results_all_sites, flag_change_type == "sGGG")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_CGO_increases$temp <- as.numeric(apply(results_CGO_increases, 1, function(x) { max( as.numeric(x[40]), as.numeric(x[41]), as.numeric(x[42]) ) }))
results_CGO_increases <- results_CGO_increases[order(results_CGO_increases$temp, decreasing=TRUE),1:75]

locations_CGO_increases <- results_CGO_increases[with(results_CGO_increases, order(chrom, as.numeric(as.character(start)))),1:3]

count_CGO_increases = nrow(results_CGO_increases)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_CGO_decreases <- subset(results_all_sites, flag_change_type == "sLLL")

# Sort by the largest effect size (based on absolute value)
#	Find the largest effect size for each site (based on absolute value) and store in a temporary column
#	Order by that temporary column and then remove
results_CGO_decreases$temp <- as.numeric(apply(results_CGO_decreases, 1, function(x) { max( abs(as.numeric(x[40])), abs(as.numeric(x[41])), abs(as.numeric(x[42])) ) }))
results_CGO_decreases <- results_CGO_decreases[order(results_CGO_decreases$temp, decreasing=TRUE),1:75]

locations_CGO_decreases <- results_CGO_decreases[with(results_CGO_decreases, order(chrom, as.numeric(as.character(start)))),1:3]

count_CGO_decreases = nrow(results_CGO_decreases)

# ================================================================================================================================================================================================================
# Differential subset 16 of 16: other changes
# ================================================================================================================================================================================================================
results_OOOO_all <- subset(results_all_sites, flag_change_type %in% c("OOOO"))
results_OOOO_all <- results_OOOO_all[with(results_OOOO_all, order(chrom, as.numeric(as.character(start)))),]

locations_OOOO_all <- results_OOOO_all[,1:3]

count_OOOO_all = nrow(results_OOOO_all)

# ================================================================================================================================================================================================================
# Finish up
# ================================================================================================================================================================================================================
cat("done.\n\n")

# ################################################################################################################################################################################################################
# END SECTION 7 OF 10: CALCULATE TOTALS AND CREATE TABLES CONTAINING SUBSETS OF THE RESULTS
# ################################################################################################################################################################################################################


# ################################################################################################################################################################################################################
# SECTION 8 OF 10: WRITE RESULTS TO FILES
#	67 files
# ================================================================================================================================================================================================================
#	All results  (1 file)
#		file 1: txt file containing all of the information			glm_analysis.all_sites.all_information.txt
#	Not sigificantly differential sites (2 files)
#		file 1: txt file containing all of the information			glm_analysis.not_significant_sites.all_information.txt
#		file 2: bed file containing location information only			glm_analysis.not_significant_sites.locations.bed
#	All differential sites (2 files)
#		file 1: txt file containing all of the information			glm_analysis.differential_sites.all_information.txt
#		file 2: bed file containing location information only			glm_analysis.differential_sites.locations.bed
#	Differential subset 1 of 16: changes in one species - human (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.Gsss.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.Gsss.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.Lsss.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.Lsss.locations.bed
#	Differential subset 2 of 16: changes in one species - chimpanzee (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.sGss.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.sGss.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.sLss.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.sLss.locations.bed
#	Differential subset 3 of 16: changes in one species - gorilla (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.ssGs.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.ssGs.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.ssLs.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.ssLs.locations.bed
#	Differential subset 4 of 16: changes in one species - orangutan (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.sssG.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.sssG.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.sssL.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.sssL.locations.bed
#	Differential subset 5 of 16: changes in one species - macaque (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.LLLL.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.LLLL.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.GGGG.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.GGGG.locations.bed
#	Differential subset 6 of 16: changes in two species - human/chimpanzee (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.GGss.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.GGss.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.LLss.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.LLss.locations.bed
#	Differential subset 7 of 16: changes in two species - human/gorilla (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.GsGs.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.GsGs.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.LsLs.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.LsLs.locations.bed
#	Differential subset 8 of 16: changes in two species - human/orangutan (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.GssG.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.GssG.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.LssL.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.LssL.locations.bed
#	Differential subset 9 of 16: changes in two species - chimpanzee/gorilla (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.sGGs.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.sGGs.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.sLLs.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.sLLs.locations.bed
#	Differential subset 10 of 16: changes in two species - chimpanzee/orangutan (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.sGsG.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.sGsG.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.sLsL.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.sLsL.locations.bed
#	Differential subset 11 of 16: changes in two species - gorilla/orangutan (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.ssGG.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.ssGG.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.ssLL.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.ssLL.locations.bed
#	Differential subset 12 of 16: changes in three species - human/chimpanzee/gorilla (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.GGGs.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.GGGs.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.LLLs.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.LLLs.locations.bed
#	Differential subset 13 of 16: changes in three species - human/chimpanzee/orangutan (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.GGsG.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.GGsG.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.LLsL.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.LLsL.locations.bed
#	Differential subset 14 of 16: changes in three species - human/gorilla/orangutan (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.GsGG.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.GsGG.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.LsLL.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.LsLL.locations.bed
#	Differential subset 15 of 16: changes in three species - chimpanzee/gorilla/orangutan (4 files)
#		file 1: increases: txt file containing all of the information		glm_analysis.sGGG.all_information.txt
#		file 2: increases: bed file containing location information only	glm_analysis.sGGG.locations.bed
#		file 3: decreases: txt file containing all of the information		glm_analysis.sLLL.all_information.txt
#		file 4: decreases: bed file containing location information only	glm_analysis.sLLL.locations.bed
#	Differential subset 16 of 16: other changes (2 files)				
#		file 1: txt file containing all of the information			glm_analysis.OOOO.all_information.txt
#		file 2: bed file containing location information only			glm_analysis.OOOO.locations.be
# ################################################################################################################################################################################################################

cat("Writing results to files ... \n")

# ================================================================================================================================================================================================================
# All results
# ================================================================================================================================================================================================================
write.table(results_all_sites, "glm_analysis.all_sites.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

# ================================================================================================================================================================================================================
# Not significantly differential sites
# ================================================================================================================================================================================================================
write.table(results_not_significant_sites,   "glm_analysis.not_significant_sites.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_not_significant_sites, "glm_analysis.not_significant_sites.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential sites
# ================================================================================================================================================================================================================
write.table(results_differential_sites,   "glm_analysis.differential_sites.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_differential_sites, "glm_analysis.differential_sites.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 1 of 16: changes in one species - human
# ================================================================================================================================================================================================================
# Increases
write.table(results_H_increases,   "glm_analysis.Gsss.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_H_increases, "glm_analysis.Gsss.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_H_decreases,   "glm_analysis.Lsss.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_H_decreases, "glm_analysis.Lsss.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 2 of 16: changes in one species - chimpanzee
# ================================================================================================================================================================================================================
# Increases
write.table(results_C_increases,   "glm_analysis.sGss.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_C_increases, "glm_analysis.sGss.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_C_decreases,   "glm_analysis.sLss.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_C_decreases, "glm_analysis.sLss.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 3 of 16: changes in one species - gorilla
# ================================================================================================================================================================================================================
# Increases
write.table(results_G_increases,   "glm_analysis.ssGs.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_G_increases, "glm_analysis.ssGs.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_G_decreases,   "glm_analysis.ssLs.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_G_decreases, "glm_analysis.ssLs.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 4 of 16: changes in one species - orangutan
# ================================================================================================================================================================================================================
# Increases
write.table(results_O_increases,   "glm_analysis.sssG.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_O_increases, "glm_analysis.sssG.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_O_decreases,   "glm_analysis.sssL.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_O_decreases, "glm_analysis.sssL.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 5 of 16: changes in one species - macaque (also includes changes in human/chimpanzee/gorilla/orangutan)
# ================================================================================================================================================================================================================
# Macaque Increases / HCGO Decreases
write.table(results_M_increases,   "glm_analysis.LLLL.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_M_increases, "glm_analysis.LLLL.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Macaque Decreases / HCGO Increases
write.table(results_M_decreases,   "glm_analysis.GGGG.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_M_decreases, "glm_analysis.GGGG.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 6 of 16: changes in two species - human/chimpanzee
# ================================================================================================================================================================================================================
# Increases
write.table(results_HC_increases,   "glm_analysis.GGss.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HC_increases, "glm_analysis.GGss.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_HC_decreases,   "glm_analysis.LLss.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HC_decreases, "glm_analysis.LLss.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 7 of 16: changes in two species - human/gorilla
# ================================================================================================================================================================================================================
# Increases
write.table(results_HG_increases,   "glm_analysis.GsGs.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HG_increases, "glm_analysis.GsGs.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_HG_decreases,   "glm_analysis.LsLs.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HG_decreases, "glm_analysis.LsLs.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 8 of 16: changes in two species - human/orangutan
# ================================================================================================================================================================================================================
# Increases
write.table(results_HO_increases,   "glm_analysis.GssG.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HO_increases, "glm_analysis.GssG.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_HO_decreases,   "glm_analysis.LssL.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HO_decreases, "glm_analysis.LssL.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 9 of 16: changes in two species - chimpanzee/gorilla
# ================================================================================================================================================================================================================
# Increases
write.table(results_CG_increases,   "glm_analysis.sGGs.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_CG_increases, "glm_analysis.sGGs.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_CG_decreases,   "glm_analysis.sLLs.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_CG_decreases, "glm_analysis.sLLs.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 10 of 16: changes in two species - chimpanzee/orangutan
# ================================================================================================================================================================================================================
# Increases
write.table(results_CO_increases,   "glm_analysis.sGsG.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_CO_increases, "glm_analysis.sGsG.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_CO_decreases,   "glm_analysis.sLsL.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_CO_decreases, "glm_analysis.sLsL.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 11 of 16: changes in two species - gorilla/orangutan
# ================================================================================================================================================================================================================
# Increases
write.table(results_GO_increases,   "glm_analysis.ssGG.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_GO_increases, "glm_analysis.ssGG.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_GO_decreases,   "glm_analysis.ssLL.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_GO_decreases, "glm_analysis.ssLL.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 12 of 16: changes in three species - human/chimpanzee/gorilla
# ================================================================================================================================================================================================================
# Increases
write.table(results_HCG_increases,   "glm_analysis.GGGs.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HCG_increases, "glm_analysis.GGGs.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_HCG_decreases,   "glm_analysis.LLLs.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HCG_decreases, "glm_analysis.LLLs.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 13 of 16: changes in three species - human/chimpanzee/orangutan
# ================================================================================================================================================================================================================
# Increases
write.table(results_HCO_increases,   "glm_analysis.GGsG.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HCO_increases, "glm_analysis.GGsG.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_HCO_decreases,   "glm_analysis.LLsL.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HCO_decreases, "glm_analysis.LLsL.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 14 of 16: changes in three species - human/gorilla/orangutan
# ================================================================================================================================================================================================================
# Increases
write.table(results_HGO_increases,   "glm_analysis.GsGG.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HGO_increases, "glm_analysis.GsGG.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_HGO_decreases,   "glm_analysis.LsLL.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_HGO_decreases, "glm_analysis.LsLL.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 15 of 16: changes in three species - chimpanzee/gorilla/orangutan
# ================================================================================================================================================================================================================
# Increases
write.table(results_CGO_increases,   "glm_analysis.sGGG.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_CGO_increases, "glm_analysis.sGGG.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Decreases
write.table(results_CGO_decreases,   "glm_analysis.sLLL.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_CGO_decreases, "glm_analysis.sLLL.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Differential subset 16 of 16 - other changes
# ================================================================================================================================================================================================================
write.table(results_OOOO_all,   "glm_analysis.OOOO.all_information.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(locations_OOOO_all, "glm_analysis.OOOO.locations.bed",       sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# ================================================================================================================================================================================================================
# Done
# ================================================================================================================================================================================================================
cat("done.\n\n")

# ################################################################################################################################################################################################################
# END SECTION 8 OF 10: WRITE RESULTS TO FILES
# ################################################################################################################################################################################################################


# ################################################################################################################################################################################################################
# SECTION 9 OF 10: CALCULATE CHECKSUM AND PRINT STATS
# ################################################################################################################################################################################################################

# ================================================================================================================================================================================================================
# Calculate checksums
# ================================================================================================================================================================================================================
# Checksum1 = "all sites" - "not significantly differential sites" - "differential sites"
checksum1 = count_all_sites - count_not_significant_sites - count_differential_sites

# Checksum2 = "all differential sites" - "each subtype"
checksum2 = count_differential_sites - count_H_all - count_C_all - count_G_all - count_O_all - count_M_all - count_HC_all - count_HG_all - count_HO_all - count_CG_all - count_CO_all - count_GO_all - count_HCG_all - count_HCO_all - count_HGO_all - count_CGO_all - count_OOOO_all

cat("=================================================================================================\n")
cat("Checksums (should be 0)\n")
cat("=================================================================================================\n")
cat("Checksum (\"all sites\" - \"not significantly differential sites\" - \"differential sites\"):", checksum1, "\n")
cat("Checksum (\"differential sites\" - each subtype):", checksum2, "\n\n")

cat("=================================================================================================\n")
cat("Counts\n")
cat("=================================================================================================\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("All sites\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("     Not significantly differential sites\t", count_not_significant_sites, "\n")
cat("     Differential sites\t\t\t\t", count_differential_sites, "\n")
cat("     Total\t\t\t\t\t", count_all_sites, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Single changes affecting one species\n")
cat("-------------------------------------------------------------------------------------------------\n")

cat("Human\n")
cat("     Increases\t", count_H_increases, "\n")
cat("     Decreases\t", count_H_decreases, "\n")
cat("     Total\t", count_H_all, "\n\n")

cat("Chimpanzee\n")
cat("     Increases\t", count_C_increases, "\n")
cat("     Decreases\t", count_C_decreases, "\n")
cat("     Total\t", count_C_all, "\n\n")

cat("Gorilla\n")
cat("     Increases\t", count_G_increases, "\n")
cat("     Decreases\t", count_G_decreases, "\n")
cat("     Total\t", count_G_all, "\n\n")

cat("Orangutan\n")
cat("     Increases\t", count_O_increases, "\n")
cat("     Decreases\t", count_O_decreases, "\n")
cat("     Total\t", count_O_all, "\n\n")

cat("Macaque\n")
cat("     Increases\t", count_M_increases, "\n")
cat("     Decreases\t", count_M_decreases, "\n")
cat("     Total\t", count_M_all, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Single changes affecting multiple species\n")
cat("-------------------------------------------------------------------------------------------------\n")

cat("Human/chimpanzee\n")
cat("     Increases\t", count_HC_increases, "\n")
cat("     Decreases\t", count_HC_decreases, "\n")
cat("     Total\t", count_HC_all, "\n\n")

cat("Human/chimpanzee/gorilla\n")
cat("     Increases\t", count_HCG_increases, "\n")
cat("     Decreases\t", count_HCG_decreases, "\n")
cat("     Total\t", count_HCG_all, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Multiple changes affecting multiple species\n")
cat("-------------------------------------------------------------------------------------------------\n")

cat("Human/gorilla\n")
cat("     Increases\t", count_HG_increases, "\n")
cat("     Decreases\t", count_HG_decreases, "\n")
cat("     Total\t", count_HG_all, "\n\n")

cat("Human/orangutan\n")
cat("     Increases\t", count_HO_increases, "\n")
cat("     Decreases\t", count_HO_decreases, "\n")
cat("     Total\t", count_HO_all, "\n\n")

cat("Chimpanzee/gorilla\n")
cat("     Increases\t", count_CG_increases, "\n")
cat("     Decreases\t", count_CG_decreases, "\n")
cat("     Total\t", count_CG_all, "\n\n")

cat("Chimpanzee/orangutan\n")
cat("     Increases\t", count_CO_increases, "\n")
cat("     Decreases\t", count_CO_decreases, "\n")
cat("     Total\t", count_CO_all, "\n\n")

cat("Gorilla/orangutan\n")
cat("     Increases\t", count_GO_increases, "\n")
cat("     Decreases\t", count_GO_decreases, "\n")
cat("     Total\t", count_GO_all, "\n\n")

cat("Human/chimpanzee/orangutan\n")
cat("     Increases\t", count_HCO_increases, "\n")
cat("     Decreases\t", count_HCO_decreases, "\n")
cat("     Total\t", count_HCO_all, "\n\n")

cat("Human/gorilla/orangutan\n")
cat("     Increases\t", count_HGO_increases, "\n")
cat("     Decreases\t", count_HGO_decreases, "\n")
cat("     Total\t", count_HGO_all, "\n\n")

cat("Chimpanzee/gorilla/orangutan\n")
cat("     Increases\t", count_CGO_increases, "\n")
cat("     Decreases\t", count_CGO_decreases, "\n")
cat("     Total\t", count_CGO_all, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Other changes\t", count_OOOO_all, "\n")
cat("-------------------------------------------------------------------------------------------------\n\n")

# ################################################################################################################################################################################################################
# END SECTION 9 OF 10: CALCULATE CHECKSUM AND PRINT STATS
# ################################################################################################################################################################################################################


# ################################################################################################################################################################################################################
# SECTION 10 OF 10: FINISH UP
# ################################################################################################################################################################################################################

cat("-------------------------------------------------------------------------------------------------\n")
cat("END:", date(), "\n")
cat("-------------------------------------------------------------------------------------------------\n")

sink()
