# ################################################################################################################################################################################################################
# GO.run_edgeR.R
#
# CC-BY Lee Edsall
#	email: edsall57@gmail.com
#	Twitter: @LeeEdsall
#	GitHub: https://github.com/ledsall
#
# This script was used to analyze data for Edsall et al. 2019 which compared DNase-seq data from 5 primates (human, chimpanzee, gorilla, orangutan, macaque)
# ################################################################################################################################################################################################################


# ################################################################################################################################################################################################################
# SECTION 1: SETUP
#       1. load libraries
#       2. set global options
#       3. set global variables
#       4. open connection to log file
# ################################################################################################################################################################################################################

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# edgeR		Empirical Analysis of Digital Gene Expression Data in R
#		Citation:  Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139-140.
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(edgeR)
library(VennDiagram)
library(gridExtra)

options(scipen=0)
pvalue_threshold = .01

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Open the log file; write date and options
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sink(file="LOG.run_edgeR.R", append=FALSE, type="output", split=TRUE)
cat("-------------------------------------------------------------------------------------------------\n")
cat("BEGIN:", date(), "\n")
cat("-------------------------------------------------------------------------------------------------\n\n")

cat("Using the following settings\n")
cat("\tLog file: LOG.run_edgeR.R\n")
cat("\tp-value threshold:",  pvalue_threshold, "\n\n")

# ################################################################################################################################################################################################################
# END SECTION 1: SETUP
# ################################################################################################################################################################################################################

# ################################################################################################################################################################################################################
# SECTION 2: READ THE GLM RESULTS FILE
# ################################################################################################################################################################################################################

cat("Reading GLM results file (glm_output_and_analysis_results_without_pairwise_analysis_results.txt)... \n")

glm_results <- as.data.frame(read.table("glm_output_and_analysis_results_without_pairwise_analysis_results.txt", sep="\t", header=TRUE))

rownames(glm_results) <- paste(glm_results$chrom, glm_results$start, glm_results$end, glm_results$PLoS_overlap, glm_results$flag_change_type, sep=":")

count_all_sites = nrow(glm_results)

cat("done.\n\n")

# ################################################################################################################################################################################################################
# END SECTION 2: READ THE GLM RESULTS FILE
# ################################################################################################################################################################################################################

# ################################################################################################################################################################################################################
# SECTION 3: RUN edgeR AND CREATE RESULTS TABLES
# ################################################################################################################################################################################################################

cat("Setting up to run edgeR ... \n")

# ================================================================================================================================================================================================================
# Set up the read counts matrix
#	Extract just the columns with the scores
#	Convert from a data.frame to a matrix
# ================================================================================================================================================================================================================

scores <- as.matrix(glm_results[,c(44:58)])

# ================================================================================================================================================================================================================
# Set up the design matrix
#	The species order is human, chimpanzee, gorilla, orangutan, macaque
#	The factors automatically get re-sorted, so the scores table does not need to be in order
# ================================================================================================================================================================================================================

species <- factor(c(2,2,2,3,3,3,4,4,4,5,5,5,1,1,1))
design <- model.matrix(~species)

# ================================================================================================================================================================================================================
# Prepare to run edgeR
# 	1. Calculate the normalization factor
#	2. Estimate the dispersion factors
# ================================================================================================================================================================================================================

y <- DGEList(counts=scores,group=species)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)

# ================================================================================================================================================================================================================
# Run edgeR
#	1. Run the glmQLFit command to generate the model that can be used to run the Quasi-likelihood F test
#	2. Run the human vs. macaque comparison
#	3. Run the chimpanzee vs. macaque comparison
#	4. Run the gorilla vs. macaque comparison
#	5. Run the orangutan vs. macaque comparison
# ================================================================================================================================================================================================================

cat("Running edgeR ... \n")

fit <- glmQLFit(y,design)

qlf.h <- glmQLFTest(fit, coef=2)
qlf.c <- glmQLFTest(fit, coef=3)
qlf.g <- glmQLFTest(fit, coef=4)
qlf.o <- glmQLFTest(fit, coef=5)

cat("done.\n\n")

# ================================================================================================================================================================================================================
# Create a master table of results
#	1. Extract the results from edgeR
#	2. Perform a BH correction
#	3. Create and update an "edgeR_flag_differential" (same idea as the GLM "flag_gateway_test")
#	4. Create and update an "edgeR_flag_change_type" (same idea as the GLM "flag_change_type" flag)
#	5. Add the GLM results table
# ================================================================================================================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract the edgeR results
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("Extracting the edgeR results ... \n")

h <- as.data.frame(qlf.h$table)
c <- as.data.frame(qlf.c$table)
g <- as.data.frame(qlf.g$table)
o <- as.data.frame(qlf.o$table)


colnames(h) <- c("h_log2FC", "h_log2CPM", "h_F", "h_uncorrected_pvalue")
colnames(c) <- c("c_log2FC", "c_log2CPM", "c_F", "c_uncorrected_pvalue")
colnames(g) <- c("g_log2FC", "g_log2CPM", "g_F", "g_uncorrected_pvalue")
colnames(o) <- c("o_log2FC", "o_log2CPM", "o_F", "o_uncorrected_pvalue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Add columns with the ln fold change to match the beta value scale
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("Calculating ln fold change values ... \n")

h$h_lnFC <- log(2^h$h_log2FC)
c$c_lnFC <- log(2^c$c_log2FC)
g$g_lnFC <- log(2^g$g_log2FC)
o$o_lnFC <- log(2^o$o_log2FC)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Perform a BH correction
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("Performing a BH correction ... \n")

uncorrected_pvalues <- as.vector(cbind(h[,4], c[,4], g[,4], o[,4]))
corrected_pvalues <- p.adjust(uncorrected_pvalues, "BH")

h$h_pvalue <- corrected_pvalues[1:89744]
c$c_pvalue <- corrected_pvalues[89745:179488]
g$g_pvalue <- corrected_pvalues[179489:269232]
o$o_pvalue <- corrected_pvalues[269233:358976]

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Combine the edgeR results tables and the GLM results table
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("Creating a master table of results ... \n")

edgeR_results <- as.data.frame(cbind(glm_results, h, c, g, o))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Add columns with the flags
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
edgeR_results$edgeR_flag_differential <- "ENS"
edgeR_results$edgeR_flag_change_type <- "ssss"


# ================================================================================================================================================================================================================
# Identify the type of change and update the edgeR flags
# ================================================================================================================================================================================================================

cat("Identifying differential sites and determining type of change ... \n")
cat("\tAnalyzing site 1 of", count_all_sites, "\n")

for (lcv in 1:count_all_sites)
{
	# Print periodic status messages
	if ( (lcv %% 1000) == 0)
	{
		cat("\tAnalyzing site", lcv, "of", count_all_sites, "\n")
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Identify differential sites
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	if (edgeR_results$h_pvalue[lcv] < .01 | edgeR_results$c_pvalue[lcv] < .01 | edgeR_results$g_pvalue[lcv] < .01 | edgeR_results$o_pvalue[lcv] < .01)  
	{
		edgeR_results$edgeR_flag_differential[lcv] <- "ED"
	}


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 1 of 15: Change in human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] < pvalue_threshold & edgeR_results$c_pvalue[lcv] > pvalue_threshold & edgeR_results$g_pvalue[lcv] > pvalue_threshold & edgeR_results$o_pvalue[lcv] > pvalue_threshold)
	{
		if (edgeR_results$h_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "Gsss"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "Lsss"	
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 2 of 15: Change in chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] > pvalue_threshold & edgeR_results$c_pvalue[lcv] < pvalue_threshold & edgeR_results$g_pvalue[lcv] > pvalue_threshold & edgeR_results$o_pvalue[lcv] > pvalue_threshold)
	{
		if (edgeR_results$c_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sGss"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "sLss"	
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 3 of 15: Change in gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] > pvalue_threshold & edgeR_results$c_pvalue[lcv] > pvalue_threshold & edgeR_results$g_pvalue[lcv] < pvalue_threshold & edgeR_results$o_pvalue[lcv] > pvalue_threshold)
	{
		if (edgeR_results$g_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "ssGs"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "ssLs"	
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 4 of 15: Change in orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] > pvalue_threshold & edgeR_results$c_pvalue[lcv] > pvalue_threshold & edgeR_results$g_pvalue[lcv] > pvalue_threshold & edgeR_results$o_pvalue[lcv] < pvalue_threshold)
	{
		if (edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sssG"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "sssL"	
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 5 of 15: Changes in HCGO
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] < pvalue_threshold & edgeR_results$c_pvalue[lcv] < pvalue_threshold & edgeR_results$g_pvalue[lcv] < pvalue_threshold & edgeR_results$o_pvalue[lcv] < pvalue_threshold)
	{
		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GGGG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GGGL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GGLG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GLGG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LGGG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GGLL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GLGL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GLLG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LGGL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LGLG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LLGG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GLLL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LGLL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LLGL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LLLG"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "LLLL"
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 6 of 15: Changes in human/chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] < pvalue_threshold & edgeR_results$c_pvalue[lcv] < pvalue_threshold & edgeR_results$g_pvalue[lcv] > pvalue_threshold & edgeR_results$o_pvalue[lcv] > pvalue_threshold)
	{
		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GGss"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GLss"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LGss"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "LLss"
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 7 of 15: Changes in human/gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] < pvalue_threshold & edgeR_results$c_pvalue[lcv] > pvalue_threshold & edgeR_results$g_pvalue[lcv] < pvalue_threshold & edgeR_results$o_pvalue[lcv] > pvalue_threshold)
	{
		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GsGs"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GsLs"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LsGs"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "LsLs"
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 8 of 15: Changes in human/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] < pvalue_threshold & edgeR_results$c_pvalue[lcv] > pvalue_threshold & edgeR_results$g_pvalue[lcv] > pvalue_threshold & edgeR_results$o_pvalue[lcv] < pvalue_threshold)
	{
		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GssG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GssL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LssG"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "LssL"
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 9 of 15: Changes in chimpanzee/gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] > pvalue_threshold & edgeR_results$c_pvalue[lcv] < pvalue_threshold & edgeR_results$g_pvalue[lcv] < pvalue_threshold & edgeR_results$o_pvalue[lcv] > pvalue_threshold)
	{
		if (edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sGGs"
			next
		}

		if (edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sGLs"
			next
		}

		if (edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sLGs"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "sLLs"
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 10 of 15: Changes in chimpanzee/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] > pvalue_threshold & edgeR_results$c_pvalue[lcv] < pvalue_threshold & edgeR_results$g_pvalue[lcv] > pvalue_threshold & edgeR_results$o_pvalue[lcv] < pvalue_threshold)
	{
		if (edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sGsG"
			next
		}

		if (edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sGsL"
			next
		}

		if (edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sLsG"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "sLsL"
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 11 of 15: Changes in gorilla/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] > pvalue_threshold & edgeR_results$c_pvalue[lcv] > pvalue_threshold & edgeR_results$g_pvalue[lcv] < pvalue_threshold & edgeR_results$o_pvalue[lcv] < pvalue_threshold)
	{
		if (edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "ssGG"
			next
		}

		if (edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "ssGL"
			next
		}

		if (edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "ssLG"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "ssLL"
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 12 of 15: Changes in human/chimpanzee/gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] < pvalue_threshold & edgeR_results$c_pvalue[lcv] < pvalue_threshold & edgeR_results$g_pvalue[lcv] < pvalue_threshold & edgeR_results$o_pvalue[lcv] > pvalue_threshold)
	{
		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GGGs"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GGLs"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LGGs"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GLLs"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LGLs"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LLGs"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "LLLs"
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 13 of 15: Changes in human/chimpanzee/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] < pvalue_threshold & edgeR_results$c_pvalue[lcv] < pvalue_threshold & edgeR_results$g_pvalue[lcv] > pvalue_threshold & edgeR_results$o_pvalue[lcv] < pvalue_threshold)
	{
		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GGsG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GGsL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LGsG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GLsL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LGsL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LLsG"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "LLsL"
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 14 of 15: Changes in human/gorilla/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] < pvalue_threshold & edgeR_results$c_pvalue[lcv] > pvalue_threshold & edgeR_results$g_pvalue[lcv] < pvalue_threshold & edgeR_results$o_pvalue[lcv] < pvalue_threshold)
	{
		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GsGG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GsGL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LsGG"
			next
		}

		if (edgeR_results$h_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "GsLL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LsGL"
			next
		}

		if (edgeR_results$h_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "LsLG"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "LsLL"
		next
	}

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test 15 of 15: Changes in chimpanzee/gorilla/orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (edgeR_results$h_pvalue[lcv] > pvalue_threshold & edgeR_results$c_pvalue[lcv] < pvalue_threshold & edgeR_results$g_pvalue[lcv] < pvalue_threshold & edgeR_results$o_pvalue[lcv] < pvalue_threshold)
	{
		if (edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sGGG"
			next
		}

		if (edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sGGL"
			next
		}

		if (edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sLGG"
			next
		}

		if (edgeR_results$c_lnFC[lcv] > 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sGLL"
			next
		}

		if (edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] > 0 & edgeR_results$o_lnFC[lcv] < 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sLGL"
			next
		}

		if (edgeR_results$c_lnFC[lcv] < 0 & edgeR_results$g_lnFC[lcv] < 0 & edgeR_results$o_lnFC[lcv] > 0)
		{
			edgeR_results$edgeR_flag_change_type[lcv] <- "sLLG"
			next
		}

		edgeR_results$edgeR_flag_change_type[lcv] <- "sLLL"
		next
	}

}  # END: for (lcv in 1:count_all_sites)

cat("done.\n\n")

# ################################################################################################################################################################################################################
# END SECTION 3: RUN edgeR AND CREATE RESULTS TABLE
# ################################################################################################################################################################################################################

# ################################################################################################################################################################################################################
# SECTION 4: WRITE RESULTS TO FILE
# ################################################################################################################################################################################################################

cat("Writing results to files ... \n")

write.table(edgeR_results, "edgeR_analysis.all_results.all_information.txt", sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)

# ################################################################################################################################################################################################################
# END SECTION 4: WRITE RESULTS TO FILE
# ################################################################################################################################################################################################################

# ################################################################################################################################################################################################################
# SECTION 5: CREATE SUBSETS
#	For each species, create 6 subsets
#		1. Increase called by both methods
#		2. Increase called by edgeR only
#		3. Increase called by the GLM only
#		4. Decrease called by both methods
#		5. Decrease called by edgeR only
#		6. Decrease called by the GLM only
#	Note: The subsets include every category where the species changed; not just the species-specific one
# ################################################################################################################################################################################################################

# ================================================================================================================================================================================================================
# Non differential sites
# ================================================================================================================================================================================================================
non_differential.both  <- subset(edgeR_results, edgeR_flag_differential == "ENS" & flag_gateway_test == "GNS")
non_differential.edgeR <- subset(edgeR_results, edgeR_flag_differential == "ENS" & flag_gateway_test == "GD")
non_differential.glm   <- subset(edgeR_results, edgeR_flag_differential == "ED"  & flag_gateway_test == "GNS")

# ================================================================================================================================================================================================================
# Differential sites
# ================================================================================================================================================================================================================
differential.both  <- subset(edgeR_results, edgeR_flag_differential == "ED"  & flag_gateway_test == "GD")
differential.edgeR <- subset(edgeR_results, edgeR_flag_differential == "ED"  & flag_gateway_test == "GNS")
differential.glm   <- subset(edgeR_results, edgeR_flag_differential == "ENS" & flag_gateway_test == "GD")

# ================================================================================================================================================================================================================
# Human increases
# ================================================================================================================================================================================================================
h_inc.both  <- subset(edgeR_results, h_pvalue < .01 & h_lnFC > 0 &  flag_change_type %in% list("Gsss","GGss","GsGs","GssG","GGGs","GGsG","GsGG","GGGG"))
h_inc.edgeR <- subset(edgeR_results, h_pvalue < .01 & h_lnFC > 0 & !flag_change_type %in% list("Gsss","GGss","GsGs","GssG","GGGs","GGsG","GsGG","GGGG"))
h_inc.glm   <- subset(edgeR_results, h_pvalue > .01 & h_lnFC > 0 &  flag_change_type %in% list("Gsss","GGss","GsGs","GssG","GGGs","GGsG","GsGG","GGGG"))

# ================================================================================================================================================================================================================
# Human decreases
# ================================================================================================================================================================================================================
h_dec.both  <- subset(edgeR_results, h_pvalue < .01 & h_lnFC < 0 &  flag_change_type %in% list("Lsss","LLss","LsLs","LssL","LLLs","LLsL","LsLL","LLLL"))
h_dec.edgeR <- subset(edgeR_results, h_pvalue < .01 & h_lnFC < 0 & !flag_change_type %in% list("Lsss","LLss","LsLs","LssL","LLLs","LLsL","LsLL","LLLL"))
h_dec.glm   <- subset(edgeR_results, h_pvalue > .01 & h_lnFC < 0 &  flag_change_type %in% list("Lsss","LLss","LsLs","LssL","LLLs","LLsL","LsLL","LLLL"))

# ================================================================================================================================================================================================================
# Chimpanzee increases
# ================================================================================================================================================================================================================
c_inc.both  <- subset(edgeR_results, c_pvalue < .01 & c_lnFC > 0 &  flag_change_type %in% list("sGss","GGss","sGGs","sGsG","GGGs","GGsG","sGGG","GGGG"))
c_inc.edgeR <- subset(edgeR_results, c_pvalue < .01 & c_lnFC > 0 & !flag_change_type %in% list("sGss","GGss","sGGs","sGsG","GGGs","GGsG","sGGG","GGGG"))
c_inc.glm   <- subset(edgeR_results, c_pvalue > .01 & c_lnFC > 0 &  flag_change_type %in% list("sGss","GGss","sGGs","sGsG","GGGs","GGsG","sGGG","GGGG"))

# ================================================================================================================================================================================================================
# Chimpanzee decreases
# ================================================================================================================================================================================================================
c_dec.both  <- subset(edgeR_results, c_pvalue < .01 & c_lnFC < 0 &  flag_change_type %in% list("sLss","LLss","sLLs","sLsL","LLLs","LLsL","sLLL","LLLL"))
c_dec.edgeR <- subset(edgeR_results, c_pvalue < .01 & c_lnFC < 0 & !flag_change_type %in% list("sLss","LLss","sLLs","sLsL","LLLs","LLsL","sLLL","LLLL"))
c_dec.glm   <- subset(edgeR_results, c_pvalue > .01 & c_lnFC < 0 &  flag_change_type %in% list("sLss","LLss","sLLs","sLsL","LLLs","LLsL","sLLL","LLLL"))

# ================================================================================================================================================================================================================
# Gorilla increases
# ================================================================================================================================================================================================================
g_inc.both  <- subset(edgeR_results, g_pvalue < .01 & g_lnFC > 0 &  flag_change_type %in% list("ssGs","GsGs","sGGs","ssGG","GGGs","GsGG","sGGG","GGGG"))
g_inc.edgeR <- subset(edgeR_results, g_pvalue < .01 & g_lnFC > 0 & !flag_change_type %in% list("ssGs","GsGs","sGGs","ssGG","GGGs","GsGG","sGGG","GGGG"))
g_inc.glm   <- subset(edgeR_results, g_pvalue > .01 & g_lnFC > 0 &  flag_change_type %in% list("ssGs","GsGs","sGGs","ssGG","GGGs","GsGG","sGGG","GGGG"))

# ================================================================================================================================================================================================================
# Gorilla decreases
# ================================================================================================================================================================================================================
g_dec.both  <- subset(edgeR_results, g_pvalue < .01 & g_lnFC < 0 &  flag_change_type %in% list("ssLs","LsLs","sLLs","ssLL","LLLs","LsLL","sLLL","LLLL"))
g_dec.edgeR <- subset(edgeR_results, g_pvalue < .01 & g_lnFC < 0 & !flag_change_type %in% list("ssLs","LsLs","sLLs","ssLL","LLLs","LsLL","sLLL","LLLL"))
g_dec.glm   <- subset(edgeR_results, g_pvalue > .01 & g_lnFC < 0 &  flag_change_type %in% list("ssLs","LsLs","sLLs","ssLL","LLLs","LsLL","sLLL","LLLL"))

# ================================================================================================================================================================================================================
# Orangutan increases
# ================================================================================================================================================================================================================
o_inc.both  <- subset(edgeR_results, o_pvalue < .01 & o_lnFC > 0 &  flag_change_type %in% list("sssG","GssG","sGsG","ssGG","GGsG","GsGG","sGGG","GGGG"))
o_inc.edgeR <- subset(edgeR_results, o_pvalue < .01 & o_lnFC > 0 & !flag_change_type %in% list("sssG","GssG","sGsG","ssGG","GGsG","GsGG","sGGG","GGGG"))
o_inc.glm   <- subset(edgeR_results, o_pvalue > .01 & o_lnFC > 0 &  flag_change_type %in% list("sssG","GssG","sGsG","ssGG","GGsG","GsGG","sGGG","GGGG"))

# ================================================================================================================================================================================================================
# Orangutan decreases
# ================================================================================================================================================================================================================
o_dec.both  <- subset(edgeR_results, o_pvalue < .01 & o_lnFC < 0 &  flag_change_type %in% list("sssL","LssL","sLsL","ssLL","LLsL","LsLL","sLLL","LLLL"))
o_dec.edgeR <- subset(edgeR_results, o_pvalue < .01 & o_lnFC < 0 & !flag_change_type %in% list("sssL","LssL","sLsL","ssLL","LLsL","LsLL","sLLL","LLLL"))
o_dec.glm   <- subset(edgeR_results, o_pvalue > .01 & o_lnFC < 0 &  flag_change_type %in% list("sssL","LssL","sLsL","ssLL","LLsL","LsLL","sLLL","LLLL"))

# ================================================================================================================================================================================================================
# Calculate counts
# ================================================================================================================================================================================================================

cat("Calculating counts ... \n")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Non differential sites
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count.non_differential.either      <- nrow(non_differential.both) + nrow(non_differential.edgeR) + nrow(non_differential.glm)
count.non_differential.both        <- nrow(non_differential.both)
count.non_differential.edgeR_only  <- nrow(non_differential.edgeR)
count.non_differential.glm_only    <- nrow(non_differential.glm)
count.non_differential.edgeR_total <- nrow(non_differential.both) + nrow(non_differential.edgeR)
count.non_differential.glm_total   <- nrow(non_differential.both) + nrow(non_differential.glm)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential sites
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count.differential.either      <- nrow(differential.both) + nrow(differential.edgeR) + nrow(differential.glm)
count.differential.both        <- nrow(differential.both)
count.differential.edgeR_only  <- nrow(differential.edgeR)
count.differential.glm_only    <- nrow(differential.glm)
count.differential.edgeR_total <- nrow(differential.both) + nrow(differential.edgeR)
count.differential.glm_total   <- nrow(differential.both) + nrow(differential.glm)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count.h_inc.either      <- nrow(h_inc.both) + nrow(h_inc.edgeR) + nrow(h_inc.glm)
count.h_inc.both        <- nrow(h_inc.both)
count.h_inc.edgeR_only  <- nrow(h_inc.edgeR)
count.h_inc.glm_only    <- nrow(h_inc.glm)
count.h_inc.edgeR_total <- nrow(h_inc.both) + nrow(h_inc.edgeR)
count.h_inc.glm_total   <- nrow(h_inc.both) + nrow(h_inc.glm)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count.h_dec.either      <- nrow(h_dec.both) + nrow(h_dec.edgeR) + nrow(h_dec.glm)
count.h_dec.both        <- nrow(h_dec.both)
count.h_dec.edgeR_only  <- nrow(h_dec.edgeR)
count.h_dec.glm_only    <- nrow(h_dec.glm)
count.h_dec.edgeR_total <- nrow(h_dec.both) + nrow(h_dec.edgeR)
count.h_dec.glm_total   <- nrow(h_dec.both) + nrow(h_dec.glm)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count.c_inc.either      <- nrow(c_inc.both) + nrow(c_inc.edgeR) + nrow(c_inc.glm)
count.c_inc.both        <- nrow(c_inc.both)
count.c_inc.edgeR_only  <- nrow(c_inc.edgeR)
count.c_inc.glm_only    <- nrow(c_inc.glm)
count.c_inc.edgeR_total <- nrow(c_inc.both) + nrow(c_inc.edgeR)
count.c_inc.glm_total   <- nrow(c_inc.both) + nrow(c_inc.glm)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count.c_dec.either      <- nrow(c_dec.both) + nrow(c_dec.edgeR) + nrow(c_dec.glm)
count.c_dec.both        <- nrow(c_dec.both)
count.c_dec.edgeR_only  <- nrow(c_dec.edgeR)
count.c_dec.glm_only    <- nrow(c_dec.glm)
count.c_dec.edgeR_total <- nrow(c_dec.both) + nrow(c_dec.edgeR)
count.c_dec.glm_total   <- nrow(c_dec.both) + nrow(c_dec.glm)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count.g_inc.either      <- nrow(g_inc.both) + nrow(g_inc.edgeR) + nrow(g_inc.glm)
count.g_inc.both        <- nrow(g_inc.both)
count.g_inc.edgeR_only  <- nrow(g_inc.edgeR)
count.g_inc.glm_only    <- nrow(g_inc.glm)
count.g_inc.edgeR_total <- nrow(g_inc.both) + nrow(g_inc.edgeR)
count.g_inc.glm_total   <- nrow(g_inc.both) + nrow(g_inc.glm)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count.g_dec.either      <- nrow(g_dec.both) + nrow(g_dec.edgeR) + nrow(g_dec.glm)
count.g_dec.both        <- nrow(g_dec.both)
count.g_dec.edgeR_only  <- nrow(g_dec.edgeR)
count.g_dec.glm_only    <- nrow(g_dec.glm)
count.g_dec.edgeR_total <- nrow(g_dec.both) + nrow(g_dec.edgeR)
count.g_dec.glm_total   <- nrow(g_dec.both) + nrow(g_dec.glm)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan increases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count.o_inc.either      <- nrow(o_inc.both) + nrow(o_inc.edgeR) + nrow(o_inc.glm)
count.o_inc.both        <- nrow(o_inc.both)
count.o_inc.edgeR_only  <- nrow(o_inc.edgeR)
count.o_inc.glm_only    <- nrow(o_inc.glm)
count.o_inc.edgeR_total <- nrow(o_inc.both) + nrow(o_inc.edgeR)
count.o_inc.glm_total   <- nrow(o_inc.both) + nrow(o_inc.glm)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan decreases
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
count.o_dec.either      <- nrow(o_dec.both) + nrow(o_dec.edgeR) + nrow(o_dec.glm)
count.o_dec.both        <- nrow(o_dec.both)
count.o_dec.edgeR_only  <- nrow(o_dec.edgeR)
count.o_dec.glm_only    <- nrow(o_dec.glm)
count.o_dec.edgeR_total <- nrow(o_dec.both) + nrow(o_dec.edgeR)
count.o_dec.glm_total   <- nrow(o_dec.both) + nrow(o_dec.glm)


# ################################################################################################################################################################################################################
# END SECTION 5: CREATE SUBSETS
# ################################################################################################################################################################################################################

# ################################################################################################################################################################################################################
# SECTION 6: CREATE PLOTS
#	1. Scatter plots for all sites (GLM beta value vs edgeR ln fold change value)
#		2 rows and 2 columns
#		row 1: human, chimpanzee
#		row 2: gorilla, orangutan
#		x-axis is the GLM beta value for the species
#		y-axis is the edgeR ln fold change value for the species
#	2. Density plots for increases (GLM beta values)
#		4 rows: human, chimpanzee, gorilla, orangutan
#		2 columns: GLM beta values for that species; GLM beta values for macaque
#		3 lines: both methods; edgeR only; glm only
#	3. Density plots for decreases (GLM beta values)
#		4 rows: human, chimpanzee, gorilla, orangutan
#		2 columns: GLM beta values for that species; GLM beta values for macaque
#		3 lines: both methods; edgeR only; glm only
#	4. Density plots for edgeR fold change and GLM beta value
#		4 rows: human, chimpanzee, gorilla, orangutan
#		2 columns: increases, decreases
#		2 lines: edgeR only; glm only
#	5. Scatter plots for increases (GLM beta value vs edgeR ln fold change value)
#		4 rows: human, chimpanzee, gorilla, orangutan
#		3 columns: both methods; edgeR only; GLM only;
#		x-axis is the GLM beta value for the species
#		y-axis is the edgeR ln fold change value for the species
#	6. Scatter plots for decreases (GLM beta value vs edgeR ln fold change value)
#		4 rows: human, chimpanzee, gorilla, orangutan
#		3 columns: both methods; edgeR only; GLM only;
#		x-axis is the GLM beta value for the species
#		y-axis is the edgeR ln fold change value for the species
#	7. Scatter plots for increases (GLM beta value for species vs macaque)
#		4 rows: human, chimpanzee, gorilla, orangutan
#		3 columns: both methods; edgeR only; GLM only;
#		x-axis is the GLM beta value for the species
#		y-axis is the GLM beta value for macaque
#	8. Scatter plots for decreases (GLM beta value for species vs macaque)
#		4 rows: human, chimpanzee, gorilla, orangutan
#		3 columns: both methods; edgeR only; GLM only;
#		x-axis is the GLM beta value for the species
#		y-axis is the GLM beta value for macaque
#	9. Venn diagrams of the counts of calls for glm only, edgeR only, both
#		5 rows x 2 columns
#			row 1: non differential sites; differential sites
#			row 2: human increases; human decreases
#			row 3: chimpanzee increases; chimpanzee decreases
#			row 4: gorilla increases; gorilla decreases
#			row 5: orangutan increases; orangutan decreases
#		Create 2 versions of the plot: one with count labels and one without
# ################################################################################################################################################################################################################

cat("Generating plots ... \n")

# ================================================================================================================================================================================================================
# Plot #1: Scatter plots for all sites (GLM beta value vs edgeR ln fold change value)
#	2 rows and 2 columns
#	row 1: human, chimpanzee
#	row 2: gorilla, orangutan
#	x-axis is the GLM beta value for the species
#	y-axis is the edgeR ln fold change value for the species
# ================================================================================================================================================================================================================

png("scatter_plots.GLM_beta_values_vs_edgeR_ln_fold_change_values.all_sites.png", height=1900, width=1600, res=300)

par(mfrow=c(2,2))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(edgeR_results$Bh, edgeR_results$h_lnFC, main="All DHS sites\n Human values",
	xlab="GLM", ylab="edgeR", xlim=c(-6,6), ylim=c(-6,6), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(edgeR_results$Bc, edgeR_results$c_lnFC, main="All DHS sites\n Chimpanzee values",
	xlab="GLM", ylab="edgeR", xlim=c(-6,6), ylim=c(-6,6), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(edgeR_results$Bg, edgeR_results$g_lnFC, main="All DHS sites\n Gorilla values",
	xlab="GLM", ylab="edgeR", xlim=c(-6,6), ylim=c(-6,6), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(edgeR_results$Bo, edgeR_results$o_lnFC, main="All DHS sites\n Orangutan values",
	xlab="GLM", ylab="edgeR", xlim=c(-6,6), ylim=c(-6,6), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev.off()

# ================================================================================================================================================================================================================
# Plot #2: Density plots for increases (GLM beta values)
#	4 rows: human, chimpanzee, gorilla, orangutan
#	2 columns: GLM beta values for that species; GLM beta values for macaque
#	3 lines: both methods; edgeR only; glm only
# ================================================================================================================================================================================================================

png("density_plots.GLM_beta_values.increases.png", height=3000, width=2400, res=300)

par(mfrow=c(4,2))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(h_inc.both$Bh), main="Increases in human (human beta values)\n green=both methods (6,842 sites)\n red=edgeR only (1,174 sites); blue=GLM only (3,246 sites)",
	xlab="GLM beta value", col="green", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(h_inc.edgeR$Bh), col="red")
lines(density(h_inc.glm$Bh), col="blue")

plot(density(h_inc.both$Bm), main="Increases in human (macaque beta values)\n green=both methods (6,842 sites)\n red=edgeR only (1,174 sites); blue=GLM only (3,246 sites)",
	xlab="GLM beta value", col="green", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(h_inc.edgeR$Bm), col="red")
lines(density(h_inc.glm$Bm), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(c_inc.both$Bc), main="Increases in chimpanzee (chimpanzee beta values)\n green=both methods (5,747 sites)\n red=edgeR only (1,099 sites); blue=GLM only (3,792 sites)",
	xlab="GLM beta value", col="green", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(c_inc.edgeR$Bc), col="red")
lines(density(c_inc.glm$Bc), col="blue")

plot(density(c_inc.both$Bm), main="Increases in chimpanzee (macaque beta values)\n green=both methods (5,747 sites)\n red=edgeR only (1,099 sites); blue=GLM only (3,792 sites)",
	xlab="GLM beta value", col="green", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(c_inc.edgeR$Bm), col="red")
lines(density(c_inc.glm$Bm), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(g_inc.both$Bg), main="Increases in gorilla (gorilla beta values)\n green=both methods (6,019 sites)\n red=edgeR only (988 sites); blue=GLM only (5,234 sites)",
	xlab="GLM beta value", col="green", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(g_inc.edgeR$Bg), col="red")
lines(density(g_inc.glm$Bg), col="blue")

plot(density(g_inc.both$Bm), main="Increases in gorilla (macaque beta values)\n green=both methods (6,019 sites)\n red=edgeR only (988 sites); blue=GLM only (5,234 sites)",
	xlab="GLM beta value", col="green", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(g_inc.edgeR$Bm), col="red")
lines(density(g_inc.glm$Bm), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(o_inc.both$Bo), main="Increases in orangutan (orangutan beta values)\n green=both methods (5,876 sites)\n red=edgeR only (1,094 sites); blue=GLM only (4,693 sites)",
	xlab="GLM beta value", col="green", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(o_inc.edgeR$Bo), col="red")
lines(density(o_inc.glm$Bo), col="blue")

plot(density(o_inc.both$Bm), main="Increases in orangutan (macaque beta values)\n green=both methods (5,876 sites)\n red=edgeR only (1,094 sites); blue=GLM only (4,693 sites)",
	xlab="GLM beta value", col="green", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(o_inc.edgeR$Bm), col="red")
lines(density(o_inc.glm$Bm), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev.off()

# ================================================================================================================================================================================================================
# Plot #3: Density plots for decreases (GLM beta values)
#	4 rows: human, chimpanzee, gorilla, orangutan
#	2 columns: GLM beta values for that species; GLM beta values for macaque
#	3 lines: both methods; edgeR only; glm only
# ================================================================================================================================================================================================================

png("density_plots.GLM_beta_values.decreases.png", height=3000, width=2400, res=300)

par(mfrow=c(4,2))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(h_dec.both$Bh), main="Decreases in human (human beta values)\n green=both methods (5,026 sites)\n red=edgeR only (805 sites); blue=GLM only (3,983 sites)",
	xlab="GLM beta value", col="green", xlim=c(-6,0), ylim=c(0,1.4), cex.main=1)
lines(density(h_dec.edgeR$Bh), col="red")
lines(density(h_dec.glm$Bh), col="blue")

plot(density(h_dec.both$Bm), main="Decreases in human (macaque beta values)\n green=both methods (5,026 sites)\n red=edgeR only (805 sites); blue=GLM only (3,983 sites)",
	xlab="GLM beta value", col="green", xlim=c(0,6), ylim=c(0,1.4), cex.main=1)
lines(density(h_dec.edgeR$Bm), col="red")
lines(density(h_dec.glm$Bm), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(c_dec.both$Bc), main="Decreases in chimpanzee (chimpanzee beta values)\n green=both methods (5,343 sites)\n red=edgeR only (799 sites); blue=GLM only (3,820 sites)",
	xlab="GLM beta value", col="green", xlim=c(-6,0), ylim=c(0,1.4), cex.main=1)
lines(density(c_dec.edgeR$Bc), col="red")
lines(density(c_dec.glm$Bc), col="blue")

plot(density(c_dec.both$Bm), main="Decreases in chimpanzee (macaque beta values)\n green=both methods (5,343 sites)\n red=edgeR only (799 sites); blue=GLM only (3,820 sites)",
	xlab="GLM beta value", col="green", xlim=c(0,6), ylim=c(0,1.4), cex.main=1)
lines(density(c_dec.edgeR$Bm), col="red")
lines(density(c_dec.glm$Bm), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(g_dec.both$Bg), main="Decreases in gorilla (gorilla beta values)\n green=both methods (5,119 sites)\n red=edgeR only (1,129 sites); blue=GLM only (3,428 sites)",
	xlab="GLM beta value", col="green", xlim=c(-6,0), ylim=c(0,1.4), cex.main=1)
lines(density(g_dec.edgeR$Bg), col="red")
lines(density(g_dec.glm$Bg), col="blue")

plot(density(g_dec.both$Bm), main="Decreases in gorilla (macaque beta values)\n green=both methods (5,119 sites)\n red=edgeR only (1,129 sites); blue=GLM only (3,428 sites)",
	xlab="GLM beta value", col="green", xlim=c(0,6), ylim=c(0,1.4), cex.main=1)
lines(density(g_dec.edgeR$Bm), col="red")
lines(density(g_dec.glm$Bm), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(o_dec.both$Bo), main="Decreases in orangutan (orangutan beta values)\n green=both methods (5,569 sites)\n red=edgeR only (1,109 sites); blue=GLM only (2,845 sites)",
	xlab="GLM beta value", col="green", xlim=c(-6,0), ylim=c(0,1.4), cex.main=1)
lines(density(o_dec.edgeR$Bo), col="red")
lines(density(o_dec.glm$Bo), col="blue")

plot(density(o_dec.both$Bm), main="Decreases in orangutan (macaque beta values)\n green=both methods (5,569 sites)\n red=edgeR only (1,109 sites); blue=GLM only (2,845 sites)",
	xlab="GLM beta value", col="green", xlim=c(0,6), ylim=c(0,1.4), cex.main=1)
lines(density(o_dec.edgeR$Bm), col="red")
lines(density(o_dec.glm$Bm), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev.off()

# ================================================================================================================================================================================================================
# Plot #4: Density plots for edgeR fold change and GLM beta value
#	4 rows: human, chimpanzee, gorilla, orangutan
#	2 columns: increases, decreases
#	2 lines: edgeR only; glm only
# ================================================================================================================================================================================================================

png("density_plots.edgeR_log_fold_and_GLM_beta_values.all.png", height=3000, width=2400, res=300)

par(mfrow=c(4,2))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(h_inc.edgeR$Bh), main="Increases in human", xlab="beta value", col="red", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(h_inc.glm$Bh), col="blue")

plot(density(h_dec.edgeR$Bh), main="Decreases in human", xlab="beta value", col="red", xlim=c(-6,0), ylim=c(0,1.4), cex.main=1)
lines(density(h_dec.glm$Bh), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(c_inc.edgeR$Bc), main="Increases in chimpanzee", xlab="beta value", col="red", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(c_inc.glm$Bc), col="blue")

plot(density(c_dec.edgeR$Bc), main="Decreases in chimpanzee", xlab="beta value", col="red", xlim=c(-6,0), ylim=c(0,1.4), cex.main=1)
lines(density(c_dec.glm$Bc), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(g_inc.edgeR$Bg), main="Increases in gorilla", xlab="beta value", col="red", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(g_inc.glm$Bg), col="blue")

plot(density(g_dec.edgeR$Bg), main="Decreases in gorilla", xlab="beta value", col="red", xlim=c(-6,0), ylim=c(0,1.4), cex.main=1)
lines(density(g_dec.glm$Bg), col="blue")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(density(o_inc.edgeR$Bo), main="Increases in orangutan", xlab="beta value", col="red", xlim=c(-1,5), ylim=c(0,1), cex.main=1)
lines(density(o_inc.glm$Bo), col="blue")

plot(density(o_dec.edgeR$Bo), main="Decreases in orangutan", xlab="beta value", col="red", xlim=c(-6,0), ylim=c(0,1.4), cex.main=1)
lines(density(o_dec.glm$Bo), col="blue")


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev.off()


# ================================================================================================================================================================================================================
# Plot #5: Scatter plots for increases (GLM beta value vs edgeR ln fold change value)
#		4 rows: human, chimpanzee, gorilla, orangutan
#		3 columns: both methods; edgeR only; GLM only;
#		x-axis is the GLM beta value for the species
#		y-axis is the edgeR ln fold change value for that species
# ================================================================================================================================================================================================================

png("scatter_plots.GLM_beta_values_vs_edgeR_ln_fold_change_values.increases.png", height=3000, width=2400, res=300)

par(mfrow=c(4,3))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(h_inc.both$Bh, h_inc.both$h_lnFC, main="Increases in human\n both methods (6,842 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

plot(h_inc.edgeR$Bh, h_inc.edgeR$h_lnFC, main="Increases in human\n edgeR only (1,174 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

plot(h_inc.glm$Bh, h_inc.glm$h_lnFC, main="Increases in human\n GLM only (3,246 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(c_inc.both$Bc, c_inc.both$c_lnFC, main="Increases in chimpanzee\n both methods (5,747 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

plot(c_inc.edgeR$Bc, c_inc.edgeR$c_lnFC, main="Increases in chimpanzee\n edgeR only (1,099 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

plot(c_inc.glm$Bc, c_inc.glm$c_lnFC, main="Increases in chimpanzee\n GLM only (3,792 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(g_inc.both$Bg, g_inc.both$g_lnFC, main="Increases in gorilla\n both methods (6,019 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

plot(g_inc.edgeR$Bg, g_inc.edgeR$g_lnFC, main="Increases in gorilla\n edgeR only (988 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

plot(g_inc.glm$Bg, g_inc.glm$g_lnFC, main="Increases in gorilla\n GLM only (5,234 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(o_inc.both$Bo, o_inc.both$o_lnFC, main="Increases in orangutan\n both methods (5,876 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

plot(o_inc.edgeR$Bo, o_inc.edgeR$o_lnFC, main="Increases in orangutan\n edgeR only (1,094 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

plot(o_inc.glm$Bo, o_inc.glm$o_lnFC, main="Increases in orangutan\n GLM only (4,693 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(0,6), ylim=c(0,6), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev.off()

# ================================================================================================================================================================================================================
# Plot #6: Scatter plots for decreases (GLM beta value vs edgeR ln fold change value)
#		4 rows: human, chimpanzee, gorilla, orangutan
#		3 columns: both methods; edgeR only; GLM only;
#		x-axis is the GLM beta value for the species
#		y-axis is the edgeR ln fold change value for that species
# ================================================================================================================================================================================================================

png("scatter_plots.GLM_beta_values_vs_edgeR_ln_fold_change_values.decreases.png", height=3000, width=2400, res=300)

par(mfrow=c(4,3))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(h_dec.both$Bh, h_dec.both$h_lnFC, main="Decreases in human\n both methods (5,026 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

plot(h_dec.edgeR$Bh, h_dec.edgeR$h_lnFC, main="Decreases in human\n edgeR only (805 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

plot(h_dec.glm$Bh, h_dec.glm$h_lnFC, main="Decreases in human\n GLM only (3,983 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(c_dec.both$Bc, c_dec.both$c_lnFC, main="Decreases in chimpanzee\n both methods (5,343 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

plot(c_dec.edgeR$Bc, c_dec.edgeR$c_lnFC, main="Decreases in chimpanzee\n edgeR only (799 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

plot(c_dec.glm$Bc, c_dec.glm$c_lnFC, main="Decreases in chimpanzee\n GLM only (3,820 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(g_dec.both$Bg, g_dec.both$g_lnFC, main="Decreases in gorilla\n both methods (5,119 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

plot(g_dec.edgeR$Bg, g_dec.edgeR$g_lnFC, main="Decreases in gorilla\n edgeR only (1,129 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

plot(g_dec.glm$Bg, g_dec.glm$g_lnFC, main="Decreases in gorilla\n GLM only (3,428 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(o_dec.both$Bo, o_dec.both$o_lnFC, main="Decreases in orangutan\n both methods (5,569 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

plot(o_dec.edgeR$Bo, o_dec.edgeR$o_lnFC, main="Decreases in orangutan\n edgeR only (1,109 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

plot(o_dec.glm$Bo, o_dec.glm$o_lnFC, main="Decreases in orangutan\n GLM only (2,845 sites)",
	xlab="GLM", ylab="edgeR", xlim=c(-6,0), ylim=c(-6,0), cex.main=1)
abline(0,1,col="red")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev.off()

# ================================================================================================================================================================================================================
# Plot #7: Scatter plots for increases (GLM beta value for species vs macaque)
#		4 rows: human, chimpanzee, gorilla, orangutan
#		3 columns: both methods; edgeR only; GLM only;
#		x-axis is the GLM beta value for the species
#		y-axis is the GLM beta value for macaque
# ================================================================================================================================================================================================================

png("scatter_plots.GLM_beta_values.species_vs_macaque.increases.png", height=3000, width=1900, res=300)

par(mfrow=c(4,3))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(h_inc.both$Bh, h_inc.both$Bm, xlab="GLM beta value (human)", main="Increases in human\n both methods (6,842 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

plot(h_inc.edgeR$Bh, h_inc.edgeR$Bm, xlab="GLM beta value (human)", main="Increases in human\n edgeR only (1,174 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

plot(h_inc.glm$Bh, h_inc.glm$Bm, xlab="GLM beta value (human)", main="Increases in human\n GLM only (3,246 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(c_inc.both$Bc, c_inc.both$Bm, xlab="GLM beta value (chimpanzee)", main="Increases in chimpanzee\n both methods (5,747 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

plot(c_inc.edgeR$Bc, c_inc.edgeR$Bm, xlab="GLM beta value (chimpanzee)", main="Increases in chimpanzee\n edgeR only (1,099 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

plot(c_inc.glm$Bc, c_inc.glm$Bm, xlab="GLM beta value (chimpanzee)", main="Increases in chimpanzee\n GLM only (3,792 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(g_inc.both$Bg, g_inc.both$Bm, xlab="GLM beta value (gorilla)", main="Increases in gorilla\nboth methods (6,019 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

plot(g_inc.edgeR$Bg, g_inc.edgeR$Bm, xlab="GLM beta value (gorilla)", main="Increases in gorilla\n edgeR only (988 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

plot(g_inc.glm$Bg, g_inc.glm$Bm, xlab="GLM beta value (gorilla)", main="Increases in gorilla\n GLM only (5,234 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(o_inc.both$Bo, o_inc.both$Bm, xlab="GLM beta value (orangutan)", main="Increases in orangutan\n both methods (5,876 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

plot(o_inc.edgeR$Bo, o_inc.edgeR$Bm, xlab="GLM beta value (orangutan)", main="Increases in orangutan\n edgeR only (1,094 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

plot(o_inc.glm$Bo, o_inc.glm$Bm, xlab="GLM beta value (orangutan)", main="Increases in orangutan\n GLM only (4,693 sites)",
	ylab="GLM beta value (macaque)", xlim=c(0,6), ylim=c(-1,5), cex.main=1, pch=".")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev.off()

# ================================================================================================================================================================================================================
# Plot #8: Scatter plots for decreases (GLM beta value for species vs macaque)
#		4 rows: human, chimpanzee, gorilla, orangutan
#		3 columns: both methods; edgeR only; GLM only;
#		x-axis is the GLM beta value for the species
#		y-axis is the GLM beta value for macaque
# ================================================================================================================================================================================================================

png("scatter_plots.GLM_beta_values.species_vs_macaque.decreases.png", height=3000, width=1900, res=300)

par(mfrow=c(4,3))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Human
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(h_dec.both$Bh, h_dec.both$Bm, xlab="GLM beta value (human)", main="Decreases in human\n both methods (5,026 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

plot(h_dec.edgeR$Bh, h_dec.edgeR$Bm, xlab="GLM beta value (human)", main="Decreases in human\n edgeR only (805 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

plot(h_dec.glm$Bh, h_dec.glm$Bm, xlab="GLM beta value (human)", main="Decreases in human\n GLM only (3,983 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chimpanzee
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(c_dec.both$Bc, c_dec.both$Bm, xlab="GLM beta value (chimpanzee)", main="Decreases in chimpanzee\n both methods (5,343 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

plot(c_dec.edgeR$Bc, c_dec.edgeR$Bm, xlab="GLM beta value (chimpanzee)", main="Decreases in chimpanzee\n edgeR only (799 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

plot(c_dec.glm$Bc, c_dec.glm$Bm, xlab="GLM beta value (chimpanzee)", main="Decreases in chimpanzee\n GLM only (3,820 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gorilla
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(g_dec.both$Bg, g_dec.both$Bm, xlab="GLM beta value (gorilla)", main="Decreases in gorilla\n both methods (5,119 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

plot(g_dec.edgeR$Bg, g_dec.edgeR$Bm, xlab="GLM beta value (gorilla)", main="Decreases in gorilla\n edgeR only (1,129 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

plot(g_dec.glm$Bg, g_dec.glm$Bm, xlab="GLM beta value (gorilla)", main="Decreases in gorilla\n GLM only (3,428 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Orangutan
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(o_dec.both$Bo, o_dec.both$Bm, xlab="GLM beta value (orangutan)", main="Decreases in orangutan\n both methods (5,569 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

plot(o_dec.edgeR$Bo, o_dec.edgeR$Bm, xlab="GLM beta value (orangutan)", main="Decreases in orangutan\n edgeR only (1,109 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

plot(o_dec.glm$Bo, o_dec.glm$Bm, xlab="GLM beta value (orangutan)", main="Decreases in orangutan\n GLM only (2,845 sites)",
	ylab="GLM beta value (macaque)", xlim=c(-6,0), ylim=c(0,6), cex.main=1, pch=".")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dev.off()


# ================================================================================================================================================================================================================
# Plot #9: Venn diagrams of the counts of calls for glm only, edgeR only, both
#	5 rows x 2 columns
#		row 1: non differential sites; differential sites
#		row 2: human increases; human decreases
#		row 3: chimpanzee increases; chimpanzee decreases
#		row 4: gorilla increases; gorilla decreases
#		row 5: orangutan increases; orangutan decreases
#	Create 2 versions of the plot: one with count labels and one without
# ================================================================================================================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate the version with count labels
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

png("venn_diagrams.increases_and_decreases.with_count_labels.png", height=3000, width=1200, res=300)

venn.non_differential <- draw.pairwise.venn(count.non_differential.edgeR_total, count.non_differential.glm_total, count.non_differential.both,
	euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=.5, ind=FALSE, margin=.1)

venn.differential <- draw.pairwise.venn(count.differential.edgeR_total, count.differential.glm_total, count.differential.both, 
	euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=.5, ind=FALSE, margin=.1)

venn.h_inc <- draw.pairwise.venn(count.h_inc.edgeR_total, count.h_inc.glm_total, count.h_inc.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=.5, ind=FALSE, margin=.1)
venn.c_inc <- draw.pairwise.venn(count.c_inc.edgeR_total, count.c_inc.glm_total, count.c_inc.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=.5, ind=FALSE, margin=.1)
venn.g_inc <- draw.pairwise.venn(count.g_inc.edgeR_total, count.g_inc.glm_total, count.g_inc.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=.5, ind=FALSE, margin=.1)
venn.o_inc <- draw.pairwise.venn(count.o_inc.edgeR_total, count.o_inc.glm_total, count.o_inc.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=.5, ind=FALSE, margin=.1)

venn.h_dec <- draw.pairwise.venn(count.h_dec.edgeR_total, count.h_dec.glm_total, count.h_dec.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=.5, ind=FALSE, margin=.1)
venn.c_dec <- draw.pairwise.venn(count.c_dec.edgeR_total, count.c_dec.glm_total, count.c_dec.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=.5, ind=FALSE, margin=.1)
venn.g_dec <- draw.pairwise.venn(count.g_dec.edgeR_total, count.g_dec.glm_total, count.g_dec.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=.5, ind=FALSE, margin=.1)
venn.o_dec <- draw.pairwise.venn(count.o_dec.edgeR_total, count.o_dec.glm_total, count.o_dec.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=.5, ind=FALSE, margin=.1)

grid.arrange(gTree(children=venn.non_differential), gTree(children=venn.differential),
	     gTree(children=venn.h_inc), gTree(children=venn.h_dec), gTree(children=venn.c_inc), gTree(children=venn.c_dec),
	     gTree(children=venn.g_inc), gTree(children=venn.g_dec), gTree(children=venn.o_inc), gTree(children=venn.o_dec), nrow=5)

dev.off()

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate the version without count labels
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

png("venn_diagrams.increases_and_decreases.without_count_labels.png", height=3000, width=1200, res=300)

venn.non_differential <- draw.pairwise.venn(count.non_differential.edgeR_total, count.non_differential.glm_total, count.non_differential.both, 
	euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=0, ind=FALSE, margin=.1)

venn.differential <- draw.pairwise.venn(count.differential.edgeR_total, count.differential.glm_total, count.differential.both, 
	euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=0, ind=FALSE, margin=.1)

venn.h_inc <- draw.pairwise.venn(count.h_inc.edgeR_total, count.h_inc.glm_total, count.h_inc.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=0, ind=FALSE, margin=.1)
venn.c_inc <- draw.pairwise.venn(count.c_inc.edgeR_total, count.c_inc.glm_total, count.c_inc.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=0, ind=FALSE, margin=.1)
venn.g_inc <- draw.pairwise.venn(count.g_inc.edgeR_total, count.g_inc.glm_total, count.g_inc.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=0, ind=FALSE, margin=.1)
venn.o_inc <- draw.pairwise.venn(count.o_inc.edgeR_total, count.o_inc.glm_total, count.o_inc.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=0, ind=FALSE, margin=.1)

venn.h_dec <- draw.pairwise.venn(count.h_dec.edgeR_total, count.h_dec.glm_total, count.h_dec.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=0, ind=FALSE, margin=.1)
venn.c_dec <- draw.pairwise.venn(count.c_dec.edgeR_total, count.c_dec.glm_total, count.c_dec.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=0, ind=FALSE, margin=.1)
venn.g_dec <- draw.pairwise.venn(count.g_dec.edgeR_total, count.g_dec.glm_total, count.g_dec.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=0, ind=FALSE, margin=.1)
venn.o_dec <- draw.pairwise.venn(count.o_dec.edgeR_total, count.o_dec.glm_total, count.o_dec.both, euler.d=TRUE, scaled=TRUE, fill=c("red","blue"), alpha=c(0.5,0.3), cex=0, ind=FALSE, margin=.1)

grid.arrange(gTree(children=venn.non_differential), gTree(children=venn.differential),
	     gTree(children=venn.h_inc), gTree(children=venn.h_dec), gTree(children=venn.c_inc), gTree(children=venn.c_dec),
	     gTree(children=venn.g_inc), gTree(children=venn.g_dec), gTree(children=venn.o_inc), gTree(children=venn.o_dec), nrow=5)


dev.off()

# ################################################################################################################################################################################################################
# END SECTION 6: CREATE PLOTS
# ################################################################################################################################################################################################################



# ################################################################################################################################################################################################################
# SECTION 7: FINISH UP
# ################################################################################################################################################################################################################

# ================================================================================================================================================================================================================
# Print counts
# ================================================================================================================================================================================================================

cat("=================================================================================================\n")
cat("Counts\n")
cat("=================================================================================================\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Non differential sites\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("      Either method\t",   count.non_differential.either, "\n")
cat("      Both methods\t",    count.non_differential.both, "\n")
cat("      edgeR only\t",      count.non_differential.edgeR_only, "\n")
cat("      GLM only\t\t",      count.non_differential.glm_only, "\n")
cat("      Total edgeR\t\t",   count.non_differential.edgeR_total, "\n")
cat("      Total GLM\t\t\t",   count.non_differential.glm_total, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Differential sites\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("      Either method\t",   count.differential.either, "\n")
cat("      Both methods\t",    count.differential.both, "\n")
cat("      edgeR only\t",      count.differential.edgeR_only, "\n")
cat("      GLM only\t\t",      count.differential.glm_only, "\n")
cat("      Total edgeR\t\t",   count.differential.edgeR_total, "\n")
cat("      Total GLM\t\t\t",   count.differential.glm_total, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Human increases\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("      Either method\t",   count.h_inc.either, "\n")
cat("      Both methods\t",    count.h_inc.both, "\n")
cat("      edgeR only\t",      count.h_inc.edgeR_only, "\n")
cat("      GLM only\t\t",      count.h_inc.glm_only, "\n")
cat("      Total edgeR\t\t",   count.h_inc.edgeR_total, "\n")
cat("      Total GLM\t\t\t",   count.h_inc.glm_total, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Human decreases\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("      Either method\t",   count.h_dec.either, "\n")
cat("      Both methods\t",    count.h_dec.both, "\n")
cat("      edgeR only\t",      count.h_dec.edgeR_only, "\n")
cat("      GLM only\t\t",      count.h_dec.glm_only, "\n")
cat("      Total edgeR\t\t",   count.h_dec.edgeR_total, "\n")
cat("      Total GLM\t\t\t",   count.h_dec.glm_total, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Chimpanzee increases\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("      Either method\t",   count.c_inc.either, "\n")
cat("      Both methods\t",    count.c_inc.both, "\n")
cat("      edgeR only\t",      count.c_inc.edgeR_only, "\n")
cat("      GLM only\t\t",      count.c_inc.glm_only, "\n")
cat("      Total edgeR\t\t",   count.c_inc.edgeR_total, "\n")
cat("      Total GLM\t\t\t",   count.c_inc.glm_total, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Chimpanzee decreases\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("      Either method\t",   count.c_dec.either, "\n")
cat("      Both methods\t",    count.c_dec.both, "\n")
cat("      edgeR only\t",      count.c_dec.edgeR_only, "\n")
cat("      GLM only\t\t",      count.c_dec.glm_only, "\n")
cat("      Total edgeR\t\t",   count.c_dec.edgeR_total, "\n")
cat("      Total GLM\t\t\t",   count.c_dec.glm_total, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Gorilla increases\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("      Either method\t",   count.g_inc.either, "\n")
cat("      Both methods\t",    count.g_inc.both, "\n")
cat("      edgeR only\t",      count.g_inc.edgeR_only, "\n")
cat("      GLM only\t\t",      count.g_inc.glm_only, "\n")
cat("      Total edgeR\t\t",   count.g_inc.edgeR_total, "\n")
cat("      Total GLM\t\t\t",   count.g_inc.glm_total, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Gorilla decreases\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("      Either method\t",   count.g_dec.either, "\n")
cat("      Both methods\t",    count.g_dec.both, "\n")
cat("      edgeR only\t",      count.g_dec.edgeR_only, "\n")
cat("      GLM only\t\t",      count.g_dec.glm_only, "\n")
cat("      Total edgeR\t\t",   count.g_dec.edgeR_total, "\n")
cat("      Total GLM\t\t\t",   count.g_dec.glm_total, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Orangutan increases\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("      Either method\t",   count.o_inc.either, "\n")
cat("      Both methods\t",    count.o_inc.both, "\n")
cat("      edgeR only\t",      count.o_inc.edgeR_only, "\n")
cat("      GLM only\t\t",      count.o_inc.glm_only, "\n")
cat("      Total edgeR\t\t",   count.o_inc.edgeR_total, "\n")
cat("      Total GLM\t\t\t",   count.o_inc.glm_total, "\n\n")

cat("-------------------------------------------------------------------------------------------------\n")
cat("Orangutan decreases\n")
cat("-------------------------------------------------------------------------------------------------\n")
cat("      Either method\t",   count.o_dec.either, "\n")
cat("      Both methods\t",    count.o_dec.both, "\n")
cat("      edgeR only\t",      count.o_dec.edgeR_only, "\n")
cat("      GLM only\t\t",      count.o_dec.glm_only, "\n")
cat("      Total edgeR\t\t",   count.o_dec.edgeR_total, "\n")
cat("      Total GLM\t\t\t",   count.o_dec.glm_total, "\n\n")


# ================================================================================================================================================================================================================
# Finish up
# ================================================================================================================================================================================================================

cat("-------------------------------------------------------------------------------------------------\n")
cat("END:", date(), "\n")
cat("-------------------------------------------------------------------------------------------------\n")

sink()

# ################################################################################################################################################################################################################
# END SECTION 7: CALCULATE COUNTS AND FINISH UP
# ################################################################################################################################################################################################################
