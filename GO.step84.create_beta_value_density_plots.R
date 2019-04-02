# ################################################################################################################################################################################################################
# GO.step84.create_beta_value_density_plots.R
#
# CC-BY Lee Edsall
#	email: le49@duke.edu
#	Twitter: @LeeEdsall
#
# This script was used to analyze data for Edsall et al. 2019 which compared DNase-seq data from 5 primates (human, chimpanzee, gorilla, orangutan, macaque)

# Create density plots of the beta values for the non-macaque species for each category of site
# ################################################################################################################################################################################################################


# theme elements to add to make a version of the legend for the ppt
# legend.title=element_text(size=14), legend.text=element_text(size=14, face="plain"), legend.background=element_rect(colour='black', fill='white', linetype='solid')

# ==========================================================================================================================================================================
# Load libraries
# ==========================================================================================================================================================================
library(ggplot2)

# ==========================================================================================================================================================================
# Set variables
# ==========================================================================================================================================================================
# Order is human, chimpanzee, gorilla, orangutan
# Colors checked for color blindness

my_pdf_line_colors <- c("black", "green3", "blue", "orange")
my_svg_line_colors <- c("black", "green3", "blue", "orange")

my_pdf_line_types <- c("solid", "longdash", "dashed", "solid")
my_svg_line_types <- c("solid", "longdash", "dashed", "solid")


my_pdf_line_widths <- c(2.5, 2.0, 1.5, 1.0)
my_svg_line_widths <- c(5.0, 4.5, 4.0, 3.0)

my_svg_text_size <- 45

# ==========================================================================================================================================================================
# Open connection to log file; write data
# ==========================================================================================================================================================================
sink(file="LOG.step85.create_beta_value_density_plots.R", append=FALSE, type="output", split=TRUE)
cat("---------------------------------------------------------------------------------------\n")
cat("BEGIN:", date(), "\n")
cat("---------------------------------------------------------------------------------------\n\n")


# ==========================================================================================================================================================================
# Load the glm results text file
# The scores have not been normalized
# There are 75 columns. The relevant columns are:
#	22: human beta value
#	25: chimp beta value
#	28: gorilla beta value
#	31: orangutan beta value
#	60: flag_change_type
# ==========================================================================================================================================================================
cat("Reading glm results file (all_DHS_sites.glm_analysis.all_sites.all_information.txt) ... ")

glm_results <- as.data.frame(read.table("all_DHS_sites.glm_analysis.all_sites.all_information.txt", sep="\t", header=TRUE))

rownames(glm_results) <- paste(glm_results$chrom, glm_results$start, glm_results$end, sep=":")

cat("done\n\n")

# ==========================================================================================================================================================================
# Generate the beta value density plots
#	Save all of them in one a pdf file and each one in its own svg file
#		The pdf file will contain a legend for each plot
#		The svg files will exclude the legends
# ==========================================================================================================================================================================

pdf("beta_value_density_plot.all.pdf", onefile=TRUE, paper="letter", bg="white", pagecentre=TRUE)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 1 of 31
# Non differential
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for sites that are non differential ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="ssss")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="ssss")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="ssss")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="ssss")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Non Differential\n(53,078 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_ssss <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Non Differential\n(53,078 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 2 of 31
# Human increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="Gsss")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="Gsss")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="Gsss")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="Gsss")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human Increases\n(2,539 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_Gsss <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human Increases\n(2,539 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 3 of 31
# Human decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="Lsss")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="Lsss")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="Lsss")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="Lsss")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human Decreases\n(209 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_Lsss <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human Decreases\n(209 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 4 of 31
# Chimpanzee increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for chimpanzee increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="sGss")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="sGss")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="sGss")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="sGss")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Chimpanzee Increases\n(1,539 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_sGss <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Chimpanzee Increases\n(1,539 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 5 of 31
# Chimpanzee decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for chimpanzee decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="sLss")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="sLss")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="sLss")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="sLss")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Chimpanzee Decreases\n(196 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_sLss <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Chimpanzee Decreases\n(196 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 6 of 31
# Gorilla increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for gorilla increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="ssGs")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="ssGs")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="ssGs")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="ssGs")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Gorilla Increases\n(2,737 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_ssGs <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Gorilla Increases\n(2,737 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 7 of 31
# Gorilla decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for gorilla decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="ssLs")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="ssLs")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="ssLs")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="ssLs")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Gorilla Decreases\n(254 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_ssLs <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Gorilla Decreases\n(254 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 8 of 31
# Orangutan increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for orangutan increases... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="sssG")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="sssG")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="sssG")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="sssG")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Orangutan Increases\n(3,759 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_sssG <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Orangutan Increases\n(3,759 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 9 of 31
# Orangutan decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for orangutan decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="sssL")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="sssL")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="sssL")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="sssL")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Orangutan Decreases\n(646 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_sssL <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Orangutan Decreases\n(646 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 10 of 31
# Macaque increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for macaque increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="LLLL")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="LLLL")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="LLLL")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="LLLL")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Macaque Increases\n(4,992 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_LLLL <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Macaque Increases\n(4,992 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 11 of 31
# Macaque decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for macaque decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="GGGG")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="GGGG")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="GGGG")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="GGGG")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Macaque Decreases\n(1,317 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_GGGG <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Macaque Decreases\n(1,317 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 12 of 31
# Human/chimpanzee increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/chimpanzee increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="GGss")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="GGss")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="GGss")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="GGss")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Chimpanzee Increases\n(1,735 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_GGss <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Chimpanzee Increases\n(1,735 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 13 of 31
# Human/chimpanzee decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/chimpanzee decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="LLss")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="LLss")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="LLss")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="LLss")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Chimpanzee Decreases\n(814 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_LLss <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Chimpanzee Decreases\n(814 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 14 of 31
# Human/gorilla increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/gorilla increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="GsGs")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="GsGs")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="GsGs")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="GsGs")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Gorilla Increases\n(997 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_GsGs <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Gorilla Increases\n(997 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 15 of 31
# Human/gorilla decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/gorilla decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="LsLs")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="LsLs")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="LsLs")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="LsLs")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Gorilla Decreases\n(221 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_LsLs <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Gorilla Decreases\n(221 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 16 of 31
# Human/orangutan increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/orangutan increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="GssG")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="GssG")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="GssG")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="GssG")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Orangutan Increases\n(997 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_GssG <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Orangutan Increases\n(997 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 17 of 31
# Human/orangutan decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/orangutan decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="LssL")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="LssL")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="LssL")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="LssL")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Orangutan Decreases\n(424 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_LssL <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Orangutan Decreases\n(424 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 18 of 31
# Chimpanzee/gorilla increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for chimpanzee and gorilla increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="sGGs")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="sGGs")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="sGGs")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="sGGs")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Chimpanzee/Gorilla Increases\n(1,568 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_sGGs <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Chimpanzee/Gorilla Increases\n(1,568 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 19 of 31
# Chimpanzee/gorilla decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for chimpanzee/gorilla decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="sLLs")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="sLLs")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="sLLs")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="sLLs")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Chimpanzee/Gorilla Decreases\n(387 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_sLLs <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Chimpanzee/Gorilla Decreases\n(387 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 20 of 31
# Chimpanzee/orangutan increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for chimpanzee/orangutan increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="sGsG")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="sGsG")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="sGsG")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="sGsG")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Chimpanzee/Orangutan Increases\n(657 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_sGsG <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Chimpanzee/Orangutan Increases\n(657 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 21 of 31
# Chimpanzee/orangutan decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for chimpanzee/orangutan decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="sLsL")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="sLsL")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="sLsL")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="sLsL")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Chimpanzee/Orangutan Decreases\n(252 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_sLsL <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Chimpanzee/Orangutan Decreases\n(252 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 22 of 31
# Gorilla/orangutan increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for gorilla/orangutan increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="ssGG")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="ssGG")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="ssGG")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="ssGG")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Gorilla/Orangutan Increases\n(2,120 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_ssGG <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Gorilla/Orangutan Increases\n(2,120 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 23 of 31
# Gorilla/orangutan decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for gorilla/orangutan decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="ssLL")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="ssLL")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="ssLL")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="ssLL")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Gorilla/Orangutan Decreases\n(517 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_ssLL <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Gorilla/Orangutan Decreases\n(517 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 24 of 31
# Human/chimpanzee/gorilla increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/chimpanzee/gorilla increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="GGGs")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="GGGs")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="GGGs")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="GGGs")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Chimpanzee/Gorilla Increases\n(1,736 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_GGGs <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Chimpanzee/Gorilla\nIncreases (1,736 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 25 of 31
# Human/chimpanzee/gorilla decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/chimpanzee/gorilla decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="LLLs")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="LLLs")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="LLLs")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="LLLs")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Chimpanzee/Gorilla Decreases\n(1,453 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_LLLs <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Chimpanzee/Gorilla\nDecreases (1,453 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 26 of 31
# Human/chimpanzee/orangutan increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/chimpanzee/orangutan increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="GGsG")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="GGsG")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="GGsG")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="GGsG")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Chimpanzee/Orangutan\nIncreases (794 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_GGsG <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Chimpanzee/Orangutan\nIncreases (794 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 27 of 31
# Human/chimpanzee/orangutan decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/chimpanzee/orangutan decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="LLsL")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="LLsL")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="LLsL")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="LLsL")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Chimpanzee/Orangutan\nDecreases (906 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_LLsL <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Chimpanzee/Orangutan\nDecreases (906 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 28 of 31
# Human/gorilla/orangutan increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/gorilla/orangutan increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="GsGG")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="GsGG")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="GsGG")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="GsGG")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Gorilla/Orangutan\nIncreases (691 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_GsGG <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Gorilla/Orangutan\nIncreases (691 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 29 of 31
# Human/gorilla/orangutan decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for human/gorilla/orangutan decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="LsLL")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="LsLL")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="LsLL")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="LsLL")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Human/Gorilla/Orangutan\nDecreases (554 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_LsLL <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Human/Gorilla/Orangutan\nDecreases (554 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 30 of 31
# Chimpanzee/gorilla/orangutan increases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for chimpanzee/gorilla/orangutan increases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="sGGG")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="sGGG")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="sGGG")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="sGGG")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Chimpanzee/Gorilla/Orangutan\nIncreases (939 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_sGGG <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Chimpanzee/Gorilla/Orangutan\nIncreases (939 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot 31 of 31
# Chimpanzee/gorilla/orangutan decreases
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Generating beta value density plot for chimpanzee/gorilla/orangutan decreases ... ")

bh <- as.data.frame(subset(glm_results, flag_change_type=="sLLL")[,22])
bc <- as.data.frame(subset(glm_results, flag_change_type=="sLLL")[,25])
bg <- as.data.frame(subset(glm_results, flag_change_type=="sLLL")[,28])
bo <- as.data.frame(subset(glm_results, flag_change_type=="sLLL")[,31])

bh$Species <- "Human"
bc$Species <- "Chimpanzee"
bg$Species <- "Gorilla"
bo$Species <- "Orangutan"

colnames(bh) <- c("beta", "Species")
colnames(bc) <- c("beta", "Species")
colnames(bg) <- c("beta", "Species")
colnames(bo) <- c("beta", "Species")

beta_values <- rbind(bh,bc,bg,bo)
beta_values$Species <- factor(beta_values$Species, levels=c("Human","Chimpanzee","Gorilla","Orangutan"))

density_plot <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_pdf_line_colors) + scale_size_manual(values=my_pdf_line_widths) + scale_linetype_manual(values=my_pdf_line_types) +
		ggtitle("Chimpanzee/Gorilla/Orangutan\nDecreases (686 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif"), plot.title=element_text(hjust=0.5, size=18, face="bold"),
		axis.text=element_text(size=14), axis.title=element_text(size=14),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12))

density_plot_sLLL <- ggplot(beta_values, aes(x=beta, color=Species, size=Species, linetype=Species)) + geom_line(stat="density") + xlab("Beta Value") + ylab("Density") + xlim(-5,5) +
		scale_color_manual(values=my_svg_line_colors) + scale_size_manual(values=my_svg_line_widths) + scale_linetype_manual(values=my_svg_line_types) +
		ggtitle("Chimpanzee/Gorilla/Orangutan\nDecreases (686 DHS Sites)") +
		theme_bw() + theme(text=element_text(family="serif", face="bold"), plot.title=element_text(hjust=0.5, size=my_svg_text_size),
		axis.text=element_text(size=my_svg_text_size, colour="black"), axis.title=element_text(size=my_svg_text_size),
		axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)), axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
		panel.grid.major=element_line(size=0.5, linetype='solid', colour="grey"), panel.grid.minor=element_line(size=0.25, linetype='solid', colour="grey"),
		legend.position="none")

print(density_plot)

dev.off()

cat("done.\n")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Close the pdf file and save the plots to png files
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Saving plots to svg files ...")

svg(file="beta_value_density_plot.ssss.svg", bg="transparent", width=10, height=10)
print(density_plot_ssss)
dev.off()

svg(file="beta_value_density_plot.Gsss.svg", bg="transparent", width=10, height=10)
print(density_plot_Gsss)
dev.off()

svg(file="beta_value_density_plot.Lsss.svg", bg="transparent", width=10, height=10)
print(density_plot_Lsss)
dev.off()

svg(file="beta_value_density_plot.sGss.svg", bg="transparent", width=10, height=10)
print(density_plot_sGss)
dev.off()

svg(file="beta_value_density_plot.sLss.svg", bg="transparent", width=10, height=10)
print(density_plot_sLss)
dev.off()

svg(file="beta_value_density_plot.ssGs.svg", bg="transparent", width=10, height=10)
print(density_plot_ssGs)
dev.off()

svg(file="beta_value_density_plot.ssLs.svg", bg="transparent", width=10, height=10)
print(density_plot_ssLs)
dev.off()

svg(file="beta_value_density_plot.sssG.svg", bg="transparent", width=10, height=10)
print(density_plot_sssG)
dev.off()

svg(file="beta_value_density_plot.sssL.svg", bg="transparent", width=10, height=10)
print(density_plot_sssL)
dev.off()

svg(file="beta_value_density_plot.GGss.svg", bg="transparent", width=10, height=10)
print(density_plot_GGss)
dev.off()

svg(file="beta_value_density_plot.LLss.svg", bg="transparent", width=10, height=10)
print(density_plot_LLss)
dev.off()

svg(file="beta_value_density_plot.GsGs.svg", bg="transparent", width=10, height=10)
print(density_plot_GsGs)
dev.off()

svg(file="beta_value_density_plot.LsLs.svg", bg="transparent", width=10, height=10)
print(density_plot_LsLs)
dev.off()

svg(file="beta_value_density_plot.GssG.svg", bg="transparent", width=10, height=10)
print(density_plot_GssG)
dev.off()

svg(file="beta_value_density_plot.LssL.svg", bg="transparent", width=10, height=10)
print(density_plot_LssL)
dev.off()

svg(file="beta_value_density_plot.sGGs.svg", bg="transparent", width=10, height=10)
print(density_plot_sGGs)
dev.off()

svg(file="beta_value_density_plot.sLLs.svg", bg="transparent", width=10, height=10)
print(density_plot_sLLs)
dev.off()

svg(file="beta_value_density_plot.sGsG.svg", bg="transparent", width=10, height=10)
print(density_plot_sGsG)
dev.off()

svg(file="beta_value_density_plot.sLsL.svg", bg="transparent", width=10, height=10)
print(density_plot_sLsL)
dev.off()

svg(file="beta_value_density_plot.ssGG.svg", bg="transparent", width=10, height=10)
print(density_plot_ssGG)
dev.off()

svg(file="beta_value_density_plot.ssLL.svg", bg="transparent", width=10, height=10)
print(density_plot_ssLL)
dev.off()

svg(file="beta_value_density_plot.GGGs.svg", bg="transparent", width=10, height=10)
print(density_plot_GGGs)
dev.off()

svg(file="beta_value_density_plot.LLLs.svg", bg="transparent", width=10, height=10)
print(density_plot_LLLs)
dev.off()

svg(file="beta_value_density_plot.GGsG.svg", bg="transparent", width=10, height=10)
print(density_plot_GGsG)
dev.off()

svg(file="beta_value_density_plot.LLsL.svg", bg="transparent", width=10, height=10)
print(density_plot_LLsL)
dev.off()

svg(file="beta_value_density_plot.GsGG.svg", bg="transparent", width=10, height=10)
print(density_plot_GsGG)
dev.off()

svg(file="beta_value_density_plot.LsLL.svg", bg="transparent", width=10, height=10)
print(density_plot_LsLL)
dev.off()

svg(file="beta_value_density_plot.sGGG.svg", bg="transparent", width=10, height=10)
print(density_plot_sGGG)
dev.off()

svg(file="beta_value_density_plot.sLLL.svg", bg="transparent", width=10, height=10)
print(density_plot_sLLL)
dev.off()

svg(file="beta_value_density_plot.GGGG.svg", bg="transparent", width=10, height=10)
print(density_plot_GGGG)
dev.off()

svg(file="beta_value_density_plot.LLLL.svg", bg="transparent", width=10, height=10)
print(density_plot_LLLL)
dev.off()

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Finish up plots
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cat("done.\n\nPlots are located in beta_value_density_plot.all.pdf and beta_value_density_plots.XXXX.svg where XXXX is the changetype flag\n\n")

# ==========================================================================================================================================================================
# Finish up script
# ==========================================================================================================================================================================

cat("---------------------------------------------------------------------------------------\n")
cat("END:", date(), "\n")
cat("---------------------------------------------------------------------------------------\n")

sink()
