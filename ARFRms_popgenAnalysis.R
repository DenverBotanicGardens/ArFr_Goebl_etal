#' ---
#' title: "ARFR_popgen.R"
#' author: "Alyson Emery"
#' output: html_document
#' ---
#+ setup, include=FALSE
knitr::opts_chunk$set(collapse = TRUE)
#' ## Packages
#' 
#' These packages are used multiple times throughout the scripts, and should be loaded prior to performing analyses. Packages used in only one script are loaded at the beginning of each section.
# load packages
library(ggplot2)
library(dplyr)
library(cowplot)

#' ## Load Data and Color Palette
#' 
#' We use the same color palette for populations throughout our R code. Here we assign hex codes for our populations and establish the order that they show in in legends throughout our plots.
az1col <- "#9E0142"
az2col <- "#D53E4F"
nm1col <- "#F46D43"
co1col <- "#FDAE61"
co2col <- "#FEE08B"
co3col <- "#FFFF85"
ut1col <- "#aabe7f"
co4col <- "#ABDDA4"
wy1col <- "#66C2A5"
wy2col <- "#3288BD"
wy3col <- "#5E4FA2"

# create a manual vector for colors for ibd plot
pop_colors <- c("AZ1" = "#9E0142",
                "AZ2" = "#D53E4F",
                "NM1" = "#F46D43",
                "CO1" = "#FDAE61",
                "CO2" = "#FEE08B",
                "CO3" = "#FFFF85",
                "UT1" = "#aabe7f",
                "CO4" = "#ABDDA4",
                "WY1" = "#66C2A5",
                "WY2" = "#3288BD",
                "WY3" = "#5E4FA2")

# set up factor levels
pop_order <- c("WY3", "WY2", "WY1", "CO4", "UT1", "CO3", "CO2", "CO1", "NM1", "AZ2", "AZ1")

# order pop colors
pop_colors <- pop_colors[pop_order]
#' Additionally, we have metadata for both our individual samples and populations to load into our R environment.
#' 
# read sample data file
sample_data <- read.csv("data/ARFRv5_strata.csv", stringsAsFactors = FALSE)

#' ## Data Filtering
#' 
#' ### Genotype (012) Matrix Assessment
#' 
#' Here we looked at heterozygosity in our dataset using a 012 matrix created in vcftools.
# set working directory and read in 012 file.
#setwd(dir = "./ARFR")
AGartemisia.012 <- read.delim("./data/AGartemisia.012", header=FALSE, row.names=1)

# create function to count the number of 1s, 2s, 0s, and -1s in the 012 file
count_values <- function(matrix) {
  result <- apply(matrix[, drop = FALSE], 2, table)
  return(result)
}

all_counts <- count_values(AGartemisia.012)

# check that the function worked. tail should reflect the number of variants in dataset, and adding the counts for each group should equal the number of individuals in the dataset.
tail(all_counts, n = 3)

# create formula to count only the values specified
count_specific_values <- function(matrix, value) {
  result <- sapply(matrix[, drop = FALSE], function(col) sum(col == value))
  return(result)
}

# count heterozygous individuals at each SNP
counts_1s <- count_specific_values(AGartemisia.012, value = 1)
head(counts_1s)

# count individuals with missing data at each SNP
counts_missing <- count_specific_values(AGartemisia.012, value = -1)
head(counts_missing)

# get the ratio of heterozygous counts to total possible counts (i.e. only from individuals with data)
# this is looking at the percentage of individuals in which that site is heterozygous.
het_ratio <- counts_1s / (96 - counts_missing)
head(het_ratio)

hetratio_df <- data.frame(Variant = colnames(AGartemisia.012), Percent_Heterozygous = het_ratio)
head(hetratio_df)

# print a histogram of heterozygous sites
het_histogram <- ggplot(data = hetratio_df, aes(x = Percent_Heterozygous)) +
  geom_histogram(fill = "blue", color = "black") +
  xlab("Percent Heterozygous") +
  ylab("Variants")

ggsave(
  filename = "Heterozygous_SNPs_Histogram.pdf",
  plot = het_histogram,
  device = "pdf",
  path = "./results/filtering/",
  scale = 1,
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 300
)
print(het_histogram)

#' We saw that there were around 100 variants that were unusually heterozygous. This led us to consider filtering for paralogs.  
#'   
#' We also used the 012 matrix to assess how well sequencing performed by comparing a sequenced sample to its replicate.
#' 
# create a dataframe for just the individual and its replicate
selected_rows <- c(1, 96)
ARFR25andRep_df <- as.data.frame(AGartemisia.012[selected_rows, ])

# change row names for readability
ARFR25andRep_df$new_row_name <- paste("ind_" , row.names(ARFR25andRep_df))
rownames(ARFR25andRep_df) <- ARFR25andRep_df$new_row_name

# transpose data and write csv
ARFR25andRep <- as.data.frame(t(ARFR25andRep_df))
head(ARFR25andRep)
write.csv(ARFR25andRep, file = "./results/replicate_check_Ind25.csv")

# compare the two columns to see how much they match
valid_rows <- ARFR25andRep$`ind_ 0` != -1 & ARFR25andRep$`ind_ 95` != -1

match_percent <- mean(ARFR25andRep$`ind_ 0`[valid_rows] == ARFR25andRep$`ind_ 95`[valid_rows]) * 100
print(match_percent)

#' Of the sites that had a sequencing read in both, 98.19% matched, indicating low sequencing error rates.  
#' 
#' We then compared the amount of missingness between our sample and its replicate to decide which to maintain in the dataset and which to remove.

# compare missingness between sample and replicate
count_missing_ARFR25 <- table(ARFR25andRep$`ind_ 0`)
print(count_missing_ARFR25["-1"])

count_missing_replicate <- table(ARFR25andRep$`ind_ 95`)
print(count_missing_replicate["-1"])
#' Because the replicate had fewer missing sites than the initial sample, we decided to keep the replicate in the dataset and remove the original.
#' 
#' ### Paralog Identification
#' 
#' To identify paralogs, we first assessed which SNPs were out of Hardy-Weinberg equilibrium. We then compared that with the mean depth at each site to decide what filter would be effective at removing paralogs but unlikely to remove non-paralogs.
#' 

## creating the dataframe

# read in mean depth file from vcftools
mean_depth <- read.csv(file = "./data/AGartemisia.ldepth.mean.csv", sep = "\t", header = TRUE)

# combine "CHROM" and "POS" columns to make a single column with a unique identifier for each variant
mean_depth$CHROM_POS <- paste(mean_depth$CHROM,mean_depth$POS)
head(mean_depth)

# read in HWE output from vcftools
hwe_output <- read.csv(file = "./data/HWE_output.csv", header = TRUE)

# combine "CHROM" and "POS" columns to make a single column with a unique identifier for each variant
hwe_output$CHROM_POS <- paste(hwe_output$CHROM,hwe_output$POS)
head(hwe_output)

# merge mean depth and HWE outputs by "CHROM_POS"
master_depth <- full_join(mean_depth,hwe_output, by = "CHROM_POS")

#' Once we had a data frame, we plotted the p-value for heterozygozity excess by mean depth at a few different scales to get a sense of where we should make a cutoff.

# plot mean depth by heterozygosity excess, no filter
plot_all_depth <- ggplot(master_depth, aes(x = MEAN_DEPTH, y=P_HET_EXCESS)) +
  geom_point()
print(plot_all_depth)
# based on this, test filter of mean depth less than 1000
criteria <- 1000
test_filter <-filter(master_depth, MEAN_DEPTH < criteria)

# plot
plot_testfilter_depth <- ggplot(test_filter, aes(x = MEAN_DEPTH, y=P_HET_EXCESS)) +
  geom_point()
print(plot_testfilter_depth)
# based on the plot above, a mean depth below approximately 500 bp is where heterozygosity excess is evenly distributed
criteria2 <- 500
master_depth_filter <- filter(master_depth, MEAN_DEPTH < criteria2)

# to see which sites were removed
master_depth_filter_junk <- filter(master_depth, MEAN_DEPTH > criteria2)

# plot
plot_filter_depth <- ggplot(master_depth_filter, aes(x = MEAN_DEPTH, y=P_HET_EXCESS)) +
  geom_point()
print(plot_filter_depth)
# plot all together and save
all_plots <- plot_grid(plot_all_depth,plot_testfilter_depth,plot_filter_depth, labels = c("A","B","C"))

ggsave2(
  filename = "Depth_Filter.pdf",
  plot = all_plots,
  device = "pdf",
  path = "./results/filtering/",
  scale = 1,
  width = 7.75,
  height = 6,
  units = c("in"),
  dpi = 300
)
print(all_plots)

#' We decided to filter out variants that had both A) a low p-value for HWE, and B) a depth higher than where variants had a evenly distributed HWE. Based on our data, we decided to filter by a mean depth of 500 bp.
#' 
#' ### Missingness in Individuals
#' 
#' Our dataset was prefiltered by SNPsaurus to only include variants with 25% missingess or less. To assess missingness in individuals, we ran --missing-ind in vcftools. We ran this after filtering for mean site depth, though not yet having removed our replicate sample from the dataset.
#' 
# read in file
miss_ind <- read.delim("./data/out.imiss")
head(miss_ind)

# plot percentage missing data for individuals
miss_ind_plot <- ggplot(data = miss_ind, aes(x = F_MISS)) +
  geom_histogram(fill = "darkgreen", color = "black") +
  ylab("Individuals")
print(miss_ind_plot)
# plot where the majority of the data is
miss_ind_plot2 <- ggplot(data = miss_ind, aes(x = F_MISS)) +
  geom_histogram(fill = "darkgreen", color = "black") +
  ylab("Individuals") +
  xlim(0.05,0.3)
print(miss_ind_plot2)

#' We recognized that there were 8 individuals with missingness higher than most of the others in our dataset. We decided to filter out individuals that had greater than 20% missing data after implementing our depth filter.

# get samples ids for individuals with >20% missingness
criteria3 <- 0.20
missind_filter_junk <- filter(miss_ind, miss_ind$F_MISS > criteria3)
ids_to_remove <- missind_filter_junk$INDV
print(ids_to_remove)

# add the sample id for the replicate individual to list
ids_to_remove <- append(ids_to_remove,"ARFR_25_sorted")
print(ids_to_remove)

# write list into a .txt file for removal in vcftools
write.table(ids_to_remove, file = "./results/filtering/ids_to_remove.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#' ## Analysis of Genetic Diversity
#' 
#' ### Nucleotide Diversity
#' 
#' Here, we calculated nucleotide diversity by population with .pi output from vcftools.

# list population-specific pi files
pilist <- c("data/v5.pi.sites/WY3.sites.pi", "data/v5.pi.sites/WY2.sites.pi", "data/v5.pi.sites/WY1.sites.pi", "data/v5.pi.sites/CO4.sites.pi", "data/v5.pi.sites/UT1.sites.pi", "data/v5.pi.sites/CO3.sites.pi", "data/v5.pi.sites/CO2.sites.pi", "data/v5.pi.sites/CO1.sites.pi", "data/v5.pi.sites/NM1.sites.pi", "data/v5.pi.sites/AZ2.sites.pi", "data/v5.pi.sites/AZ1.sites.pi")

# initialize an empty vector to store the means
sp_popmeans <- numeric(length(pilist))
sp_popvariance <- numeric(length(pilist))
sp_popsdev <- numeric(length(pilist))

# loop through each file
for (i in 1:length(pilist)) {
  # read the file
  sitepi <- read.delim(pilist[i])
  
  # calculate the mean of the desired column (replace "column_name" with the actual column name)
  sp_popmeans[i] <- mean(sitepi$PI, na.rm = TRUE)
  sp_popvariance[i] <- var(sitepi$PI, na.rm = TRUE)
  sp_popsdev[i] <- sd(sitepi$PI, na.rm = TRUE)
  
  # print the mean for this file
  cat("Mean of", pilist[i], ":", sp_popmeans[i], "\n")
}

print(sp_popmeans)
print(sp_popvariance)
print(sp_popsdev)

# make a data frame with these stats
ARFR_sitepi.df <- data.frame(Population = pop_order, Site_Pi_Mean = sp_popmeans, Site_Pi_Variance = sp_popvariance, Site_Pi_SD = sp_popsdev)

# print the dataframe
print(ARFR_sitepi.df)
#' These results are visualized alongside the inbreeding coefficient F in the section below.
#' 
#' ### Inbreeding Coefficient F
#' 
#' Calculating the inbreeding coefficient F using output from vcftools "het".
#' 
# Read your .het output file
AGartemisia.het <- read.delim("data/AGartemisia.v5.het")
head(AGartemisia.het)

# add pop codes to het file from sample data
AGartemisia.het$Population <- sample_data$pop_abbr

# check to make sure they align correctly
head(AGartemisia.het)

# list table of means of F by population
ARFR_F_popmeans <- aggregate(AGartemisia.het$F, list(AGartemisia.het$Population), FUN=mean)

# list table of standard deviations of F by population
ARFR_F_popsdev <- aggregate(AGartemisia.het$F, list(AGartemisia.het$Population), FUN=sd)

# lists table of standard deviations of F by population
ARFR_F_popvariance <- aggregate(AGartemisia.het$F, list(AGartemisia.het$Population), FUN=var)

# merge these into a single data frame
ARFR_F_df <- merge(merge(ARFR_F_popmeans, ARFR_F_popvariance, by = "Group.1"), ARFR_F_popsdev, by = "Group.1")
names(ARFR_F_df) <- c("Population", "F_Mean", "F_Variance", "F_SD") # List new column names in the same order as existing columns
head(ARFR_F_df)
#' ### Plotting Genetic Diversity
#' We merged our inbreeding coefficient dataframe with our nucleotide diversity dataframe to plot them together.

# merge by population in both dataframes
ARFR_df <- merge(ARFR_F_df, ARFR_sitepi.df, by = "Population")
head(ARFR_df)

# export table to csv
write.csv(ARFR_df, file = "results/NeutralDiversityStats.csv", row.names = FALSE)

# plot F by nucleotide diversity
ARFR_df$Population <- factor(ARFR_df$Population, levels = pop_order)

FbyPI <- ggplot(data = ARFR_df, aes(Site_Pi_Mean, F_Mean)) +
  geom_pointrange(aes(ymin = F_Mean - F_Variance, ymax = F_Mean + F_Variance, color = Population)) +
  geom_point(shape = 1, colour = "black", size = 3) +
  scale_color_manual(values = c(wy3col, wy2col, wy1col, co4col, ut1col, co3col, co2col, co1col, nm1col, az2col, az1col)) +
  xlab("π Mean") +
  ylab("F Mean") +
  theme_bw() +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.95))

# save and print plot
ggsave2(
  filename = "FbyPI.pdf",
  plot = FbyPI,
  device = "pdf",
  path = "results/analysisfigures/",
  scale = 1,
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 300
)

print(FbyPI)
#' 
#' ## Population Differentiation
#' ### F~ST~ and Mantel Test for Isolation by Distance
#' 
#' Using dartR, we calculate Weir and Cockerham's Pairwise F~ST~ while simultaneously evaluating isolation by distance with a Mantel test. We start by creating a genlight object to use in the program and running the model.
#' 
#' 
# load packages
suppressPackageStartupMessages(library(dartR))
library(vcfR)
library(reshape2)

# read vcf file and convert to genlight
vcf <- read.vcfR("data/AGartemisia_v5.recode.vcf")
gl <- vcfR2genlight(vcf)

# add pop info to df
pop(gl) <- as.factor(sample_data$pop_abbr)

# ensure order of coordinates matches order of individuals
sample_data <- sample_data[match(indNames(gl), sample_data$vcf_id), ]

# attach coordinates to the genlight object
gl@other$latlon <- sample_data[, c("lat", "lon")]

# running ibd with populations
ibd_result <- gl.ibd(
  x = gl,
  distance = "Fst",
  coordinates = "latlon",
  paircols = "pop",
  plot.out = FALSE
)

#' We extracted F~ST~ from our data list to create a heat map. The plot was customized so that the gray scale would reflect the range of our F~ST~ data.

# convert Dgen and Dgeo to matrices
Dgen_mat <- as.matrix(ibd_result$Dgen)
Dgeo_mat <- as.matrix(ibd_result$Dgeo)
pops <- rownames(Dgen_mat)
n <- length(pops)

# view and save raw fst matrix
print(Dgen_mat)
write.csv(Dgen_mat, file = "./data/fstmatrix.csv")

# make correlation data frame
Fstcorrplot_df <- melt(Dgen_mat)
colnames(Fstcorrplot_df) <- c("X", "Y", "Fst")
Fstcorrplot_df$X <- factor(Fstcorrplot_df$X, levels = c("WY3", "WY2", "WY1", "CO4", "UT1", "CO3", "CO2", "CO1", "NM1", "AZ1", "AZ2")) # levels are for population labels
Fstcorrplot_df$Y <- factor(Fstcorrplot_df$Y, levels = c("AZ2","AZ1","NM1","CO1","CO2","CO3","UT1","CO4","WY1","WY2","WY3")) # levels are for populations labels

# plot fst heatmap using ggplot2
fst_plot <- ggplot(Fstcorrplot_df, aes(x = X, y = Y, fill = Fst)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "white", high = "#151414", 
                       limit = c(0,0.081), 
                       name="Fst",
                       breaks=seq(-1, 1, by=.02)) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  coord_fixed() +
  labs(x = NULL, y = NULL)

# save and print
ggsave(
  filename = "fst.pdf",
  plot = fst_plot,
  device = "pdf",
  path = "./results/analysisfigures/",
  scale = 1,
  width = 7,
  height = 5,
  units = c("in"),
  dpi = 300
)
print(fst_plot)

#' Next, we created a data frame to extract the statistics from the mantel test and linear model, and to create a plot for isolation by distance.
# get all off-diagonal combinations (pop-pop)
pair_idx <- which(row(Dgen_mat) != col(Dgen_mat), arr.ind = TRUE)

# make isolation by distance dataframe
ibd_df <- data.frame(
  pop1 = pops[pair_idx[, 1]],
  pop2 = pops[pair_idx[, 2]],
  dist_gen = Dgen_mat[pair_idx],
  dist_geo = Dgeo_mat[pair_idx]
)

# apply this order to both pop columns
ibd_df$pop1 <- factor(ibd_df$pop1, levels = pop_order)
ibd_df$pop2 <- factor(ibd_df$pop2, levels = pop_order)

# finding regression formula
model <- lm(dist_gen ~ dist_geo, data = ibd_df)
slope <- coef(model)[2]
intercept <- coef(model)[1]

# format outputs
slope_formatted <- formatC(slope, format = "e", digits = 2)  # e.g., 1.23e-05
intercept_formatted <- round(intercept, 3)

# calculate r squared for linear model
r2 <- summary(model)$r.squared
print(r2)

# significance of the linear model
p_val <- summary(model)$coefficients[2, 4]
print(p_val)

#' At this point, we can see that our R squared value for the linear regression model is quite low - explaining less than 10% of the variation in our data.

# extract mantel correlation coefficient
mantel_r <- ibd_result$mantel$statistic

# compute mantel r²
mantel_r2 <- mantel_r^2

# extract mantel test p-value
mantel_p <- ibd_result$mantel$signif

# print mantel statistics
print(mantel_r)
print(mantel_r2)
print(mantel_p)

# format a label for the plot
label_text <- paste0("y = ", slope_formatted, "x + ", intercept_formatted,
                     "\nR² = ", round(r2, 4), 
                     "\np = ", format.pval(mantel_p, digits = 3, eps = .001))

# plot and save
ibd_plot <- ggplot(ibd_df, aes(x = dist_geo, y = dist_gen)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray70", linewidth = 1) +
  # Outline for large points
  geom_point(aes(color = pop1), size = 4.6, stroke = 0, shape = 16, alpha = 1, color = "black") +
  geom_point(aes(color = pop1), size = 4, shape = 16, alpha = 1, show.legend = TRUE) +
  # Outline for small points
  geom_point(aes(color = pop2), size = 2.6, stroke = 0, shape = 16, alpha = 1, color = "black") +
  geom_point(aes(color = pop2), size = 2, shape = 16, alpha = 1, show.legend = FALSE) +
  scale_color_manual(values = pop_colors, name = "Population") +
  labs(
    x = "Geographic Distance (m)",
    y = "Genetic Distance (Fst)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  annotate("text", 
           x = Inf, y = Inf, 
           label = label_text, 
           hjust = 1, vjust = 1.1, 
           size = 4, fontface = "italic")

# save and print to screen
ggsave(
  filename = "ibd.pdf",
  plot = ibd_plot,
  device = "pdf",
  path = "./results/analysisfigures/",
  scale = 1,
  width = 8,
  height = 5,
  units = c("in"),
  dpi = 300
)
print(ibd_plot)
#' 
#' ## Population Structure
#' ### PCA in SNPrelate
#' Here we used SNP data to perform a principal component analysis using our vcf file and the R package SNPRelate.
#' 
# load packages
library(gdsfmt)
library(SNPRelate)

# make .gds file out of the vcf (saves file to directory)
snpgdsVCF2GDS("./data/AGartemisia_v5.recode.vcf", "./data/ARFR.gds", method="biallelic.only")

# summarize
snpgdsSummary("./data/ARFR.gds")

# open the gds file
genofile <- snpgdsOpen(file = "./data/ARFR.gds")

# prune snps for ld
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only = FALSE, maf = NaN, missing.rate = NaN, verbose = FALSE)
str(snpset, list.len = 5)

# get all selected snp id
snpset.id <- unlist(unname(snpset))
head(snpset.id)

# run PCA on genofile with snpset.id you created
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = FALSE)

# calculate variance explained by each component
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
head(sample.id)

# ensure the order of sample ids is as the same as order in sample_data for adding population codes
head(cbind(sample.id, sample_data$pop_abbr))

# make a data.frame
pca.df <- data.frame(sample.id = pca$sample.id,
                        Population = factor(sample_data$pop_abbr)[match(pca$sample.id, sample.id)],
                        EV1 = pca$eigenvect[,1],    # the first eigenvector
                        EV2 = pca$eigenvect[,2],    # the second eigenvector
                        EV3 = pca$eigenvect[,3],
                        EV4 = pca$eigenvect[,4],
                        stringsAsFactors = FALSE)
head(pca.df)

# save dataframe
write.csv(pca.df, file = "./results/pca_table.csv")

# order populations within dataframe
pca.df$Population <- factor(pca.df$Population, levels = pop_order)

# create plot
pc1pc2_plot <- ggplot(data = pca.df, mapping = aes(x = EV1, y = EV2)) +
  geom_point(aes(color=Population)) +
  scale_y_reverse() +
  coord_fixed() +
  theme_bw() +
  ylab("Principal Component 2 (2.20% variance)") +
  xlab("Principal Component 1 (2.74% variance)") +
  scale_color_manual(values = c(wy3col, wy2col, wy1col, co4col, ut1col, co3col, co2col, co1col, nm1col, az2col, az1col)) +
  guides(fill = guide_legend(title = "Population"))

# save and print to screen
ggsave(
  filename = "PC1PC2.pdf",
  plot = pc1pc2_plot,
  device = "pdf",
  path = "./results/analysisfigures/",
  scale = 1,
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 300
)
print(pc1pc2_plot)
# plot PC 3 and PC 4
pc3pc4_plot <- ggplot(data = pca.df, mapping = aes(x = EV3, y = EV4)) +
  geom_point(aes(color=Population)) +
  coord_fixed() +
  theme_bw() +
  ylab("Principal Component 4 (1.92% variance)") +
  xlab("Principal Component 3 (1.97% variance)") +
  scale_color_manual(values = c(wy3col, wy2col, wy1col, co4col, ut1col, co3col, co2col, co1col, nm1col, az2col, az1col)) +
  guides(fill = guide_legend(title = "Population"))

# save and print to screen
ggsave(
  filename = "PC3PC4.pdf",
  plot = pc3pc4_plot,
  device = "pdf",
  path = "./results/analysisfigures/",
  scale = 1,
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 300
)
print(pc3pc4_plot)

#' ### STRUCTURE in pophelper
#' 
#' We ran an admixture analysis using the program STRUCTURE, and evaluate and visualize the results here. We evaluate the results using the Evanno method to test for the most likely number of groups (K) when comparing across runs. For both the Evanno method and plot visualization, we use the package "pophelper", which can be obtained via "install_github('royfrancis/pophelper')".
#' 
# load packages
library(pophelper)
library(gridExtra)

# load in structure results files
sfiles <- list.files(path = "./data/structureresults", full.names = T)
slist <- readQ(files = sfiles, filetype = "structure", indlabfromfile=TRUE)
head(slist[[1]])

# extract individual order from files
indorder <- slist[[1]]
write.csv(indorder, file = "./data/structureindorder.csv")

# assess k values using evanno method
sr1 <- summariseQ(tabulateQ(slist))
write.csv(evannoMethodStructure(sr1), "./results/EvannoMethodStructure.csv", na = "")
evanno_plot <- evannoMethodStructure(data = sr1, exportplot = F, returnplot = T, returndata = F, basesize = 12, linesize = 0.7, xaxisbreaks = 1:11)
# arrange, save, and print
emp <- arrangeGrob(evanno_plot)
ggsave("./results/analysisfigures/evanno_plot.pdf", emp, width = 8, height = 8)
grid.arrange(evanno_plot)

#' The Evanno method supported K=1 as the most likely number of groups, with Ks 2-4 as the next most likely. Once we had decided on the most supported K values with the Evanno method, we merged and plotted results for these STRUCTURE runs.
#' 
# align k values in the file list
slist <- alignK(slist)
length(slist)
#> [1] 17

# merge files into lists for analysis and plotting
K2_slist <- mergeQ(slist[c(6,7,8,9,10)])
K3_slist <- mergeQ(slist[c(11,12,13,14,15)])
K4_slist <- mergeQ(slist[c(16,17,18,19,20)])
K5_slist <- mergeQ(slist[c(21,22,23,24,25)])
c_slist <- c(K2_slist, K3_slist, K4_slist, K5_slist)

# add pop information to dataframe
prefixes <- substr(rownames(indorder), 1, 11) # assuming indorder row names are truncated
indorder$pop_abbr <- sample_data$pop_abbr[match(prefixes, substr(sample_data$vcf_id, 1, 11))]
indorder$id <- sample_data$id[match(prefixes, substr(sample_data$vcf_id, 1, 11))]
head(indorder)

# create data frame of only pop labels in correct order
pop_lab <- indorder[,11,drop=FALSE]

# sort the groups by latitude using subsetgrp and plot
structure_plot <- plotQ(alignK(c_slist[c(1,2,3)]), imgoutput = "join", clustercol = c("#B0ADAB", "#E1E0E0", "#787675", "#333130"), returnplot = T, exportplot = F, basesize = 11, grplab = pop_lab, subsetgrp = pop_order, showindlab = F, useindlab = T, grplabsize = 3, splab = c("K=2","K=3", "K=4"), indlabsize = 5, showlegend = FALSE)
# arrange, save, and print
sp <- arrangeGrob(
  structure_plot$plot[[1]])
ggsave("./results/analysisfigures/structure_plot.pdf", sp, width = 8, height = 5)
grid.arrange(structure_plot$plot[[1]])