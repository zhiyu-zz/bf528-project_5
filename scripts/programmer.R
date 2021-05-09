# set working directory
setwd("/projectnb2/bf528/users/dreadlocks/project_5/sarah/")
outdir <- "/projectnb2/bf528/users/dreadlocks/project_5/sarah/DESeq2_res/"

library(tidyverse)
library(stringr)
library(ggplot2)
library(DESeq2)
library(apeglm)
library(grid)
library(gridExtra)
library(EnhancedVolcano)
# Combine counts files from each sample into a single comma-delimited text file

# list files in the "feature" folder with filenames that end with ".txt"
file_names <- list.files("feature", pattern = ".txt$", full.names = T)
sample_names <- str_extract(file_names, "SRR\\d*") # extract sample names 

# get gene id column
gene_id <- read.table(file_names[1], header=T)
# initialize count matrix with gene id
count_matrix <- gene_id["Geneid"]


# read txt files and combine count from each sample to new column
for (i in 1:length(file_names)){
  sp <- read.table(file_names[i], header = T)
  sp <- sp[,7] # select last column
  count_matrix <- bind_cols(count_matrix, sp)
  colnames(count_matrix)[i+1] <- sample_names[i] # assign sample id to each count column. i+1 because first column is gene id
}

mcm <- melt(count_matrix)
head(count_matrix)


boxplot(count_matrix[,-1], pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, xlab="Counts", outline = F,
        col="blue")
# histogram of counts per gene (total)
par(mfrow=c(2,1))
hist(as.matrix(count_matrix[,-1]), col="blue", border="white",
     breaks=30000, xlim=c(0,1500), main="Counts per gene",
     xlab="Counts", ylab="Number of genes", 
     las=1, cex.axis=0.7)
epsilon <- 1 # pseudo-count to avoid problems with log(0)
hist(as.matrix(log2(count_matrix[,-1] + epsilon)), breaks=100, col="blue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

#-----------------------------DESeq2-------------------------------------------#
# metadata
samples <- read.csv("/projectnb/bf528/users/dreadlocks/project_3/curator/sample/toxgroup_6_rna_info.csv")
# count matrix for control groups
control_matrix <- read.csv("/project/bf528/project_3/samples/control_counts.csv")

### Data preprocessing 

# check if gene ids of treatment and control are in the same order 
all(count_matrix$Geneid == control_matrix$Geneid)
# subset control count matrix
control_matrix <- control_matrix[, samples$Run[samples$mode_of_action == "Control"]]
# combine treatment and control count matrices 
count_matrix <- cbind(count_matrix, control_matrix) # 9 treatment samples, 6 control samples
# remove genes which count equals zero across all samples 
count_matrix <- subset(count_matrix, rowSums(count_matrix == 0) == 0) 

# check if count matrix and metadata are consistent in sample order
rownames(samples) <- samples$Run
rownames(count_matrix) <- count_matrix$Geneid
count_matrix <- count_matrix[,2:ncol(count_matrix)] # delete gene id column 
count_matrix <- count_matrix[, rownames(samples)] # reorder count matrix 
all(rownames(samples) == colnames(count_matrix))
  
### separate three treatments and corresponding controls in metadata 
treatment_AhR <- samples[samples$chemical=="3-METHYLCHOLANTHRENE",]
control_AhR <- samples[samples$vehicle=="CMC_.5_%" & samples$mode_of_action=="Control",]
AhR<- rbind(treatment_AhR, control_AhR)
AhR$mode_of_action <- factor(AhR$mode_of_action, levels = c("Control", "AhR")) # convert moa levels to factor 

treatment_CAR <- samples[samples$chemical=="FLUCONAZOLE",]
control_CAR <- samples[samples$vehicle=="CORN_OIL_100_%" & samples$mode_of_action=="Control",]
CAR<- rbind(treatment_CAR, control_CAR)
CAR$mode_of_action <- factor(CAR$mode_of_action, levels = c("Control", "CAR/PXR"))

treatment_PPARA <- samples[samples$chemical=="PIRINIXIC_ACID",]
PPARA <- rbind(treatment_PPARA, control_AhR) # treatment 1 and 3 uses the same vehicle
PPARA$mode_of_action <- factor(PPARA$mode_of_action, levels = c("Control", "PPARA"))

### Run DESeq2 for three chemical treatments 
#-------------------------------AhR-----------------------------------#
dds_AhR <- DESeqDataSetFromMatrix(
  countData = count_matrix[,rownames(AhR)],
  colData = AhR,
  design= ~ mode_of_action
)
dds_AhR <- DESeq(dds_AhR)
res05_AhR <- results(dds_AhR, alpha = 0.05)
summary(res05_AhR)
resultsNames(dds_AhR)
deg_AhR <- subset(res05_AhR_shrink, padj <= 0.05)
dim(deg_AhR)
# Log fold change shrinkage 
# "apeglm" is the adaptive Student's t prior shrinkage estimator from the 'apeglm' package
res05_AhR_shrink <- lfcShrink(dds_AhR, coef = 2, type="apeglm")
summary(res05_AhR_shrink)
# Report the top 10 DE genes sorted by adjusted p-values 
top10_res05_AhR_shrink <- res05_AhR_shrink[order(res05_AhR_shrink$padj, decreasing = F),][1:10,]
#-------------------------------CAR/PXR-----------------------------------#
dds_CAR <- DESeqDataSetFromMatrix(
  countData = count_matrix[,rownames(CAR)],
  colData = CAR,
  design= ~ mode_of_action
)
dds_CAR <- DESeq(dds_CAR)
res05_CAR <- results(dds_CAR, alpha = 0.05)
deg_CAR <- subset(res05_CAR_shrink, padj <= 0.05)
dim(deg_CAR)
res05_CAR_shrink <- lfcShrink(dds_CAR, coef = 2, type="apeglm")
summary(res05_CAR_shrink)
top10_res05_CAR_shrink <- res05_CAR_shrink[order(res05_CAR_shrink$padj, decreasing = F),][1:10,]
#-------------------------------PPARA-----------------------------------#
dds_PPARA <- DESeqDataSetFromMatrix(
  countData = count_matrix[,rownames(PPARA)],
  colData = PPARA,
  design= ~ mode_of_action
)
dds_PPARA <- DESeq(dds_PPARA)
res05_PPARA <- results(dds_PPARA, alpha = 0.05)
deg_PPARA <- subset(res05_PPARA_shrink, padj <= 0.05)
dim(deg_PPARA)
  
res05_PPARA_shrink <- lfcShrink(dds_PPARA, coef = 2, type="apeglm")
summary(res05_PPARA_shrink)
top10_res05_PPARA_shrink <- res05_PPARA_shrink[order(res05_PPARA_shrink$padj, decreasing = F),][1:10,]

# Generate normalized counts matrices
ncm_AhR <- counts(dds_AhR, normalized=T)
ncm_CAR <- counts(dds_CAR, normalized=T)
ncm_PPARA <- counts(dds_PPARA, normalized=T)

# write results to csv
write.csv(deg_AhR,
          file = paste0(outdir, "deg_AhR.csv", sep = ""), row.names = T)
write.csv(top10_res05_AhR_shrink,
          file = paste0(outdir, "top10_res05_AhR_shrink.csv", sep = ""), row.names = T)
write.csv(ncm_AhR,
          file = paste0(outdir, "normalized_count_matrix_AhR.csv", sep = ""), row.names = T)

write.csv(deg_CAR,
          file = paste0(outdir, "deg_CAR.csv", sep = ""), row.names = T)
write.csv(top10_res05_CAR_shrink,
          file = paste0(outdir, "top10_res05_CAR_shrink.csv", sep = ""), row.names = T)
write.csv(ncm_CAR,
          file = paste0(outdir, "normalized_count_matrix_CAR.csv", sep = ""), row.names = T)

write.csv(deg_PPARA,
          file = paste0(outdir, "deg_PPARA.csv", sep = ""), row.names = T)
write.csv(top10_res05_PPARA_shrink,
          file = paste0(outdir, "top10_res05_PPARA_shrink.csv", sep = ""), row.names = T)
write.csv(ncm_PPARA,
          file = paste0(outdir, "normalized_count_matrix_PPARA.csv", sep = ""), row.names = T)
#------------------------------Plotting-----------------------------------#
df_AhR <- as.data.frame(res05_AhR_shrink)
df_CAR <- as.data.frame(res05_CAR_shrink)
df_PPARA <- as.data.frame(res05_PPARA_shrink)
# Create histograms of fold change values from the significant DE genes
par(mfrow=c(1,3))
hist(df_AhR$log2FoldChange[res05_AhR_shrink$padj <= 0.05], 
     xlab = "Log2FC", probability = T, main = "AhR (3ME)", col=NULL, border="blue", cex.lab=1.5, cex.axis=1.5)
hist(df_CAR$log2FoldChange[res05_CAR_shrink$padj <= 0.05], 
     xlab = "Log2FC", probability = T, main = "CAR/PXR (FLU)", col=NULL, border="blue", cex.lab=1.5, cex.axis=1.5)
hist(df_PPARA$log2FoldChange[res05_PPARA_shrink$padj <= 0.05], 
     xlab = "Log2FC", probability = T, main = "PPARA (PIR)", col=NULL, border="blue", cex.lab=1.5, cex.axis=1.5)

# Also create scatter plots of fold change vs nominal p-value
p1 <- EnhancedVolcano(res05_AhR_shrink,
                x = 'log2FoldChange',
                y = 'pvalue',
                lab = rownames(res05_AhR_shrink),
                FCcutoff = log2(1.5),
                pCutoff = 0.05,
                labSize = 0,
                title = "3ME (AhR)",
                titleLabSize = 15,
                subtitle = NULL,
                legendLabSize = 11,
                legendIconSize = 2,
                legendLabels = c('NS', expression(Log[2]~FC), "p-value", 'both'),
                axisLabSize = 14)

p2 <- EnhancedVolcano(res05_CAR_shrink,
                      lab = rownames(res05_AhR_shrink),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      FCcutoff = log2(1.5),
                      pCutoff = 0.05,
                      labSize = 0,
                      title = "FLU (CAR/PXR)",
                      titleLabSize = 15,
                      subtitle = NULL,
                      legendLabSize = 11,
                      legendIconSize = 2,
                      legendLabels = c('NS', expression(Log[2]~FC), "p-value", 'both'),
                      axisLabSize = 14)
p3 <- EnhancedVolcano(res05_PPARA_shrink,
                      lab = rownames(res05_AhR_shrink),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      FCcutoff = log2(1.5),
                      pCutoff = 0.05,
                      labSize = 0,
                      title = "PIR (PPARA)",
                      titleLabSize = 15,
                      subtitle = NULL,
                      legendLabSize = 11,
                      legendIconSize = 2,
                      legendLabels = c('NS', expression(Log[2]~FC), "p-value", 'both'),
                      axisLabSize = 14)
grid.arrange(p1, p2, p3, nrow=1)


                   