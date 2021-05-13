mydir = '/projectnb/bf528/users/dreadlocks/project_3/analyst/'
outdir = '/projectnb/bf528/users/dreadlocks/project_3/analyst/output/'
setwd(mydir)
library(limma)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
#----------------------Load Files---------------------------------#
# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_6_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)

#------------------------Subset------------------------------------------#
# sample subsets based on treatment
# find appropiate controls for each treatment 

treatment1 <- samples[samples$chemical=="3-METHYLCHOLANTHRENE",]
control1 <- samples[samples$vehicle=="CMC_.5_%" & samples$chemical=="Control",]
group1<- rbind(treatment1, control1)

treatment2 <- samples[samples$chemical=="FLUCONAZOLE",]
control2 <- samples[samples$vehicle=="CORN_OIL_100_%" & samples$chemical=="Control",]
group2<- rbind(treatment2, control2)

treatment3 <- samples[samples$chemical=="PIRINIXIC_ACID",]
group3 <- rbind(treatment3, control1) # treatment 1 and 3 uses the same vehicle

# rma expression matrix subset
rmasub1 <- rma[paste0('X',group1$array_id)]
rmasub2 <- rma[paste0('X',group2$array_id)]
rmasub3 <- rma[paste0('X',group3$array_id)]

#-----------------------------Design Matrix------------------------------#
# construct a design matrix modeling treatment vs control for use by limma
design1 <- model.matrix(~factor(group1$chemical, levels = c("Control", "3-METHYLCHOLANTHRENE")))
design2 <- model.matrix(~factor(group2$chemical, levels = c("Control", "FLUCONAZOLE")))
design3 <- model.matrix(~factor(group3$chemical, levels = c("Control", "PIRINIXIC_ACID")))
colnames(design1) <- c('Intercept','3-METHYLCHOLANTHRENE')
colnames(design2) <- c('Intercept','FLUCONAZOLE')
colnames(design3) <- c('Intercept', "PIRINIXIC_ACID")

#---------------------------# run limma---------------------------------#
findDE <- function(x,y, filename){
  fit <- eBayes(lmFit(x,y))
  deres <- topTable(fit, coef=2, n=nrow(x), adjust='BH', sort.by = "P", p.value = 0.05)
 
  write.csv(deres, paste(outdir, filename, sep = ''))
  write.csv(deres[1:10,], paste(outdir, 'top10', filename, sep = ''))
  return(deres)
}

deres_METHYLCHOLANTHRENE <- findDE(rmasub1, design1, "METHYLCHOLANTHRENE.csv")
deres_FLUCONAZOLE<- findDE(rmasub2, design2, "FLUCONAZOLE.csv")
deres_PIRINIXIC_ACID <- findDE(rmasub3, design3, "PIRINIXIC_ACID.csv")


##--------------------------plotting----------------------------##
deres_METHYLCHOLANTHRENE$pass_filter <- "NO"
deres_METHYLCHOLANTHRENE$pass_filter[deres_METHYLCHOLANTHRENE$logFC > 1.5 & deres_METHYLCHOLANTHRENE$P.Val < 0.05] <- "UP"
deres_METHYLCHOLANTHRENE$pass_filter[deres_METHYLCHOLANTHRENE$logFC <= -1.5 & deres_METHYLCHOLANTHRENE$P.Val < 0.05] <- "DOWN"
deres_METHYLCHOLANTHRENE$pass_filter <- factor(deres_METHYLCHOLANTHRENE$pass_filter)

deres_FLUCONAZOLE$pass_filter <- "NO"
deres_FLUCONAZOLE$pass_filter[deres_FLUCONAZOLE$logFC > 1.5 & deres_FLUCONAZOLE$P.Val < 0.05] <- "UP"
deres_FLUCONAZOLE$pass_filter[deres_FLUCONAZOLE$logFC <= 1.5 & deres_FLUCONAZOLE$P.Val < 0.05] <- "DOWN"
deres_FLUCONAZOLE$pass_filter <- factor(deres_FLUCONAZOLE$pass_filter)

deres_PIRINIXIC_ACID$pass_filter <- "NO"
deres_PIRINIXIC_ACID$pass_filter[deres_PIRINIXIC_ACID$logFC > 1.5 & deres_PIRINIXIC_ACID$P.Val < 0.05] <- "UP"
deres_PIRINIXIC_ACID$pass_filter[deres_PIRINIXIC_ACID$logFC <= 1.5 & deres_PIRINIXIC_ACID$P.Val < 0.05] <- "DOWN"
deres_PIRINIXIC_ACID$pass_filter <- factor(deres_PIRINIXIC_ACID$pass_filter)

# Histograms of fold change values
hist1 <- ggplot(data = deres_METHYLCHOLANTHRENE, aes(x=logFC)) + geom_histogram(color="blue", fill=NA, bins = 15) + ggtitle("3-Methylcholanthrene") + theme_minimal() + theme(plot.title = element_text(size=10))
hist2 <- ggplot(data = deres_FLUCONAZOLE, aes(x=logFC)) + geom_histogram(color="blue", fill=NA, bins = 15) + ggtitle("fluconazole") + theme_minimal() + theme(plot.title = element_text(size=10))
hist3 <- ggplot(data = deres_PIRINIXIC_ACID, aes(x=logFC)) + geom_histogram(color="blue", fill=NA, bins = 15) + ggtitle("pirinixic acid") + theme_minimal() + theme(plot.title = element_text(size=10))
grid.arrange(hist1, hist2, hist3, nrow=1, top = "Histograms of Log2 Fold Change") 

par(mfrow=c(1,3))
hist(deres_METHYLCHOLANTHRENE$logFC, xlab = "Log2FC", probability = T, main = "3-METHYLCHOLANTHRENE")
hist(deres_FLUCONAZOLE$logFC, xlab = "Log2FC", probability = T, main = "FLUCONAZOLE")
hist(deres_PIRINIXIC_ACID$logFC, xlab = "Log2FC", probability = T, main = "PIRINIXIC_ACID")


# scatter plots of fold change vs nominal p-value
sc1 <- ggplot(data = deres_METHYLCHOLANTHRENE, aes(x=logFC, y=-log10(P.Value))) + geom_point(shape=1) + ggtitle("3-Methylcholanthrene") + theme_minimal() + theme(plot.title = element_text(face = "bold", size=10)) 
sc2 <- ggplot(data =deres_FLUCONAZOLE, aes(x=logFC, y=-log10(P.Value))) + geom_point(shape=1) + ggtitle("fluconazole") + theme_minimal() + theme(plot.title = element_text(face = "bold", size=10))
sc3 <- ggplot(data =deres_PIRINIXIC_ACID, aes(x=logFC, y=-log10(P.Value))) + geom_point(shape=1) + ggtitle("pirinixic acid") + theme_minimal() + theme(plot.title = element_text(face = "bold", size=10))
grid.arrange(sc1, sc2, sc3, nrow=1)


