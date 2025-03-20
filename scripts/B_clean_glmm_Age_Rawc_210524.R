# Script: Calculate DEGs using glmmSeq
# Date: 20 Dec 2021
#Updated:21 May 2024
##set seed added

# Descriptions:
# This program calculates DEGs using raw counts data
# 
# glmmSeq allows model to include random effects, like patient ID and analyse longitudinal data.
#rm(list=ls())

##### LIBRARIES ####
library(glmmSeq)
library(limma)
library(edgeR)
library(tidyverse)
library(kableExtra)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)


#### Vd2 analysis ####

#### INPUTS ####

#Reading clean expected counts data see script "data_exploration.R" to
#reproduce these files.

counts = read.delim(file = "output_data/counts_Age_Vd2.txt", sep = "\t",
                        header=T, row.names=1)
meta = read.delim(file = "output_data/metadata_Age_Vd2.txt", sep = "\t",
                      header=T, row.names=1)

meta$state_bool <- factor(meta$state_bool, levels = c("unstim","stim"), labels=c("unstim","stim") )
meta$age_bool <- factor(meta$age_bool, levels = c("Child","Adult"), labels=c("Child","Adult") )

meta$state_bool <- relevel(meta$state_bool, ref = "unstim")
meta$age_bool <- relevel(meta$age_bool, ref = "Child")

###MODEL FITTING#####
#Using negative binomial models requires gene dispersion estimates to be made.
#Calculating dispersion with DESeq.
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~state_bool*age_bool)

dds <- DESeq(dds)

dispersions <- setNames(dispersions(dds), rownames(counts))

#DESeq2 dispersion estimates are inversely related to the mean and directly 
#related to variance. Based on this relationship, the dispersion is higher 
#for small mean counts and lower for large mean counts.
res <- results(dds)
res

res <- res[order(res$padj),]
res

pdf(file = "graphs/Vd2 Dispersal estimates_DESeq2.pdf")
plotDispEsts(dds, main="Vd2 Dispersal estimates_DESeq2")
dev.off()

sum(res$padj < 0.05, na.rm=TRUE)

pdf(file = "graphs/Vd2 MA Plot_DESeq2.pdf")
plotMA(res, ylim=c(-2,2), main="Vd2 MA Plot_DESeq2")
dev.off()

#An MA-plot is a scatter plot of log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis)
## The MA-plot shows the log2 fold changes attributable to a given variable over the mean of normalized counts 
## for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. 
## Points which fall out of the window are plotted as open triangles pointing either up or down.


##Estimating size factors 
sizeFactors <- estimateSizeFactorsForMatrix(counts)

###MODEL###
##Running the model state accounting for age
set.seed(1234)
glmm_results <- glmmSeq(~state_bool*age_bool + (1 | samples),
                        countdata = counts,
                        metadata = meta,
                        dispersion = dispersions,
                        progress=TRUE,
                        cores = 8)

#Errors in 5 gene(s): ENSG00000139044, ENSG00000139998, 
#ENSG00000140287, ENSG00000160113, ENSG00000163393
glmm_results <- glmmQvals(glmm_results, cutoff = 0.05, verbose = TRUE)  

# state_bool
# ----------
#   Not Significant     Significant 
# 6378            8404 
# 
# age_bool
# --------
#   Not Significant     Significant 
# 13908             874 
# 
# state_bool:age_bool
# -------------------
#   Not Significant     Significant 
# 11791            2991 

###RESULTS##
##Getting data for plots
stats = data.frame(glmm_results@stats)
predict = data.frame(glmm_results@predict)
stats$gene.id <- rownames(stats)
predict$gene.id <- rownames(predict)

plotdata0 <- left_join(predict, stats, by = "gene.id")
plotdata0$ENSEMBL <- plotdata0$gene.id

##Changing the ensembl ids for gene names
keys(org.Hs.eg.db, keytype="ENSEMBL")[1:10]

anno <- AnnotationDbi::select(org.Hs.eg.db,keys=plotdata0$ENSEMBL,
                              columns=c("SYMBOL","GENENAME"),
                              keytype="ENSEMBL")
# Have a look at the annotation
head(anno)

dup_ids <- anno$ENSEMBL[duplicated(anno$ENSEMBL)]

filter(anno, ENSEMBL %in% dup_ids) %>% 
  arrange(ENSEMBL) %>% head

anno <- AnnotationDbi::select(org.Hs.eg.db,keys=plotdata0$ENSEMBL,
                              columns=c("ENSEMBL","SYMBOL","GENENAME"),
                              keytype="ENSEMBL") %>% filter(!duplicated(ENSEMBL))
  
dim(anno)
#14782     3
##Calculating FC and subsetting for plotting
plotdata <- left_join(plotdata0, anno,by="ENSEMBL")
head(plotdata)

plotdata1 = plotdata %>% mutate(
  #State_AdultFC=y_Adultstim-y_Adult_unstim
  State_AdultFC = log2(plotdata[, 4]+1) - log2(plotdata[, 3]+1),
  #State_ChildFC=y_Childstim-Child_unstim
  State_ChildFC = log2(plotdata[, 2]+1) - log2(plotdata[, 1]+1),
  #Age_stimFC=y_Adultstim - y_Child_stim
  Age_stim = log2(plotdata[, 4]+1) - log2(plotdata[, 2]+1),
  #Age_unstimFC=y_Adultunstim - y_Childunstim
  Age_unstim = log2(plotdata[, 3]+1) - log2(plotdata[, 1]+1))


plotdata2 <- plotdata1[,c(38:44,35:37)]


##Coding for variable SigG
## age = Sig for age but NS for state or state:age
## state = Sig for state but NS for age or state:age
## state:age = Sig for state:age
## ns = NS for age, state or state:age
plotdata2$SigG[plotdata2$qvals.state_bool>= 0.05 & plotdata2$qvals.state_bool.age_bool>=0.05 & plotdata2$qvals.age_bool>=0.05 ] <- "ns"
plotdata2$SigG[plotdata2$qvals.state_bool>= 0.05 & plotdata2$qvals.state_bool.age_bool>=0.05 & plotdata2$qvals.age_bool<=0.05 ] <- "age"
plotdata2$SigG[plotdata2$qvals.state_bool<= 0.05 & plotdata2$qvals.state_bool.age_bool>=0.05 & plotdata2$qvals.age_bool>=0.05 ] <- "state"
plotdata2$SigG[plotdata2$qvals.state_bool<= 0.05 & plotdata2$qvals.state_bool.age_bool>=0.05 & plotdata2$qvals.age_bool<=0.05 ] <- "state&age"
plotdata2$SigG[plotdata2$qvals.state_bool.age_bool<0.05] <- "state:age"
table(plotdata2$SigG, useNA = "ifany")

#age        ns     state state:age state&age 
#57      5959      5704      2991        71


tapply(plotdata2$qvals.state_bool.age_bool, plotdata2$SigG, summary)

# $`state:age`
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000e+00 1.895e-05 1.139e-02 1.543e-02 2.627e-02 4.997e-02 

tapply(plotdata2$qvals.state_bool, plotdata2$SigG, summary)
# $state
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000e+00 1.510e-06 2.989e-04 6.590e-03 7.154e-03 4.996e-02 
tapply(plotdata2$qvals.age_bool, plotdata2$SigG, summary)
# $age
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.007988 0.014822 0.018585 0.029009 0.049102 

##Coding for variable SigG
plotdata2$SigGbyAge[plotdata2$SigG == "state:age" & abs(plotdata2$State_AdultFC) > abs(plotdata2$State_ChildFC)] <- "blue_adult_SG_state:age"
plotdata2$SigGbyAge[plotdata2$SigG == "state:age" & abs(plotdata2$State_AdultFC) < abs(plotdata2$State_ChildFC)] <- "yellow_child_SG_state:age"
table(plotdata2$SigGbyAge)

#blue_adult_SG_state:age yellow_child_SG_state:age 
#294                      2697 

##Getting significant gene lists
plotdata3 <- subset(plotdata2, SigG != "ns")
tapply(plotdata3$qvals.state_bool.age_bool, plotdata3$SigG, summary)
#$`state:age`
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000e+00 1.895e-05 1.139e-02 1.543e-02 2.627e-02 4.997e-02 
##Saving counts and metadata
write.table(plotdata1, file = "output_data/glmmseq_Age_Vd2_output.txt", sep = "\t",
            row.names = TRUE)
write.table(plotdata2, file = "output_data/glmmseq_Age_Vd2_output_summary.txt", sep = "\t",
            row.names = TRUE)
write.table(plotdata3, file = "output_data/glmmseq_Age_Vd2_output_SiGeneList.txt", sep = "\t",
            row.names = TRUE)

##This will clean up the local enviroment to start the monocyte analysis from scratch.
rm(list=ls())
#### Monocytes analysis ####

#### INPUTS ####

#Reading clean expected counts data see script "data_exploration.R" to
#reproduce these files.
counts= read.delim(file = "output_data/counts_Age_monos.txt", sep = "\t",
                          header=T, row.names=1)
meta= read.delim(file = "output_data/metadata_Age_monos.txt", sep = "\t",
                        header=T, row.names=1)

###MODEL FITTING#####
#Using negative binomial models requires gene dispersion estimates to be made.
#Calculating dispersion with DESeq.

meta$state_bool <- factor(meta$state_bool, levels = c("unstim","stim"), labels=c("unstim","stim") )
meta$age_bool <- factor(meta$age_bool, levels = c("Child","Adult"), labels=c("Child","Adult") )

meta$state_bool <- relevel(meta$state_bool, ref = "unstim")
meta$age_bool <- relevel(meta$age_bool, ref = "Child")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~state_bool*age_bool)

dds <- DESeq(dds)
dispersions <- setNames(dispersions(dds), rownames(counts))

#DESeq2 dispersion estimates are inversely related to the mean and directly 
#related to variance. Based on this relationship, the dispersion is higher 
#for small mean counts and lower for large mean counts.
res <- results(dds)
res

res <- res[order(res$padj),]
res

pdf(file="graphs/Monocytes Dispersal estimates_DESeq2.pdf")
plotDispEsts(dds, main="Monocytes Dispersal estimates_DESeq2")
dev.off()

sum(res$padj < 0.05, na.rm=TRUE)

pdf(file="graphs/Monocytes MA Plot_DESeq2.pdf")
plotMA(res, ylim=c(-2,2), main="Monocytes MA Plot_DESeq2")
dev.off()

#An MA-plot is a scatter plot of log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis)
## The MA-plot shows the log2 fold changes attributable to a given variable over the mean of normalized counts 
## for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. 
## Points which fall out of the window are plotted as open triangles pointing either up or down.

##Estimating size factors 
sizeFactors <- estimateSizeFactorsForMatrix(counts)

###MODEL###
set.seed(1234)
##Running the model state accounting for age
glmm_results <- glmmSeq(~ state_bool * age_bool + (1 | samples),
                        countdata = counts,
                        metadata = meta,
                        dispersion = dispersions,
                        progress=TRUE,
                        cores = 8)
#Errors in 1 gene(s): ENSG00000146453
glmm_results <- glmmQvals(glmm_results, cutoff = 0.05, verbose = TRUE)  


# state_bool
# ----------
#   Not Significant     Significant 
# 5472            9801 
# 
# age_bool
# --------
#   Not Significant     Significant 
# 14469             804 
# 
# state_bool:age_bool
# -------------------
#   Not Significant     Significant 
# 13002            2271 


###RESULTS##

###RESULTS##
##Getting data for plots
stats = data.frame(glmm_results@stats)
predict = data.frame(glmm_results@predict)
stats$gene.id <- rownames(stats)
predict$gene.id <- rownames(predict)

##Changing the ensembl ids for gene names
plotdata0 <- left_join(predict, stats, by = "gene.id")
dim(plotdata0)
#15273    38
plotdata0$ENSEMBL <- plotdata0$gene.id

keys(org.Hs.eg.db, keytype="ENSEMBL")[1:10]

anno <- AnnotationDbi::select(org.Hs.eg.db,keys=plotdata0$ENSEMBL,
                              columns=c("SYMBOL","GENENAME"),
                              keytype="ENSEMBL")
# Have a look at the annotation
head(anno)

dup_ids <- anno$ENSEMBL[duplicated(anno$ENSEMBL)]

filter(anno, ENSEMBL %in% dup_ids) %>% 
  arrange(ENSEMBL) %>% head

anno <- AnnotationDbi::select(org.Hs.eg.db,keys=plotdata0$ENSEMBL,
                              columns=c("ENSEMBL","SYMBOL","GENENAME"),
                              keytype="ENSEMBL") %>% filter(!duplicated(ENSEMBL))
  
dim(anno)
#15273     3
##Calculating FC and subsetting data for plots
plotdata <- left_join(plotdata0, anno,by="ENSEMBL")
head(plotdata)

plotdata1 = plotdata %>% mutate(
  #State_AdultFC=y_Adultstim-y_Adult_unstim
  State_AdultFC = log2(plotdata[, 4]+1) - log2(plotdata[, 3]+1),
  #State_ChildFC=y_Childstim-Child_unstim
  State_ChildFC = log2(plotdata[, 2]+1) - log2(plotdata[, 1]+1),
  #Age_stimFC=y_Adultstim - y_Child_stim
  Age_stim = log2(plotdata[, 4]+1) - log2(plotdata[, 2]+1),
  #Age_unstimFC=y_Adultunstim - y_Childunstim
  Age_unstim = log2(plotdata[, 3]+1) - log2(plotdata[, 1]+1))


plotdata2 <- plotdata1[,c(38:44,35:37)]

##Coding for variable SigG
## age = Sig for age but NS for state or state:age
## state = Sig for state but NS for age or state:age
## state:age = Sig for state:age
## ns = NS for age, state or state:age
# state&age = Sig for state and age but not for interaction
plotdata2$SigG[plotdata2$qvals.state_bool>= 0.05 & plotdata2$qvals.state_bool.age_bool>=0.05 & plotdata2$qvals.age_bool>=0.05 ] <- "ns"
plotdata2$SigG[plotdata2$qvals.state_bool>= 0.05 & plotdata2$qvals.state_bool.age_bool>=0.05 & plotdata2$qvals.age_bool<=0.05 ] <- "age"
plotdata2$SigG[plotdata2$qvals.state_bool<= 0.05 & plotdata2$qvals.state_bool.age_bool>=0.05 & plotdata2$qvals.age_bool>=0.05 ] <- "state"
plotdata2$SigG[plotdata2$qvals.state_bool<= 0.05 & plotdata2$qvals.state_bool.age_bool>=0.05 & plotdata2$qvals.age_bool<=0.05 ] <- "state&age"
plotdata2$SigG[plotdata2$qvals.state_bool.age_bool<0.05] <- "state:age"
table(plotdata2$SigG, useNA = "ifany")
#age        ns     state state:age state&age 
#167      4711      7660      2271       464 

tapply(plotdata2$qvals.state_bool.age_bool, plotdata2$SigG, summary)
#$`state:age`
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.009385 0.020698 0.022144 0.034025 0.049928 

tapply(plotdata2$qvals.state_bool, plotdata2$SigG, summary)
#$state
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000e+00 0.000e+00 6.591e-05 6.223e-03 5.420e-03 4.998e-02

tapply(plotdata2$qvals.age_bool, plotdata2$SigG, summary)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#4.600e-07 2.019e-02 3.053e-02 2.892e-02 4.112e-02 4.978e-02

##Coding for variable SigG
plotdata2$SigGbyAge[plotdata2$SigG == "state:age" & abs(plotdata2$State_AdultFC) > abs(plotdata2$State_ChildFC)] <- "blue_adult_SG_state:age"
plotdata2$SigGbyAge[plotdata2$SigG == "state:age" & abs(plotdata2$State_AdultFC) < abs(plotdata2$State_ChildFC)] <- "yellow_child_SG_state:age"
table(plotdata2$SigGbyAge)
#blue_adult_SG_state:age yellow_child_SG_state:age 
#1262                      1009

##Getting significant gene lists
plotdata3 <- subset(plotdata2, SigG != "ns")
tapply(plotdata3$qvals.state_bool.age_bool, plotdata3$SigG, summary)
##Saving counts and metadata
write.table(plotdata1, file = "output_data/glmmseq_Age_Monocytes_output.txt", sep = "\t",
            row.names = TRUE)
write.table(plotdata2, file = "output_data/glmmseq_Age_Monocytes_output_summary.txt", sep = "\t",
            row.names = TRUE)
write.table(plotdata3, file = "output_data/glmmseq_Age_Monocytes_output_SiGeneList.txt", sep = "\t",
            row.names = TRUE)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

