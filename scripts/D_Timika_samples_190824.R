#Author: Z Pava

library('edgeR')
library('BiocParallel')
library('tidyverse')
library('statmod')
library(here)
##### LIBRARIES glmseq####
library(glmmSeq)
library(limma)
library(edgeR)
library(tidyverse)
library(kableExtra)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
#####PCA ###
library(tidyverse)
library(here)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(edgeR)
library(DESeq2)  # For normalization
library(factoextra) #for PCA visualization

############INPUT###########

##Getting the HTSeq raw data counts
counts0 <- read_tsv(here("data/GSE154317_raw_counts_GRCh38.p13_NCBI.tsv"))
dim(counts0)
head(counts0)
counts0 <- as.data.frame(counts0)

##changing the gsm number for samples id obtained from ncbi
samples_id <- read.csv(here("data/samples_id_timika_ncbi.csv"))
head(samples_id)
rownames(counts0) <- counts0$GeneID

# Rename columns using the key
counts <- counts0 %>%
  dplyr::select(-GeneID)%>%
  rename_with(~ samples_id$sample[match(., samples_id$ncbi)], everything())
head(counts)

##Getting gene names for gene ids
counts$ENTREZID <- rownames(counts)

#### Changing ENSEMBL names to gene IDS ####

# Sample data
plotdata0 <- counts
dim(plotdata0)
# 39376    47

# Annotate ENSEMBL IDs with gene symbols and gene names
anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              keys = plotdata0$ENTREZID,
                              columns = c("SYMBOL", "GENENAME"),
                              keytype = "ENTREZID")

# Have a look at the annotation
head(anno)

# Handle duplicates
# Identify duplicate ENSEMBL IDs
dup_ids <- anno$ENTREZID[duplicated(anno$SYMBOL)]

# Filter out duplicated ENSEMBL IDs
anno <- anno %>% filter(!is.na(SYMBOL))

# Merge plotdata0 with the annotation data using  IDs
plotdata <- left_join(plotdata0, anno, by = "ENTREZID")

# Average counts for duplicate gene symbols
plotdata_avg <- plotdata %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), ~ round(mean(.x, na.rm = TRUE)))) %>%
  ungroup() %>% filter(!is.na(SYMBOL)) %>%
  as.data.frame(row.names = TRUE)

dim(plotdata_avg)
# 37765    47

##Adding rownames
rownames(plotdata_avg) <- plotdata_avg$SYMBOL

counts_sF <- plotdata_avg %>%
  dplyr::select(-SYMBOL)%>%
  as.data.frame()

head(counts_sF)
dim(counts_sF)
#37765    46

##selecting only endemic samples
counts_F <- counts_sF %>%
  dplyr::select(-contains("CHMI"))
dim(counts_F)
#37765    36
##creating metadata
samplesname <- as.character(colnames(counts_F))
meta1<- as.data.frame(counts_F[1:36, c(1:length(counts_F))], row.names = samplesname)
meta1$iud <- rownames(meta1)
meta <- as.data.frame(meta1[1:36, 37], row.names = samplesname)
colnames(meta)[1] = "uid"
head(meta)

#Extracting meta data from uid
meta %>% extract(col = uid, into = "samples", regex = "(Adult\\d+|Child\\d+)", remove = FALSE)-> meta
meta %>% extract(col = uid, into = "age", regex = "(Adult|Child)", remove = FALSE) -> meta
meta %>% extract(col = uid, into = "state", regex = "(day0|day28)", remove = FALSE) -> meta

##give proper names
meta$state <- gsub("day0", "stim", meta$state)
meta$state <- gsub("day28", "unstim", meta$state)
##creating groups
meta$group <- paste0(meta$age, "_", meta$state)

##Creating booleans

meta = meta %>% mutate(
  state_bool = case_when(state == "unstim" ~ "0",
                         state == "stim" ~ "1"),
  age_bool = case_when(age == "Child" ~ "0",
                       age == "Adult" ~ "1"),
  group_bool = case_when(group == "Child_stim" ~ "3",
                         group == "Adult_stim" ~ "1",
                         group == "Child_unstim" ~ "2",
                         group == "Adult_unstim" ~ "0")
)

meta$state_bool <- factor(meta$state_bool, levels = c(0,1), labels=c("unstim","stim") )
meta$age_bool <- factor(meta$age_bool, levels = c(0,1), labels=c("Child","Adult") )


##### Doing QC ###

##Create DGEList object
d0 <- DGEList(counts_F)

##Calculate normalization factors:
#calcNormFactors doesn’t normalize the data, it just calculates 
#normalization factors for use downstream.
d0 <- calcNormFactors(d0)
d0

####Filtering low-expressed genes####
#Removing genes that are lowly expressed
#Here we perform a minimal pre-filtering to keep only rows that have 
#at least 10 reads total. Note that more strict filtering to increase 
#power is automatically applied via independent filtering on the mean of 
#normalized counts within the results function

keep <- rowSums(d0[[1]]) >= 10
d <- d0[keep,]
dim(d)
#23316    36

#A mean-difference plot (MD-plot) is a plot of log-intensity ratios (differences) 
#versus log-intensity averages (means). Useful to detect outliers.

for (i in 1:36) {
  pdf(file= paste("graphs/MD_TimikaMonos_",i,".pdf", sep = ""))
  plotMD(cpm(d[[1]]), column=i, xlab = "Average log-expression",
         ylab = "Expression log-ratio (this sample vs others)") 
  abline(h=0,col="red",lty=2,lwd=2)
  dev.off()
}

group <- factor(meta$group) 
group

##Multidimentional scaling (MDS) plot
#Visualizes the differences between the expression profiles 
#of different samples in two dimensions

pdf(file="graphs/MDS_TimikaMonos_pergroup.pdf")
points <- c(0:35) 
colors <- rep(c("blue","darkgreen","red", "purple"),2)
plotMDS(cpm(d[[1]]), col=colors[group], pch=points[group], xlab = "leading FC dim1", ylab = "leading FC dim2") 
legend("topleft",legend=levels(group),pch=points,col=colors,ncol=2)
dev.off()

##Saving counts and metadata
d[[1]] -> counts_timika


write.table(counts_timika, file = "output_data/counts_Age_Timika_monos.txt", sep = "\t",
            row.names = TRUE)
write.table(meta, file = "output_data/meta_Age_Timika_monos.txt", sep = "\t",
            row.names = TRUE)

#### GLMMSEQ ####

#Reading clean expected counts data see script "data_exploration.R" to
#reproduce these files.

counts = read.delim(file = here("output_data/counts_Age_Timika_monos.txt"), sep = "\t", header=T, row.names=1)
meta = read.delim(file = here("output_data/meta_Age_Timika_monos.txt"), sep = "\t", header=T, row.names=1)

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

pdf(file = "graphs/Timika_monocytes_Dispersal estimates_DESeq2.pdf")
plotDispEsts(dds, main="Timika_monocytes Dispersal estimates_DESeq2")
dev.off()

sum(res$padj < 0.05, na.rm=TRUE)

pdf(file = "graphs/Timika_monocytes Plot_DESeq2.pdf")
plotMA(res, ylim=c(-7,7), main="Timika_monocytes MA Plot_DESeq2")
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
                        cores = 32)

glmm_results <- glmmQvals(glmm_results, cutoff = 0.05, verbose = TRUE)  

# Errors in 9 gene(s): FABP4, GSTM5, IGHD4-4, IGKV1-13, IGLV3-1, KIF1A, MIR650, RNF17, SMIM43> 
#   > glmm_results <- glmmQvals(glmm_results, cutoff = 0.05, verbose = TRUE)  
# 
# state_bool
# ----------
#   Not Significant     Significant 
# 22802             505 
# 
# age_bool
# --------
#   Not Significant     Significant 
# 22943             364 
# 
# state_bool:age_bool
# -------------------
#   Not Significant     Significant 
# 23177             130 
###RESULTS##
##Getting data for plots
stats = data.frame(glmm_results@stats)
predict = data.frame(glmm_results@predict)
stats$gene.id <- rownames(stats)
predict$gene.id <- rownames(predict)

plotdata0 <- left_join(predict, stats, by = "gene.id")
head(plotdata0)


plotdata1 = plotdata0 %>% mutate(
  #State_AdultFC=y_Adultstim-y_Adult_unstim
  State_AdultFC = log2(plotdata0[, 4]+1) - log2(plotdata0[, 3]+1),
  #State_ChildFC=y_Childstim-Child_unstim
  State_ChildFC = log2(plotdata0[, 2]+1) - log2(plotdata0[, 1]+1),
  #Age_stimFC=y_Adultstim - y_Child_stim
  Age_stim = log2(plotdata0[, 4]+1) - log2(plotdata0[, 2]+1),
  #Age_unstimFC=y_Adultunstim - y_Childunstim
  Age_unstim = log2(plotdata0[, 3]+1) - log2(plotdata0[, 1]+1))

colnames(plotdata1)

plotdata2 <- plotdata1[,c(13,38:41,35:37)]
head(plotdata2)

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
#244     22521       379       130        33


tapply(plotdata2$qvals.state_bool.age_bool, plotdata2$SigG, summary)

#$`state:age`
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000000 0.0000000 0.0000000 0.0046379 0.0001136 0.0498185 

tapply(plotdata2$qvals.state_bool, plotdata2$SigG, summary)
#$state
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.005415 0.020782 0.020585 0.035638 0.049974
tapply(plotdata2$qvals.age_bool, plotdata2$SigG, summary)
#$age
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.007337 0.021599 0.023035 0.040052 0.049761  

##Coding for variable SigG
plotdata2$SigGbyAge[plotdata2$SigG == "state:age" & abs(plotdata2$State_AdultFC) > abs(plotdata2$State_ChildFC)] <- "blue_adult_SG_state:age"
plotdata2$SigGbyAge[plotdata2$SigG == "state:age" & abs(plotdata2$State_AdultFC) < abs(plotdata2$State_ChildFC)] <- "yellow_child_SG_state:age"
table(plotdata2$SigGbyAge)

#blue_adult_SG_state:age yellow_child_SG_state:age 
#49                        81

##Getting significant gene lists
plotdata3 <- subset(plotdata2, SigG != "ns")
tapply(plotdata3$qvals.state_bool.age_bool, plotdata3$SigG, summary)
#$`state:age`
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000000 0.0000000 0.0000000 0.0046379 0.0001136 0.0498185 
##Saving counts and metadata
write.table(plotdata1, file = here("output_data/glmmseq_Timika_monos_output.txt"), sep = "\t",
            row.names = TRUE)
write.table(plotdata2, file = here("output_data/glmmseq_Timika_monos_output_summary.txt"), sep = "\t",
            row.names = TRUE)
write.table(plotdata3, file = here("output_data/glmmseq_Timika_monos_output_SiGeneList.txt"), sep = "\t",
            row.names = TRUE)

####Doing PCAs ####
##selecting genes age and state:age interaction
monos_de_genes <- read.table(here("output_data/glmmseq_Timika_monos_output_SiGeneList.txt"))
monos_de_age <- monos_de_genes%>%
  filter(SigG!= "state")
dim(monos_de_age)
#407 10
keep_monos <- monos_de_age$gene.id
# Selecting de for each celltype
counts_pca <- counts_F[keep_monos, , drop=FALSE]
dim(counts_pca)
#407  36
####Doing a PCA?
group <- factor(meta$group) 
group
#create deseq object
dds <- DESeqDataSetFromMatrix(countData = counts_pca,
                              colData = meta,
                              design = ~state_bool*age_bool)
# Normalize the data using DESeq2's variance stabilizing transformation (VST) or regularized log transformation (rlog)
vsd <- vst(dds, blind = TRUE)  # Variance stabilizing transformation

# Extract normalized expression data
normalized_counts <- assay(vsd)

# Perform PCA using prcomp
pca <- prcomp(t(normalized_counts), scale. = TRUE)  # transpose to have samples as rows

#### Visualize the PCA results####
# Create a data frame for ggplot2
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Sample = colnames(counts_F), Group = meta$group, Age = meta$age)

# Plot
# Ensure the Age column is a factor with levels "Adults" and "Child"
pca_data$Age <- factor(pca_data$Age, levels = c("Adult", "Child"))
#col by sample
colt <- c("Adult1.day28" = "lightblue","Adult1.day0" = "darkblue", "Child1.day28" = "#FF9999", "Child1.day0" = "darkred",
          "Adult2.day28" = "lightblue","Adult2.day0" = "darkblue", "Child2.day28" = "#FF9999", "Child2.day0" = "darkred",
          "Adult3.day28" = "lightblue","Adult3.day0" = "darkblue", "Child3.day28" = "#FF9999", "Child3.day0" = "darkred",
          "Adult4.day28" = "lightblue","Adult4.day0" = "darkblue", "Child4.day28" = "#FF9999", "Child4.day0" = "darkred",
          "Adult5.day28" = "lightblue","Adult5.day0" = "darkblue", "Child5.day28" = "#FF9999", "Child5.day0" = "darkred",
          "Adult6.day28" = "lightblue","Adult6.day0" = "darkblue", "Child6.day28" = "#FF9999", "Child6.day0" = "darkred",
          "Adult7.day28" = "lightblue","Adult7.day0" = "darkblue", "Child7.day28" = "#FF9999", "Child7.day0" = "darkred",
          "Adult8.day28" = "lightblue","Adult8.day0" = "darkblue", "Child8.day28" = "#FF9999", "Child8.day0" = "darkred",
          "Adult10.day28" = "lightblue","Adult10.day0" = "darkblue")
# Create the PCA plot colored by Age
dds <- DESeqDataSetFromMatrix(countData = counts_pca,
                              colData = DataFrame(condition = factor(meta$group)),
                              ~ condition)
# Normalize the data using DESeq2's variance stabilizing transformation (VST) or regularized log transformation (rlog)

##too few genes
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)  # Variance stabilizing transformation

# Extract normalized expression data
normalized_counts <- assay(vsd)

# Perform PCA using prcomp
pca <- prcomp(t(normalized_counts), scale. = TRUE)  # transpose to have samples as rows

#### Visualize the PCA results####
# Create a data frame for ggplot2
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Sample = colnames(counts_pca), Group = meta$group, Age = meta$age)
# Plot
# Ensure the Age column is a factor with levels "Adults" and "Child"
pca$Age <- factor(pca$Age, levels = c("Adult", "Child"))
pca$Sample <- factor(pca$Sample)

colz <- c("Adult_unstim" = "lightblue","Adult_stim" = "darkblue", "Child_unstim" = "#FF9999", "Child_stim" = "darkred")
#
pca_deage <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(summary(pca)$importance[2,1] * 100, 2), "% variance")) +
  ylab(paste0("PC2: ", round(summary(pca)$importance[2,2] * 100, 2), "% variance")) +
  scale_color_manual(values = colz) +
  theme_bw()

# Plot the PCA
pca_deage
# Save it
ggsave(pca_deage, filename = here("graphs/PCA_deg_Timika_Age_groups_Monos.pdf"),  width = 7, height = 5, dpi = 600,scale = 1)

#extracting information
#Here we’ll show how to calculate the PCA results for variables: coordinates, cos2 and contributions:
#var.coord = loadings * the component standard deviations
#var.cos2 = var.coord^2
#var.contrib. The contribution of a variable to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component)

# Helper function 
#::::::::::::::::::::::::::::::::::::::::
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- pca$rotation
sdev <- pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
head(var.coord[, 1:4])

# Compute Cos2
#::::::::::::::::::::::::::::::::::::::::
var.cos2 <- var.coord^2
head(var.cos2[, 1:4])

# Compute contributions
#::::::::::::::::::::::::::::::::::::::::
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[, 1:4])

####Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.####
barp_var <- fviz_eig(pca)
ggsave(barp_var, filename = here("graphs/barplot_expvardim_Timikamonos.190824.pdf"),  width = 7, height = 5, dpi = 600,scale = 1)

####Graph biologically interesting genes contributors of variance.####
int_genes <- c("S100A16", "CABP4", "PIGQ", "ELAPOR1", "CXCL11", "APOA2",
               "ARG1", "BCL2L1","FCGR2B", "HLA-DRA", "PIK3CD", "IL27RA" )
####Graph biologically interesting genes contributors of variance.####
int_genes_db <- var.cos2[int_genes,]

pca_var <- fviz_pca_var(pca,
                        col.var = "contrib", # Color by contributions to the PC
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE,     # Avoid text overlapping
                        select.var = list(names=int_genes)
)
pca_var
pca_var2 <- fviz_pca_var(pca, col.var = "contrib", # Color by contributions to the PC
                         axes = c(1,2),
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE,     # Avoid text overlapping
                         select.var = list(contrib =20)
)
pca_var2
##cos2 -1 how well a variable is represented in pca. closer to 1 the better
pca_var3 <- fviz_pca_var(pca, col.var = "cos2", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE,     # Avoid text overlapping
                         select.var = list(cos2 = 20)
)
pca_var3
ggsave(pca_var, filename = here("graphs/PCA_selectvar_TimikaMonos.190824.pdf"),  width = 7, height = 5, dpi = 600,scale = 1)
ggsave(pca_var2, filename = here("graphs/PCA_varcontrib_TimikaMonos.190824.pdf"),  width = 7, height = 5, dpi = 600,scale = 1)
ggsave(pca_var3, filename = here("graphs/PCA_varcos2_TimikaMonos.190824.pdf"),  width = 7, height = 5, dpi = 600,scale = 1)


