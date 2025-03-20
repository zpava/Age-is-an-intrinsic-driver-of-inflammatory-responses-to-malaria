#Author: Z Pava

library('edgeR')
library('BiocParallel')
library('tidyverse')
library('statmod')

############INPUT###########

##Getting the HTSeq raw data counts
counts0 <- read.csv("data/Boyle_-_RNASeq_NextSeq550_HTSeq.genes.RawCounts.csv", header=T, row.names = 1)

samplesname <- as.character(colnames(counts0[, c(8:47)]))
meta1<- as.data.frame(counts0[1:40, c(8:length(counts0))], row.names = samplesname)
meta1$iud <- rownames(meta1)
meta <- as.data.frame(meta1[1:40, 41], row.names = samplesname)
colnames(meta)[1] = "uid"

#Extracting meta data from uid
meta %>% extract(col = uid, into = "samples", regex = "(RDH\\d+)", remove = FALSE) -> meta
meta %>% extract(col = uid, into = "age", regex = "(Adult|Child\\w)", remove = FALSE) -> meta
meta %>% extract(col = uid, into = "celltype", regex = "(Vd2|Class\\w)", remove = FALSE) -> meta
meta %>% extract(col = uid, into = "state", regex = "(Exviv|4hpRB\\w)", remove = FALSE) -> meta

meta$state <- gsub("Exviv", "unstim", meta$state)
meta$state <- gsub("4hpRBC", "stim", meta$state)
meta$celltype <- gsub("Classm", "Monocytes", meta$celltype)
meta$age <- gsub("Child4", "Child", meta$age)
meta$age <- gsub("ChildE", "Child", meta$age)
meta$group <- factor(paste0(meta$age, sep = ".", meta$state))
meta = meta %>% mutate(
  state_bool = case_when(state == "unstim" ~ "0",
                         state == "stim" ~ "1"),
  age_bool = case_when(age == "Child" ~ "0",
                       age == "Adult" ~ "1"),
  age_state_bool = case_when(age == "Child" & state == "stim" ~ "3",
                             age == "Adult" & state == "stim" ~ "1",
                             age == "Child" & state == "unstim" ~ "2",
                             age == "Adult" & state == "unstim" ~ "0")
)

meta$state_bool <- factor(meta$state_bool, levels = c(0,1), labels=c("unstim","stim") )
meta$age_bool <- factor(meta$age_bool, levels = c(0,1), labels=c("Child","Adult") )

#meta = arrange(meta, age, state)

#Subsetting per celltype.

meta_Vd2 = meta %>% filter(grepl(pattern = "Vd2", x = celltype))  
meta_monos = meta %>% filter(grepl(pattern = "Monocytes", x = celltype))

##counts
counts0$ENSEMBL <- rownames(counts0)
counts1 <- counts0 %>% filter(grepl(pattern = "protein_coding", x = Gene.Biotype))
##Removing duplicates for now. 
#counts1 <- counts1[!(duplicated(counts1$Associated.Gene.Name)), ]
dim(counts1)
#20356 48
counts1 <- counts1[!(duplicated(counts1$ENSEMBL)), ]

##change rownames
#rownames(counts1) <- counts1$Associated.Gene.Name

##Subsetting per cell type and sorting. rownames(metadata) must be in the 
#same order as colnames(counts)

counts_Vd2 <- counts1[, grep("Vd2$", colnames(counts1))]

counts_monos <- counts1[, grep("monos$", colnames(counts1))]

######Vd2 data######
# 
# ###Quality check
counts <- counts_Vd2
meta <- meta_Vd2

##Create DGEList object
d0 <- DGEList(counts)

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


#A mean-difference plot (MD-plot) is a plot of log-intensity ratios (differences) 
#versus log-intensity averages (means). Useful to detect outliers.

for (i in 1:20) {
  pdf(file= paste("graphs/MD_Vd2_",i,".pdf", sep = ""))
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

pdf(file="graphs/MDS_Vds_pergroup.pdf")
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19) 
colors <- rep(c("blue","darkgreen","red", "purple"),2)
plotMDS(cpm(d[[1]]), col=colors[group], pch=points[group], xlab = "leading FC dim1", ylab = "leading FC dim2") 
legend("center",legend=levels(group),pch=points,col=colors,ncol=2)
dev.off()

##Plot with samples ID
sam <- as.data.frame(colnames(counts_Vd2))
sam %>% extract(col = `colnames(counts_Vd2)`, into = "samples", regex = "(RDH\\d+)", remove = FALSE) -> sam
sam %>% extract(col = `colnames(counts_Vd2)`, into = "state", regex = "(Exviv|4hpRB\\w)", remove = FALSE) -> sam
sam$state <- gsub("Exviv", "unstim", sam$state)
sam$state <- gsub("4hpRBC", "stim", sam$state)
sam$plotid = paste0(sam$samples, "_", sam$state)
head(sam)

plot_counts <- counts
colnames(plot_counts) <- sam$plotid
sams <- factor(sam$plotid)

pdf(file="graphs/MDS_Vds_persample.pdf")
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19) 
colors <- rep(c("blue","darkgreen","red","purple","orange","brown","salmon","black","turquoise","darkgrey"),2)
plotMDS(cpm(d[[1]]), col=colors[sams], pch=points[sams], xlab = "leading FC dim1", ylab = "leading FC dim2") 
legend("center",legend=levels(sams),pch=points,col=colors,ncol=4)
dev.off()

##Saving counts and metadata
d[[1]] -> counts_Vd2_2
meta -> meta_Vd2_2

write.table(counts_Vd2_2, file = "output_data/counts_Age_Vd2.txt", sep = "\t",
            row.names = TRUE)
write.table(meta_Vd2_2, file = "output_data/metadata_Age_Vd2.txt", sep = "\t",
            row.names = TRUE)


######Monocyte data######

###Quality check
counts <- counts_monos
meta <- meta_monos

##Create DGEList object
d0 <- DGEList(counts)

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

#A mean-difference plot (MD-plot) is a plot of log-intensity ratios (differences) 
#versus log-intensity averages (means). Useful to detect outliers.

for (i in 1:20) {
  pdf(file= paste("graphs/MD_Monocytes_",i,".pdf", sep = ""))
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

pdf(file="graphs/MDS_Monocytes_pergroup.pdf")
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19) 
colors <- rep(c("blue","darkgreen","red", "purple"),2)
plotMDS(cpm(d[[1]]), col=colors[group], pch=points[group], xlab = "leading FC dim1", ylab = "leading FC dim2") 
legend("topleft",legend=levels(group),pch=points,col=colors,ncol=2)
dev.off()

##Plot with samples ID
sam <- as.data.frame(colnames(counts))
sam %>% extract(col = `colnames(counts)`, into = "samples", regex = "(RDH\\d+)", remove = FALSE) -> sam
sam %>% extract(col = `colnames(counts)`, into = "state", regex = "(Exviv|4hpRB\\w)", remove = FALSE) -> sam
sam$state <- gsub("Exviv", "unstim", sam$state)
sam$state <- gsub("4hpRBC", "stim", sam$state)
sam$plotid = paste0(sam$samples, "_", sam$state)
head(sam)

plot_counts <- counts
colnames(plot_counts) <- sam$plotid
sams <- factor(sam$plotid)

pdf(file="graphs/MDS_Monocytes_persample.pdf")
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19) 
colors <- rep(c("blue","darkgreen","red","purple","orange","brown","salmon","black","turquoise","darkgrey"),2)
plotMDS(cpm(d[[1]]), col=colors[sams], pch=points[sams], xlab = "leading FC dim1", ylab = "leading FC dim2") 
legend("center",legend=levels(sams),pch=points,col=colors,ncol=4)
dev.off()

##Saving counts and metadata
d[[1]] -> counts_monos
meta -> meta_monos

write.table(counts_monos, file = "output_data/counts_Age_monos.txt", sep = "\t",
            row.names = TRUE)
write.table(meta_monos, file = "output_data/metadata_Age_monos.txt", sep = "\t",
            row.names = TRUE)

