##Analysis of the Age-monocyte bulk transcriptional data
##PBMCs from n=5 malaria naive adults and n=5 malaria naive children were stimulated with parasites- 4h, 1:1
##unstimulated and stimulated monocytes cells were sorted and sequenced
##data was analysed with glm to take into account repeated measures and stimulation


#library(readxl)
library(ggplot2)
library(tidyr)
#library(plyr)
library(dplyr)
#library(ggpubr)
#library(reshape2)
#library(gridExtra)
#library(corrplot)
#library(RColorBrewer)
#library(tidyverse)
#library(pheatmap)
#library(scales)
#library(amap)
#library(Hmisc)
#library(corrplot)
#library(grid)
#library(gridExtra)
#library(umap)
#library(Hmisc)
#library(viridis)
#library(viridisLite)
#library(pgirmess)
#library(ggiraph)
#library(ggiraphExtra)
#library(FSA)
#library(useful)
#library(factoextra)
#library(ggforce)
#library(concaveman)
#library(rstatix)
#library(factoextra)
#library(FactoMineR)
#library(multipanelfigure)

##here I am using this data to understand the differences between malaria naive children and adults after pRBC stim. 
##I downloaded data from L drive on the 25September23

setwd("F:/2024/Age/BulkRNAseq/monocytes")

##start with the mono data

glmmseq_age_mono <- read.csv("age.monos.2024.csv")

glmmseq_age_mono <- read.table(here("output_data/glmmseq_Age_Monocytes_output.txt"), header = TRUE)

str(glmmseq_age_mono)

glmmseq_age_mono$Sig_state[glmmseq_age_mono$qvals.state_bool<=0.05] <- 1
glmmseq_age_mono$Sig_state[glmmseq_age_mono$qvals.state_bool>0.05] <- 0

glmmseq_age_mono$Sig_age[glmmseq_age_mono$qvals.age_bool<=0.05] <- 1
glmmseq_age_mono$Sig_age[glmmseq_age_mono$qvals.age_bool>0.05] <- 0

glmmseq_age_mono$Sig_state.age[glmmseq_age_mono$qvals.state_bool.age_bool<=0.05] <- 1
glmmseq_age_mono$Sig_state.age[glmmseq_age_mono$qvals.state_bool.age_bool>0.05] <- 0


glmmseq_age_mono$Sig_cat[glmmseq_age_mono$Sig_state == 1 & glmmseq_age_mono$Sig_age == 0 & glmmseq_age_mono$Sig_state.age==0 ] <- "state_only"
glmmseq_age_mono$Sig_cat[glmmseq_age_mono$Sig_state == 0 & glmmseq_age_mono$Sig_age == 1 & glmmseq_age_mono$Sig_state.age==0 ] <- "age_only"
glmmseq_age_mono$Sig_cat[glmmseq_age_mono$Sig_state == 0 & glmmseq_age_mono$Sig_age == 0 & glmmseq_age_mono$Sig_state.age==1 ] <- "interaction_only"
  
glmmseq_age_mono$Sig_cat[glmmseq_age_mono$Sig_state == 1 & glmmseq_age_mono$Sig_age == 1 & glmmseq_age_mono$Sig_state.age==0 ] <- "stateANDage"
glmmseq_age_mono$Sig_cat[glmmseq_age_mono$Sig_state == 1 & glmmseq_age_mono$Sig_age == 0 & glmmseq_age_mono$Sig_state.age==1 ] <- "stateANDinteraction"
glmmseq_age_mono$Sig_cat[glmmseq_age_mono$Sig_state == 0 & glmmseq_age_mono$Sig_age == 1 & glmmseq_age_mono$Sig_state.age==1 ] <- "ageANDinteraction"

glmmseq_age_mono$Sig_cat[glmmseq_age_mono$Sig_state == 1 & glmmseq_age_mono$Sig_age == 1 & glmmseq_age_mono$Sig_state.age==1 ] <- "ALL"


glmmseq_age_mono$Sig_cat[glmmseq_age_mono$Sig_state ==0 & glmmseq_age_mono$Sig_age == 0 & glmmseq_age_mono$Sig_state.age==0 ] <- "none"

table(glmmseq_age_mono$Sig_cat)

#new analysis-- will clean up script
#age_only   ageANDinteraction                 ALL    interaction_only                none          state_only         stateANDage 
#167                  35                      138               559                 4711                7660                 464 
#stateANDinteraction 
#1539 


#############################################################################
#####scatterplot to describe these data###################################################

##note- in this new analysis- we are including STATEandAGE genes as interaction genes.

##think about the interesting genes- removed the none genes
glmmseq_age_mono_sig <- glmmseq_age_mono %>%
  subset(Sig_cat!= "none")

glmmseq_age_mono_sig$Sig_cat2 <- glmmseq_age_mono_sig$Sig_cat
glmmseq_age_mono_sig$Sig_cat2[glmmseq_age_mono_sig$Sig_cat=="interaction_only"] <- "interaction" 
glmmseq_age_mono_sig$Sig_cat2[glmmseq_age_mono_sig$Sig_cat=="stateANDinteraction"] <- "interaction" 
glmmseq_age_mono_sig$Sig_cat2[glmmseq_age_mono_sig$Sig_cat=="ageANDinteraction"] <- "interaction" 
glmmseq_age_mono_sig$Sig_cat2[glmmseq_age_mono_sig$Sig_cat=="ALL"] <- "interaction"
glmmseq_age_mono_sig$Sig_cat2[glmmseq_age_mono_sig$Sig_cat=="stateANDage"] <- "interaction"

table(glmmseq_age_mono_sig$Sig_cat2)
#age_only  state_only interaction 
#167        7660        2735 

glmmseq_age_mono_sig$Sig_cat2 <- factor(glmmseq_age_mono_sig$Sig_cat2, levels= c("age_only", "state_only", "interaction"))


custom.colours <- c("Adult" = '#3399FF', "Child" = '#FF3333')
custom.colours2 <- c("age_only" = '#CCCCCC', "interaction" = '#333333', "state_only" ='#999999')



#Adding values to the quadrants
quad_count <- glmmseq_age_mono_sig %>%
  count(Sig_cat2, right = State_AdultFC > 0, top = State_ChildFC > 0) %>%
  mutate(State_AdultFC = 10 * (right - 0.5), State_ChildFC = 10 * (top - 0.5))

##Data by State
pdf("Scatter_Stim_numbers_2024.pdf", h=8.25 , w=11.75)
g1 <- ggplot(glmmseq_age_mono_sig, aes(State_AdultFC, State_ChildFC))+
  geom_point(aes(colour=Sig_cat2))+
  facet_grid(.~Sig_cat2)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  ylab("Delta STIM Child")+
  xlab("Delta STIM Adult")+
  scale_colour_manual(values = custom.colours2) +
  geom_text(data = quad_count, aes(label = n, x = State_AdultFC, y = State_ChildFC), size = 4) + 
  theme_bw()
print(g1)
dev.off()

##Data by Age
quad_count2 <- glmmseq_age_mono_sig %>%
  count(Sig_cat2, right = Age_unstim > 0, top = Age_stim > 0) %>%
  mutate(Age_unstim = 5 * (right - 0.5), Age_stim = 5 * (top - 0.5))

pdf("Scatter_Age_numbers_2024.pdf", h=8.25 , w=11.75)
g2 <- ggplot(glmmseq_age_mono_sig, aes(Age_unstim, Age_stim))+
  geom_point(aes(colour=Sig_cat2))+
  facet_grid(.~Sig_cat2)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  ylab("Delta Age Stim")+
  xlab("Delta Age unstim")+
  scale_colour_manual(values = custom.colours2) +
  geom_text(data = quad_count2, aes(label = n, x = Age_unstim, y = Age_stim), size = 4) +
  theme_bw()
print(g2)
dev.off()


#write this files to CSV so dont have to run the above significance code again.
write.csv(glmmseq_age_mono_sig, "monos.sig.2024.csv")

#read in CSV file
glmmseq_age_mono_sig <- read.csv("monos.sig.2024.csv")

#First will look at baseline responses- so this is genes which were significant for age

####################################################################
#so now I want to subset based on genes which are significant for age
#these include ALL, age_only, age and interaction, interaction only 

####glmmseq_age_mono <- glmmseq_age_mono.clean #What is glmseq_age_mono.clean?####
glmmseq_age_mono.clean <- glmmseq_age_mono_sig %>%
                        filter(Sig_cat == "ALL"|
                               Sig_cat == "age_only"|
                               Sig_cat == "ageANDinteraction"|
                                 Sig_cat == "interaction_only")
glmmseq_age_mono <- glmmseq_age_mono.clean
##how to think about this data-- selecting genes which are significant for age 
glmmseq_age_mono$Age <- NA

glmmseq_age_mono$Age[glmmseq_age_mono$Sig_cat=="age_only"] <- "Age" 
glmmseq_age_mono$Age[glmmseq_age_mono$Sig_cat=="ageANDinteraction"] <- "Age" 
glmmseq_age_mono$Age[glmmseq_age_mono$Sig_cat=="ALL"] <- "Age" 

sig.age <- subset(glmmseq_age_mono, Age == "Age")
dim(sig.age)
####I got 340 genes? Zuly 18112024####

#320 genes which are significant for Age.
write.csv(sig.age, "sig.age.csv")

#selecting genes based on if they are higher in adults or children at unstim
sig.age$Delta.higher.unstim[sig.age$Age_unstim > 0] <- "Adult"
sig.age$Delta.higher.unstim[sig.age$Age_unstim < 0] <- "child"

################################################################################################
##Plots of data
##thinking we could also do a scatter plot of the age-significant genes- as it includes not only the "age-only" group
#only scatter by age because it makes sense
custom.colours2 <- c("age_only" = '#999999', "interaction" = '#333333', "state_only" = '#000000')

##Data by Age-- adding values to quadrants
quad_count2 <- sig.age %>%
  count(Sig_cat2, right = Age_unstim > 0, top = Age_stim > 0) %>%
  mutate(Age_unstim = 5 * (right - 0.5), Age_stim = 5 * (top - 0.5))

pdf("Scatter_sig.age.pdf", h=8.25 , w=11.75)
g2 <- ggplot(sig.age, aes(Age_unstim, Age_stim))+
  geom_point(aes(colour=Sig_cat2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  ylab("Delta Age Stim")+
  xlab("Delta Age unstim")+
  scale_colour_manual(values = custom.colours2) +
  geom_text(data = quad_count2, aes(label = n, x = Age_unstim, y = Age_stim), size = 4) +
  theme_bw()
print(g2)
dev.off()

## just make a simple bar plot to show the number of genes which are higher at baseline between children and adults.
table(sig.age$Delta.higher.unstim)
#Adult Child 
#241   79
colnames(sig.age)
barplot1 <- sig.age[,c(13,51)]

#sum up all the genes in each group
colSums(table(barplot1))
#make it into a dataframe
barplot1 <- as.data.frame(colSums(table(barplot1)))
#rename column
names(barplot1)[names(barplot1) == "colSums(table(barplot1))"] <- "DEG.number"
barplot1$Group = rownames(barplot1)

barplot1$Group <- factor(barplot1$Group, levels = c("child", "Adult"))

summary.bar <- ggplot(barplot1, aes(x = Group, y = DEG.number)) +
  geom_bar(stat = "identity", fill = "#333333")+
  theme_bw() +
  geom_text(aes(label = DEG.number), hjust = 0.5, size = 6, nudge_y = 5) +
  labs(title = "number of DEGs higher in each group at baseline", x = "Group", y = "DEG.number")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
summary.bar

#for the age significant group
#not sure if this will be interesting but potentially look at the Age_FC across genes at baseline and stim.
#what it shows is that the majority of genes which are significant for age have higher gene expression in adults before and after stim.

age.violin <- sig.age %>%
  dplyr::select(SYMBOL, qvals.state_bool, qvals.age_bool, qvals.state_bool.age_bool, State_AdultFC, State_ChildFC, Age_stim, Age_unstim, Sig_cat, Age) %>%
  dplyr::rename(unstim = Age_unstim, stim = Age_stim)

##make long
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)

age.violin.long <- gather(age.violin, 
                          stim, 
                          age_FC,
                          unstim:stim, 
                          factor_key = TRUE)

library(rstatix)
library(reshape2)
library(ggpubr)

boxplot.colour <- c("unstim" = "#D1D3D4", "stim" = "#58595B")
my_comparisons <- list(c("Adult", "child"))

age.violin <- ggplot(age.violin.long, aes(x = stim, y = age_FC, fill = stim)) +
  geom_violin() +
  scale_fill_manual(values=c('#989898', '#333333'))+
  geom_boxplot(width = 0.2, fill = boxplot.colour, color = "black") +
  labs(title = "distribution of genes at baseline and after stim", x = "Age", y = "Age_FC") +
  theme_bw()
age.violin

#I think we present these data as a bar plot now.
# Sort the dataframe based on FDR in descending order

#need to just select unstim data
adult.baseline <- subset(age.violin.long, stim == "unstim")

order.BH <- adult.baseline%>% arrange(desc(qvals.age_bool))

#selecting the top 100 genes
subset.long<- order.BH[271:320, ]

#set a FC threshold
qval.threshold <- 0.01

# Subset the DataFrame based on the qvals.age_bool
qval_0.01 <- subset(order.BH, qvals.age_bool <= qval.threshold)

## This tells us if the gene was higher in adults (positive) or higher in children (negative) before and after stim
g6 <- ggplot(subset.long, aes (x=reorder(SYMBOL, age_FC), y=age_FC, fill=stim))+
  geom_bar(stat="identity", colour= "black", position =position_dodge(), width=0.7)+
  scale_fill_manual(values=c('#99999999', "#E69F00"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  labs(title = "Monocytes Baseline gene expression top50", x = "Gene symbol", y = "Age_FC")
print(g6)

##IPA pathway plots-- Age differences at baseline
setwd("F:/2023/Age_bulk/IPA/Baseline")

##start with the mono data

B.pathways<- read.csv("Baseline.pathways.csv")

#remove NAs
B.pathways_noNA <- B.pathways%>%
  filter(z.score!= "#NUM!")

#make numeric values
B.pathways_noNA$z.score<- as.numeric(B.pathways_noNA$z.score)
B.pathways_noNA$pvalue<- as.numeric(B.pathways_noNA$pvalue)

# Sort the dataframe based on FDR in descending order
order.pvalue <- B.pathways_noNA %>% arrange(desc(pvalue))
#selecting the significant pathways
threshold <- 1.3
baseline1.3<- subset(order.pvalue, pvalue > threshold)

## This tells us if the gene was higher in adults (positive) or higher in children (negative) before and after stim
g6 <- ggplot(baseline1.3, aes (x=reorder(Canonical.Pathways, z.score ), y=z.score))+
  geom_bar(stat="identity", colour= "black", position =position_dodge(), width=0.7)+
  scale_fill_manual(values=c("#010101"))+
  theme_bw()+
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  labs(title = "Baseline pathways-significant", x = "canonical pathway", y = "z.score")
print(g6)

#now do the same for the upstream regulators
B.upstream<- read.csv("Baseline.upstream.csv")
B.upstream_noNA <- B.upstream%>%
  filter(z.score!= "NA")

#selected which have a predicted activation state
B.state <- B.upstream_noNA [(B.upstream_noNA$Predicted.State %in% c("Activated", "Inhibited")), ]

g7 <- ggplot(B.state, aes (x=reorder(Upstream.Regulator, z.score), y=z.score))+
  geom_bar(stat="identity", colour= "black", position =position_dodge(), width=0.7)+
  scale_fill_manual(values=c("#010101"))+
  theme_bw()+
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  labs(title = "Baseline upstream regulators", x = "upstream regulators", y = "z.score")
print(g7)
####################################################################################################################################
####################################################################################################################################
#########################################################################################################
##Now we want to look at interaction genes-- so these will make up most of the analysis in the main figure

table(glmmseq_age_mono_sig$Sig_cat2)
#age_only  state_only interaction 
#167        7660        2735            

# 
##this tells us if the slope is of the change in gene expression is bigger or smaller in children and adults
glmmseq_age_mono_sig$stim_FC_bigger_in <- ifelse(glmmseq_age_mono_sig$State_AdultFC > glmmseq_age_mono_sig$State_ChildFC, "Adult", "Child")

#now I want to create a column to look at which genes are higher in adults or children at baseline (Age_unstim)
glmmseq_age_mono_sig$unstim_bigger_in <- ifelse(glmmseq_age_mono_sig$Age_unstim > 0, "Adult", "Child")

#what about genes which are higher in adults and children after stim (Age_stim)
glmmseq_age_mono_sig$stim_bigger_in <- ifelse(glmmseq_age_mono_sig$Age_stim > 0, "Adult", "Child")

#subset for only interaction genes
interaction <- subset(glmmseq_age_mono_sig, Sig_cat2 == "interaction")

#this tells me that the majority of genes have a greater change in adults-- the slope is bigger in adults
table(interaction$stim_FC_bigger_in)
#Adult Child 
#2184   551 

#Overall I want to see the number of genes that go down or up in the INTERACTION genes
#the sope of the change is bigger in adults but this will now tell me the direction of that change
interaction$stim_goes[interaction$State_AdultFC >0 & interaction$State_ChildFC > 0] <- "Up"
interaction$stim_goes[interaction$State_AdultFC <0 & interaction$State_ChildFC < 0] <- "Down"
interaction$stim_goes[interaction$State_AdultFC >0 & interaction$State_ChildFC < 0] <- "Up.adult"
interaction$stim_goes[interaction$State_AdultFC <0 & interaction$State_ChildFC > 0] <- "Up.child"

#this is telling me that within the interaction genes, how many go down, up and opposite directions for children and adults
table(interaction$stim_goes)
#Down       Up     Up.adult Up.child 
#897        814      896      128 

#this tells us the number of genes which start higher in children or adults
table(interaction$unstim_bigger_in)
#Adult Child 
#699  2036

#this tells us the number of genes which end higher in children or adults-AFTER STIM
table(interaction$stim_bigger_in)
#Adult Child 
#2490   254

#this breaks down the direction of those changes between children and adults.
table(interaction$unstim_bigger_in, interaction$stim_bigger_in, interaction$stim_goes)
#Down
          #stim.Adult stim.Child
#unstim.Adult    367      61
#unstim.Child   448       21

#Up
          #stim.Adult stim.Child
#unstim.Adult    123        15
#unstim.Child   619         57

#Up.adults
         #stim.Adult stim.Child
#unstim.Adult    6         0
#unstim.Child   876        14

#Up.child
         #stim.Adult stim.Child
#unstim.Adult    51        76
#unstim.Child    0         1

#############################################################################################################
#here- we are trying to understand these data more-- looking at genes which are always up in adult or children
#now I am going to subset all the genes which are MORE up in adults when compared to children.
#For adults- Up.up (19), down.up (603), Up. Adult.up.up (5), Adult.down.up (857). Total up genes- 1484 genes -- ((1484)/2219)*100 = 67%

Up.adult <-interaction[!(interaction$stim_goes %in% c("Down", "Up.child")), ]
#subset further to choose the genes which are up in adults after stim
#these are the genes to look at in IPA/STRING-
Up.adult2 <-Up.adult[!(Up.adult$stim_bigger_in %in% c("Child")), ]

#export this list of genes to CSV
write.csv(Up.adult2, "monos.adult.up.2024.csv")

#now I am going to subset all the genes which are MORE up in children when compared to adults.
Up.child <-interaction[!(interaction$stim_goes %in% c("Down", "Up.adult")), ]
#these are the genes to look at in IPA/STRING-- or futherdown...
Up.child2 <-Up.child[!(Up.child$stim_bigger_in %in% c("Adult")), ]

write.csv(Up.child2, "monos.child.up.2024.csv")

#####################################################################################
#not really doing the down genes
#subset for genes which are down in adults
Down.adult <-interaction[!(interaction$stim_goes %in% c("Up", "Up.adult")), ]
#subset for genes which are down in children
Down.child <-interaction[!(interaction$stim_goes %in% c("Up", "Up.child")), ]
####################################################################################

setwd("F:/2024/Age/BulkRNAseq/monocytes")
adult.up <- read.csv("monos.adult.up.2024.csv")

#too many genes in the adults to visulise- so make an FC cutoff
#set a FC threshold
#different options to set a cut off
state.threshold <- 1.5
qval.threshold <- 0.001
cutoff <- 3.450000e-295
FC.pval <- subset(FC.1.5, qvals.state_bool.age_bool < qval.threshold)

# Subset the DataFrame based on the qvals.age_bool

FC.1.5 <- subset(Up.adult2, State_AdultFC > state.threshold)
#73 genes which with a P-value <0.001

adult.sig <- subset(Up.adult2, qvals.state_bool < qval.threshold)
#82 observations, note there are some NAs
#na.omit

adult.sig <- subset(Up.adult2, qvals.state_bool.age_bool < qval.threshold)

child.sig <- subset(Up.child2, qvals.state_bool.age_bool < qval.threshold)
dim(child.sig)
#82 observations, note there are some NAs
#na.omit

child.sig.noNA <- na.omit(child.sig)

#create a threhold
FC.1.5.Child <- subset(Up.child2, State_ChildFC > state.threshold)
dim(FC.1.5.Child)
#29 observations

##BAR GRAPHS OF THESE

age <- child.sig.noNA  %>%
  select(SYMBOL, qvals.state_bool, qvals.age_bool, qvals.state_bool.age_bool, State_AdultFC, State_ChildFC, Age_stim, Age_unstim, Sig_cat, Sig_cat2) %>%
  dplyr::rename(unstim = Age_unstim, stim = Age_stim)

state <- child.sig.noNA %>%
  select(SYMBOL, qvals.state_bool, qvals.age_bool, qvals.state_bool.age_bool, State_AdultFC, State_ChildFC,  Age_stim, Age_unstim, Sig_cat, Sig_cat2) %>%
  dplyr::rename(Adult = State_AdultFC, Child = State_ChildFC)

##make long
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)

age.long <- gather(age, 
                   stim, 
                   age_FC,
                   unstim:stim, 
                   factor_key = TRUE)

state.long <- gather(state, 
                     Age, 
                     State_FC,
                     Adult:Child, 
                     factor_key = TRUE)


## genes that are only differnt between age (but not stim)
g6 <- ggplot(age.long, aes (x=reorder(SYMBOL, age_FC), y=age_FC, fill=stim))+
  geom_bar(stat="identity", colour= "black", position =position_dodge(), width=0.7)+
  scale_fill_manual(values=c('#99999999', "#E69F00"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  xlab("FC (Adult / Child)")
print(g6)

g6.1 <- ggplot(state.long, aes(x=(reorder(SYMBOL, ((Age_stim+Age_unstim)/2))), y=State_FC, fill=Age))+
  geom_bar(stat="identity", colour= "black", position =position_dodge(), width=0.7)+
  scale_fill_manual(values=c('lightblue', 'red'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  xlab("FC (Stim / Unstim)")
print(g6.1)

###################################################################################################################################
############################################################################################################

##Used Zuly's code to make the barplot below

#going to make a bar plot to summarise these interaction genes in their groups
#I found this quite a work around to get the dataframe into a format which was capable of makeing this barpot
#creating a column grouping the different variables of interest
interaction$dir <- NA
interaction$dir <- paste0(interaction$unstim_bigger_in,"_",interaction$stim_bigger_in,"_",interaction$stim_goes)
colnames(interaction)

#here I took the ENSMBL number instead of the gene name
barplot1 <- interaction[,c(1,23)]

#sum up all the genes in each group
colSums(table(barplot1))
#make it into a dataframe
barplot1 <- as.data.frame(colSums(table(barplot1)))
#rename column
names(barplot1)[names(barplot1) == "colSums(table(barplot1))"] <- "DEG.number"
barplot1$Group = rownames(barplot1)



summary.bar <- ggplot(barplot1, aes(x = reorder(Group, -DEG.number), y = DEG.number)) +
  geom_bar(stat = "identity", fill = "#333333")+
  theme_bw() +
  geom_text(aes(label = DEG.number), hjust = 0.5, size = 6, nudge_y = 40) +
  labs(title = "interaction.genes_MonoAgeBulk", x = "Group", y = "DEG.number")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
summary.bar

write.csv(barplot1, "monocytes.barplot.csv")

###################################################################################
#using zuly's code for the barplot

monos_filter <- interaction

#Creating the groups for Barplot
monos_filter$groupsAge<-NA
monos_filter$groupsAge[monos_filter$Age_unstim<0 & monos_filter$Age_stim<0 & monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC<0] <- "A"#A
monos_filter$groupsAge[monos_filter$Age_unstim<0 & monos_filter$Age_stim<0 & monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC>0] <- "B"#B
monos_filter$groupsAge[monos_filter$Age_unstim<0 & monos_filter$Age_stim<0 & monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC<0] <- "C"#C
monos_filter$groupsAge[monos_filter$Age_unstim<0 & monos_filter$Age_stim<0 & monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC>0] <- "D"#D
monos_filter$groupsAge[monos_filter$Age_unstim>0 & monos_filter$Age_stim<0 & monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC>0] <- "E"#E
monos_filter$groupsAge[monos_filter$Age_unstim>0 & monos_filter$Age_stim<0 & monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC<0] <- "F"#F
monos_filter$groupsAge[monos_filter$Age_unstim>0 & monos_filter$Age_stim<0 & monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC<0] <- "G"#G
monos_filter$groupsAge[monos_filter$Age_unstim>0 & monos_filter$Age_stim>0 & monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC<0] <- "H"#H
monos_filter$groupsAge[monos_filter$Age_unstim>0 & monos_filter$Age_stim>0 & monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC<0] <- "I"#I
monos_filter$groupsAge[monos_filter$Age_unstim<0 & monos_filter$Age_stim>0 & monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC<0] <- "J"#J
monos_filter$groupsAge[monos_filter$Age_unstim<0 & monos_filter$Age_stim>0 & monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC>0] <- "K"#K
monos_filter$groupsAge[monos_filter$Age_unstim<0 & monos_filter$Age_stim>0 & monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC>0] <- "L"#L
monos_filter$groupsAge[monos_filter$Age_unstim>0 & monos_filter$Age_stim>0 & monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC>0] <- "M"#M
monos_filter$groupsAge[monos_filter$Age_unstim>0 & monos_filter$Age_stim>0 & monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC>0] <- "N"#N
table(monos_filter$groupsAge)
 #A   B   C   D   E   F   G   H   I   J   K   L   M   N 
#21  14   1  57  15  76  61  51 367 448 876 619 123   6 

##Creating direction groups
monos_filter$Direction_Ch_Ad<-NA
monos_filter$Direction_Ch_Ad[monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC>0] <- "UP_UP"
monos_filter$Direction_Ch_Ad[monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC<0] <- "DOWN_DOWN"
monos_filter$Direction_Ch_Ad[monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC>0] <- "DOWN_UP"
monos_filter$Direction_Ch_Ad[monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC<0] <- "UP_DOWN"
table(monos_filter$Direction_Ch_Ad)

#DOWN_DOWN   DOWN_UP   UP_DOWN     UP_UP 
#   897       896       128       814 

##Creating variable higher expression pre-stim
monos_filter$HigherExp_prestim<-NA
monos_filter$HigherExp_prestim[monos_filter$Age_unstim<0] <- "Children"
monos_filter$HigherExp_prestim[monos_filter$Age_unstim>0] <- "Adults"
table(monos_filter$HigherExp_prestim)

#Adults Children 
#699     2036 

##Creating variable higher expression post-stim
monos_filter$HigherExp_poststim<-NA
monos_filter$HigherExp_poststim[monos_filter$Age_stim<0] <- "Children"
monos_filter$HigherExp_poststim[monos_filter$Age_stim>0] <- "Adults"
table(monos_filter$HigherExp_poststim)

#Adults Children 
#2490      245 


##Suggested barplot code for monocytes

monos_bardata <- as.data.frame(table(monos_filter$groupsAge,monos_filter$HigherExp_prestim, monos_filter$HigherExp_poststim,monos_filter$Direction_Ch_Ad))
monos_bardata1 <- monos_bardata[!monos_bardata$Freq==0,]
head(monos_bardata1)

colnames(monos_bardata1)[1]="Group"
colnames(monos_bardata1)[2]="HigherExp_prestim"
colnames(monos_bardata1)[3]="HigherExp_poststim"
colnames(monos_bardata1)[4]="Direction_Ch_Ad"
colnames(monos_bardata1)[5]="No of DEGs"
head(monos_bardata1)

monos_bardata1$Group <- factor(monos_bardata1$Group, levels=c("A", "B", "C", "D","E", "F", "G", "H", "I", "J", "K", "L", "M", "N"))

monos_bardata2 <- monos_bardata1 %>% arrange(Group)
head(monos_bardata2)

##Adding coloring info
monos_bardata2%>%
  mutate(
    GeneSelect= case_when(HigherExp_poststim=="Children"& Direction_Ch_Ad=="UP_UP"|HigherExp_poststim=="Children"&Direction_Ch_Ad=="UP_DOWN"~"Child",
                          HigherExp_poststim=="Adults"& Direction_Ch_Ad=="UP_UP"|HigherExp_poststim=="Adults"& Direction_Ch_Ad=="DOWN_UP"~"Adult",
                          TRUE~"Not selected"))-> monos_bardata3


table(monos_bardata3$GeneSelect)
#Adult        Child Not selected 
#4            4            6 

# Barplot with table
#FF9999 = child
#lightblue =Adult
##For children we selected the genes that were upregulated in children
##Direction_Ch_Ad (UP_UP, UP_DOWN) and Exp_poststim was higher in children(Children)
##For children we selected the genes that were upregulated in children
##Direction_Ch_Ad (UP_UP, UP_DOWN) and Exp_poststim was higher in children(Children)

boxplot.colour <- c("Child"="#FF9999", "Not selected"="black","Adult"="lightblue")

##Simple barplot for the presentation
summary.bar <- ggplot(monos_bardata3, aes(x = Group, y =`No of DEGs`, fill = GeneSelect)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values=boxplot.colour)+
  geom_text(aes(label=`No of DEGs`), vjust=-0.5)+
  theme_bw() +
  labs(title = " ", x = "Group", y = "Interaction genes")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
summary.bar
##################################################################################################################################
#violin plots to show the distribution of the slope in gene expression change after pRBC stim- for UP and DOWN genes 
#subset for up genes

library(rstatix)
library(reshape2)
library(ggpubr)

Up <-interaction[!(interaction$stim_goes %in% c("Down", "Up.child", "Up.adult")), ]

table(Up$stim_FC_bigger_in)
#Adult Child 
#695   119 

#So this is telling me that the slope of the change in gene expression is greater in adults most of the time

#subset for down genes
Down <-interaction[!(interaction$stim_goes %in% c("Up", "Up.child", "Up.adult")), ]

#create a column to identify which have greater negative number.
Down$stim_FC_smaller_in[Down$State_AdultFC > Down$State_ChildFC] <- "Child"
Down$stim_FC_smaller_in[Down$State_AdultFC < Down$State_ChildFC] <- "Adult"

table(Down$stim_FC_smaller_in)
#Adult Child 
#304   593 

# this shows the majority of genes that are going down are also going further down in children compared to adults
# thinking that a violin plot might look good...

as.matrix(names(Up))

#change to dataframes labelled 'Up' or 'Down' to make the violin plots
state <- Up %>%
  select(SYMBOL, qvals.state_bool, qvals.age_bool, qvals.state_bool.age_bool, State_AdultFC, State_ChildFC,  Age_stim, Age_unstim, Sig_cat, Sig_cat2) %>%
  dplyr::rename(Adult = State_AdultFC, Child = State_ChildFC)

##make long
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)

state.long <- gather(state, 
                     Age, 
                     State_FC,
                     Adult:Child, 
                     factor_key = TRUE)


boxplot.colour <- c("Adult" = "light blue", "Child" = "#FF9999")
my_comparisons <- list(c("Adult", "Child"))
state.long$Age <- factor(state.long$Age, levels = c("Child", "Adult"))

Up.violin <- ggplot(state.long, aes(x = Age, y = State_FC, fill = Age)) +
  geom_violin() +
  scale_fill_manual(values=c('#3399FF', '#FF3333'))+
  stat_compare_means(comparisons = my_comparisons, paired=FALSE, method="wilcox",
                     symnum.args = list(cutpoints = c(0,0.01, 0.05, Inf), 
                                        symbols = c("p < 0.01","p < 0.05", "ns"))) +
  geom_boxplot(width = 0.2, fill = boxplot.colour, color = "black") +
  labs(title = "Slope of the change, Up genes", x = "Age", y = "State_FC") +
  theme_bw()
Up.violin

Down.violin <- ggplot(state.long, aes(x = Age, y = State_FC, fill = Age)) +
  geom_violin() +
  scale_fill_manual(values=c('#3399FF', '#FF3333')) +
  stat_compare_means(comparisons = my_comparisons, paired=FALSE, method="wilcox", 
                     symnum.args = list(cutpoints = c(0,0.01, 0.05, Inf),
                                        symbols = c("p < 0.01","p < 0.05", "ns"))) +
  geom_boxplot(width = 0.2, fill = boxplot.colour, color = "black") +
  labs(title = "Slope of the change, Down genes", x = "Age", y = "State_FC") +
  theme_bw()
Down.violin

###########################################################################################################################
#IPA analysis of interaction genes.
#ADULTS:
#Genes which went UP after stim and were higher in adults overall after stim were input into IPA. We used the STATE_AGE_qval and Adult_STATEFC as the 
#log ratio. With the same genes, we then used the Child_STATEFC to compare the adult and children pathways enriched.

#CHILDREN
#Genes which went UP after stim and were higher in children overall after stim were input into IPA. We used the STATE_AGE_qval and Child_STATEFC as the 
#log ratio. With the same genes, we then used the Adult_STATEFC to compare the children and adults pathways enriched.
#note- only 9 pathways enriched for this analysis

library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(tidyr)

#set my working environment
#these might change depending on the comparison and where the CSV is located.
#I will make sure a copy of the CSV is saved in the same folder as this R-script and associated R-Markdown files

#PATHWAYS-IPA

setwd("F:/2024/Age/BulkRNAseq/monocytes/IPA/Adult.up")

df <-read.csv("Monos.AdultUp.pathways.csv")

##remember you might have to make values numeric

##make dataframe long--- this dataframe is for the POSITIVE z-scores
#gathers z.scores into one column and adds a new column for age

df.long <- pivot_longer(df, cols=c(2:3), names_to = "age", values_to= "z.score")
#gathers FDRs into one column and adds a new column for age
df.long2 <- pivot_longer(df, cols=c(4:5), names_to = "age", values_to= "BH-Pvalue")
#changes the names of the group to just Adult and Child
df.long$age <- ifelse(df.long$age == "Adult.zscore", "Adult", "Child")
df.long2$age <- ifelse(df.long2$age == "Adult.BH", "Adult", "Child")
#selects the columns we want to join together
test <- df.long[,c(1,4:5)]
test2 <- df.long2[,c(1,4:5)]
#joins the columns together
final_long <- left_join(test, test2)


#remove NAs
pathways_noNA <- final_long%>%
  filter(z.score!= "N/A")

#make numeric values
pathways_noNA$z.score<- as.numeric(pathways_noNA$z.score)
pathways_noNA$`BH-Pvalue`<- as.numeric(pathways_noNA$`BH-Pvalue`)

# Sort the dataframe based on FDR in descending order
order.BH <- pathways_noNA %>% arrange(desc(`BH-Pvalue`))

#selecting the top 60 pathways-- THIS IS FOR THE ADULT PATHWAYS
subset.long<- order.BH[1:120, ]

subset.long$ordering <- ifelse(subset.long$age == "Adult", subset.long$z.score, NA)
subset.long <-subset.long %>%
  arrange(ordering, Canonical.Pathways) %>%
  mutate(Canonical.Pathways = factor(Canonical.Pathways, levels = unique(Canonical.Pathways)))

#I think we present these data as a bar plot now.
pdf("Child.higher.pathways.ordered.pdf", h=8.25 , w=11.75)
graph1 <- ggplot(subset.long, aes(x = Canonical.Pathways, y= z.score, fill=age)) +
  geom_col(position = "dodge", width = 0.7)+
  scale_fill_manual(values=c("#5690CC", "#EF3D39")) +
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "z-score", x= "Canonical pathways", title = "Genes which are always higher in children after stim")
graph1
dev.off()


#UPSTREAM REGULATORS-IPA
#ADULTS- Upstream regulators split into "cytokines, transcription factors and transmembrane proteins"

## Cytokines- upstream- adults ##############################################################################################
#set my working environment
#updating this code so we only include significant upstream regulators
setwd("F:/2024/Age/BulkRNAseq/monocytes/IPA/Adult.up")

Cyto.list <- read.csv("Cytokine.table.csv")
Cyto.comp <- read.csv("monos.AdultUp.upstream.csv")

matching_cyto <- Cyto.comp[Cyto.comp$Upstream.Regulator %in% Cyto.list$Upstream.Regulator, ]

#remove NAs
matching_cyto_noNA <- matching_cyto %>%
  filter(Child.zscore!= "N/A" & Adult.zscore!= "N/A")

##make dataframe long--- this dataframe is for the POSITIVE z-scores
#gathers z.scores into one column and adds a new column for age
cyto.long <- pivot_longer(matching_cyto_noNA, cols=c(2:3), names_to = "age", values_to= "z.score")
#gathers FDRs into one column and adds a new column for age
cyto.long2 <- pivot_longer(matching_cyto_noNA, cols=c(4:5), names_to = "age", values_to= "BH-Pvalue")
#changes the names of the group to just Adult and Child
cyto.long$age <- ifelse(cyto.long$age == "Adult.zscore", "Adult", "Child")
cyto.long2$age <- ifelse(cyto.long2$age == "Adult.BH", "Adult", "Child")
#selects the columns we want to join together
test <- cyto.long[,c(1,4:5)]
test2 <- cyto.long2[,c(1,4:5)]
#joins the columns together
final_long <- left_join(test, test2)

#make numeric values
final_long$z.score<- as.numeric(final_long$z.score)
final_long$`BH-Pvalue`<- as.numeric(final_long$`BH-Pvalue`)

# Sort the dataframe based on FDR in descending order
order.BH <- final_long %>% arrange(desc(`BH-Pvalue`))

order.BH$ordering <- ifelse(order.BH$age == "Adult", order.BH$z.score, NA)
order.BH  <- order.BH  %>%
  arrange(ordering, Upstream.Regulators) %>%
  mutate(Upstream.Regulators = factor(Upstream.Regulators, levels = unique(Upstream.Regulators)))

#bargraph, pathways visulalisation
pdf("Adult.higher.cytokines.sig.2024.pdf", h=8.25 , w=11.75)
graph <-ggplot(order.BH, aes(x = Upstream.Regulators, y= z.score, fill=age)) + 
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values=c("#5690CC", "#EF3D39"))+
  scale_x_discrete(position = "bottom")+
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "Z-score", x= "Cytokines", title = "Adult Up genes: Upstream regulators- Cytokines")
graph
dev.off()

## Transcription factor comparison ############################################################################################

#set my working environment

TF.list <- read.csv("TF.table.csv")
TF.comp <- read.csv("monos.AdultUp.upstream.csv")

matching_TF <- TF.comp[TF.comp$Upstream.Regulator %in% TF.list$Upstream.Regulator, ]

#remove NAs
matching_TF_noNA <- matching_TF %>%
  filter(Child.zscore!= "N/A" & Adult.zscore!= "N/A")

##make dataframe long--- this dataframe is for the POSITIVE z-scores
#gathers z.scores into one column and adds a new column for age
TF.long <- pivot_longer(matching_TF_noNA, cols=c(2:3), names_to = "age", values_to= "z.score")
#gathers FDRs into one column and adds a new column for age
TF.long2 <- pivot_longer(matching_TF_noNA, cols=c(4:5), names_to = "age", values_to= "BH-Pvalue")
#changes the names of the group to just Adult and Child
TF.long$age <- ifelse(TF.long$age == "Adult.zscore", "Adult", "Child")
TF.long2$age <- ifelse(TF.long2$age == "Adult.BH", "Adult", "Child")
#selects the columns we want to join together
test <- TF.long[,c(1,4:5)]
test2 <- TF.long2[,c(1,4:5)]
#joins the columns together
final_long <- left_join(test, test2)

#make numeric values
final_long$z.score<- as.numeric(final_long$z.score)
final_long$`BH-Pvalue`<- as.numeric(final_long$`BH-Pvalue`)

# Sort the dataframe based on FDR in descending order
order.BH <- final_long %>% arrange(desc(`BH-Pvalue`))

order.BH$ordering <- ifelse(order.BH$age == "Adult", order.BH$z.score, NA)
order.BH  <- order.BH %>%
  arrange(ordering, Upstream.Regulators) %>%
  mutate(Upstream.Regulators = factor(Upstream.Regulators, levels = unique(Upstream.Regulators)))


#bargraph, pathways visulalisation
pdf("Adult.higher.TF.sig.2024.pdf", h=8.25 , w=11.75)
graph <-ggplot(order.BH, aes(x = Upstream.Regulators, y= z.score, fill=age)) + 
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values=c("#5690CC", "#EF3D39"))+
  scale_x_discrete(position = "bottom")+
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "Z-score", x= "TF", title = "Adult Up genes: Upstream regulators- Transcripton factors")
graph
dev.off()

#############################################################################################################
#only use this code if you are separating the plots between activated and inhibited upstream regulators
#filter positive z-score by adults-- "activated TFs"-- if making plots separately
matchingTF_activated <- matching_TF_noNA %>%
  filter(Child.zscore > 0) %>%
  group_by(Upstream.Regulators)

#filter positive z-score by adults-- "inhibited TFs"-- if making plots separately
matchingTF_down <- matching_TF_noNA %>%
  filter(Child.zscore < 0) %>%
  group_by(Upstream.Regulators)


##Transmembrane receptors comparison#####################################################################################

#set my working environment

TM.list <- read.csv("TM.table.csv")
TM.comp <- read.csv("monos.AdultUp.upstream.csv")

matching_TM <- TM.comp[TM.comp$Upstream.Regulator %in% TM.list$Upstream.Regulator, ]

#remove NAs
matching_TM_noNA <- matching_TM %>%
  filter(Child.zscore!= "N/A" & Adult.zscore!= "N/A")

##make dataframe long--- 
#gathers z.scores and BH P-values into one column and adds a new column for age
TM.long <- pivot_longer(matching_TM_noNA, cols=c(2:3), names_to = "age", values_to= "z.score")
#gathers FDRs into one column and adds a new column for age
TM.long2 <- pivot_longer(matching_TM_noNA, cols=c(4:5), names_to = "age", values_to= "BH-Pvalue")
#changes the names of the group to just Adult and Child
TM.long$age <- ifelse(TM.long$age == "Adult.zscore", "Adult", "Child")
TM.long2$age <- ifelse(TM.long2$age == "Adult.BH", "Adult", "Child")
#selects the columns we want to join together
test <- TM.long[,c(1,4:5)]
test2 <- TM.long2[,c(1,4:5)]
#joins the columns together
final_long <- left_join(test, test2)

#make numeric values
final_long$z.score<- as.numeric(final_long$z.score)
final_long$`BH-Pvalue`<- as.numeric(final_long$`BH-Pvalue`)

# Sort the dataframe based on FDR in descending order
order.BH <- final_long %>% arrange(desc(`BH-Pvalue`))

order.BH$ordering <- ifelse(order.BH$age == "Adult", order.BH$z.score, NA)
order.BH  <- order.BH  %>%
  arrange(ordering, Upstream.Regulators) %>%
  mutate(Upstream.Regulators = factor(Upstream.Regulators, levels = unique(Upstream.Regulators)))

#bargraph, pathways visulalisation
pdf("Adult.higher.TM.sig.2024.pdf", h=8.25 , w=11.75)
graph <-ggplot(order.BH, aes(x = Upstream.Regulators, y= z.score, fill=age)) + 
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values=c("#5690CC", "#EF3D39"))+
  scale_x_discrete(position = "bottom")+
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "Z-score", x= "Transmembrane proteins", title = "Upstream regulators- TM")
graph
dev.off()

##############################################################################################
###########################################################################################################################
#Children- pathways graph-- this will go into supplementary Figure 2

setwd("F:/2024/Age/BulkRNAseq/monocytes/IPA/Child.up")

df <-read.csv("Monos.ChildUp.pathways.sig.csv")

##remember you might have to make values numeric

##make dataframe long--- this dataframe is for the POSITIVE z-scores
#gathers z.scores into one column and adds a new column for age

df.long <- pivot_longer(df, cols=c(2:3), names_to = "age", values_to= "z.score")
#gathers FDRs into one column and adds a new column for age
df.long2 <- pivot_longer(df, cols=c(4:5), names_to = "age", values_to= "BH-Pvalue")
#changes the names of the group to just Adult and Child
df.long$age <- ifelse(df.long$age == "Adult.zscore", "Adult", "Child")
df.long2$age <- ifelse(df.long2$age == "Adult.BH", "Adult", "Child")
#selects the columns we want to join together
test <- df.long[,c(1,4:5)]
test2 <- df.long2[,c(1,4:5)]
#joins the columns together
final_long <- left_join(test, test2)


#remove NAs
pathways_noNA <- final_long%>%
  filter(z.score!= "N/A")

#make numeric values
pathways_noNA$z.score<- as.numeric(pathways_noNA$z.score)
pathways_noNA$`BH-Pvalue`<- as.numeric(pathways_noNA$`BH-Pvalue`)

# Sort the dataframe based on FDR in descending order
order.BH <- pathways_noNA %>% arrange(desc(`BH-Pvalue`))

order.BH$ordering <- ifelse(order.BH$age == "Child", order.BH$z.score, NA)
order.BH <-order.BH%>%
  arrange(ordering, Canonical.Pathways) %>%
  mutate(Canonical.Pathways = factor(Canonical.Pathways, levels = unique(Canonical.Pathways)))

#I think we present these data as a bar plot now.
pdf("Child.higher.pathways.ordered.pdf", h=8.25 , w=11.75)
graph1 <- ggplot(order.BH, aes(x = Canonical.Pathways, y= z.score, fill=age)) +
  geom_col(position = "dodge", width = 0.7)+
  scale_fill_manual(values=c("#5690CC", "#EF3D39")) +
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "z-score", x= "Canonical pathways", title = "Genes which are always higher in children after stim")
graph1
dev.off()

#CHILDREN- Upstream regulators- note there were that many significantly enriched, therefore we didnt have to split into groups.

setwd("F:/2024/Age/BulkRNAseq/monocytes/IPA/Child.up")

#CYTOKINES

#data is significantly enriched upstream regulators-comparing children and adults
Child.table <- read.csv("ChildUP.cytokine.csv")
Child.comp <- read.csv("Monos.ChildUp.upstream.csv")

matching_child <- Child.comp[Child.comp$Upstream.Regulators %in% Child.table$Upstream.Regulator, ]

matching_child_noNA <- matching_child %>%
   filter(Child.zscore!= "N/A" & Adult.zscore!= "N/A")

##make dataframe long--- 
#gathers z.scores/BH-P.values into one column and adds a new column for age
child.long <- pivot_longer(matching_child_noNA, cols=c(2:3), names_to = "age", values_to= "zscore")
#gathers FDRs into one column and adds a new column for age
child.long2 <- pivot_longer(matching_child_noNA, cols=c(4:5), names_to = "age", values_to= "BH_P.value")
#changes the names of the group to just Adult and Child
child.long$age <- ifelse(child.long$age == "Adult.zscore", "Adult", "Child")
child.long2$age <- ifelse(child.long2$age == "Adult.BH", "Adult", "Child")
#selects the columns we want to join together
test <- child.long[,c(1,4:5)]
test2 <- child.long2[,c(1,4:5)]
#joins the columns together
final_long <- left_join(test, test2)

#make numeric values
final_long$zscore<- as.numeric(final_long$zscore)
final_long$BH_P.value <- as.numeric(final_long$BH_P.value)

# Sort the dataframe based on FDR in descending order
order.BH <- final_long %>% arrange(desc(BH_P.value))

order.BH$ordering <- ifelse(order.BH$age == "Child", order.BH$zscore, NA)
order.BH  <- order.BH  %>%
  arrange(ordering, Upstream.Regulators) %>%
  mutate(Upstream.Regulators = factor(Upstream.Regulators, levels = unique(Upstream.Regulators)))

#bargraph, pathways visulalisation
pdf("child.higher.upstream.cytokine.pdf", h=8.25 , w=11.75)
graph2 <-ggplot(order.BH, aes(x = Upstream.Regulators, y= zscore, fill=age)) + 
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values=c("#5690CC", "#EF3D39"))+
  scale_x_discrete(position = "bottom")+
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "z-score", x= "Upstream regulators", title = "Cytokines upstream regulators- BH sig")
graph2
dev.off()

#Transcription Factors

#data is significantly enriched upstream regulators-comparing children and adults
Child.table <- read.csv("ChildUP.TF.csv")
Child.comp <- read.csv("Monos.ChildUp.upstream.csv")

matching_child <- Child.comp[Child.comp$Upstream.Regulators %in% Child.table$Upstream.Regulator, ]

matching_child_noNA <- matching_child %>%
  filter(Child.zscore!= "N/A" & Adult.zscore!= "N/A")

##make dataframe long--- 
#gathers z.scores/BH-P.values into one column and adds a new column for age
child.long <- pivot_longer(matching_child_noNA, cols=c(2:3), names_to = "age", values_to= "zscore")
#gathers FDRs into one column and adds a new column for age
child.long2 <- pivot_longer(matching_child_noNA, cols=c(4:5), names_to = "age", values_to= "BH_P.value")
#changes the names of the group to just Adult and Child
child.long$age <- ifelse(child.long$age == "Adult.zscore", "Adult", "Child")
child.long2$age <- ifelse(child.long2$age == "Adult.BH", "Adult", "Child")
#selects the columns we want to join together
test <- child.long[,c(1,4:5)]
test2 <- child.long2[,c(1,4:5)]
#joins the columns together
final_long <- left_join(test, test2)

#make numeric values
final_long$zscore<- as.numeric(final_long$zscore)
final_long$BH_P.value <- as.numeric(final_long$BH_P.value)

# Sort the dataframe based on FDR in descending order
order.BH <- final_long %>% arrange(desc(BH_P.value))

order.BH$ordering <- ifelse(order.BH$age == "Child", order.BH$zscore, NA)
order.BH  <- order.BH  %>%
  arrange(ordering, Upstream.Regulators) %>%
  mutate(Upstream.Regulators = factor(Upstream.Regulators, levels = unique(Upstream.Regulators)))

#bargraph, pathways visulalisation
pdf("child.higher.upstream.TF.pdf", h=8.25 , w=11.75)
graph2 <-ggplot(order.BH, aes(x = Upstream.Regulators, y= zscore, fill=age)) + 
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values=c("#5690CC", "#EF3D39"))+
  scale_x_discrete(position = "bottom")+
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "z-score", x= "TF", title = "TF upstream regulators- BH sig")
graph2
dev.off()

#Transmembrane receptors

#data is significantly enriched upstream regulators-comparing children and adults
Child.table <- read.csv("ChildUP.TM.csv")
Child.comp <- read.csv("Monos.ChildUp.upstream.csv")

matching_child <- Child.comp[Child.comp$Upstream.Regulators %in% Child.table$Upstream.Regulator, ]

matching_child_noNA <- matching_child %>%
  filter(Child.zscore!= "N/A" & Adult.zscore!= "N/A")

##make dataframe long--- 
#gathers z.scores/BH-P.values into one column and adds a new column for age
child.long <- pivot_longer(matching_child_noNA, cols=c(2:3), names_to = "age", values_to= "zscore")
#gathers FDRs into one column and adds a new column for age
child.long2 <- pivot_longer(matching_child_noNA, cols=c(4:5), names_to = "age", values_to= "BH_P.value")
#changes the names of the group to just Adult and Child
child.long$age <- ifelse(child.long$age == "Adult.zscore", "Adult", "Child")
child.long2$age <- ifelse(child.long2$age == "Adult.BH", "Adult", "Child")
#selects the columns we want to join together
test <- child.long[,c(1,4:5)]
test2 <- child.long2[,c(1,4:5)]
#joins the columns together
final_long <- left_join(test, test2)

#make numeric values
final_long$zscore<- as.numeric(final_long$zscore)
final_long$BH_P.value <- as.numeric(final_long$BH_P.value)

# Sort the dataframe based on FDR in descending order
order.BH <- final_long %>% arrange(desc(BH_P.value))

order.BH$ordering <- ifelse(order.BH$age == "Child", order.BH$zscore, NA)
order.BH  <- order.BH  %>%
  arrange(ordering, Upstream.Regulators) %>%
  mutate(Upstream.Regulators = factor(Upstream.Regulators, levels = unique(Upstream.Regulators)))

#bargraph, pathways visulalisation
pdf("child.higher.upstream.TM.pdf", h=8.25 , w=11.75)
graph2 <-ggplot(order.BH, aes(x = Upstream.Regulators, y= zscore, fill=age)) + 
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values=c("#5690CC", "#EF3D39"))+
  scale_x_discrete(position = "bottom")+
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "z-score", x= "TM", title = "TM upstream regulators- BH sig")
graph2
dev.off()

##################################################################################################################################
##################################################################################################################################
#Baseline and after stim comparison - Venn diagrams for children and adults

#So now I am thinking about when we do the Venn diagram of genes which are up in adults at baseline and up in adults after stim,
#we see CXCL9/CXCL10 fall out as the top genes which are in both groups (only 17 genes for in this category).

#So I want to do the same for children and see what genes are in the shared group-- I am hoping that SLAMF7 is there but I dont think
#it will be-- note children first in this code and then adults
#install.packages("VennDiagram")
library(VennDiagram)
library(ggplot2)
library(tidyr)
library(dplyr)

##children comparison code
setwd("F:/2023/Age_bulk")
Stimhigh.child <- read.csv("child.up.csv")
baseline.child <- read.csv(("baseline.child.csv"))


# Define your sets
set1 <- baseline.child
set2 <- Stimhigh.child

set1 <- set1["SYMBOL"]
set2 <- set2["SYMBOL"]

#maybe here I need to make a dataframe not a list
set1_clean <- na.omit(set1)
set2_clean <- na.omit(set2)

set1 <- set1[!is.na(set1$SYMBOL), ]
set2 <- set2[!is.na(set2$SYMBOL), ]

# Create the Venn diagram
venn.plot <- venn.diagram(x = list(set1, set2),
                          category.names = c("Baseline", "Higher after stim"),
                          filename = NULL, 
                          output = TRUE,
                          area.proportional = TRUE, 
                          category.col = c("red", "blue"),
                          scaled = TRUE,
                          fontface = "bold")

# Display the Venn diagram
grid.draw(venn.plot)

# Find the intersections
intersection_1_2 <- intersect(set1, set2)

# Create a data frame to represent the intersections
intersection_df <- data.frame(
  Sets = c("Set1", "Set2", "Set1 & Set2"),
  Size = c(length(set1), length(set2),
           length(intersection_1_2)))
print(intersection_df)

#made dataframe of all shared genes
shared.genes <- data.frame(gene = intersection_1_2)

#but now what I need is to come back and have all the other information of those genes
shared.genes <- subset(Stimhigh.child, SYMBOL %in% intersection_1_2)

##BAR GRAPHS OF THESE
shared.state <- shared.genes %>%
  select(SYMBOL, qvals.state_bool, qvals.age_bool, qvals.state_bool.age_bool, State_AdultFC, State_ChildFC,  Age_stim, Age_unstim, Sig_cat, Sig_cat2) %>%
  dplyr::rename(Adult = State_AdultFC, Child = State_ChildFC)

#make long
shared.state.long <- gather(shared.state, 
                            Age, 
                            State_FC,
                            Adult:Child, 
                            factor_key = TRUE)


shared.state.long$ordering <- ifelse(shared.state.long$Age == "Child", shared.state.long$State_FC, NA)
shared.state.long  <- shared.state.long  %>%
  arrange(ordering, SYMBOL) %>%
  mutate(SYMBOL = factor(SYMBOL, levels = unique(SYMBOL)))

#bargraph, pathways visulalisation
pdf("child.sharedgenes.fromVENN.pdf", h=8.25 , w=11.75)
graph4 <-ggplot(shared.state.long, aes(x = SYMBOL, y= State_FC, fill=Age)) + 
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values=c("#5690CC", "#EF3D39"))+
  scale_x_discrete(position = "bottom")+
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "State_FC", x= "Gene", title = "Genes are shared between baseline and higher in adults after stim")
graph4
dev.off()

#to look at the genes which are uniquely up in children after stim (105 genes)

# Get the genes unique to each set
unique_set1 <- setdiff(set1_clean$SYMBOL, set2_clean$SYMBOL)
unique_set2 <- setdiff(set2_clean$SYMBOL, set1_clean$SYMBOL)

# Create data frames for unique genes
unique_set1_df <- data.frame(SYMBOL = unique_set1)
unique_set2_df <- data.frame(SYMBOL = unique_set2)

# If you want to subset the original data with unique genes
unique_set1_data <- subset(baseline.child, SYMBOL %in% unique_set1)
unique_set2_data <- subset(Stimhigh.child, SYMBOL %in% unique_set2)

##BAR GRAPHS OF THESE
shared.state <- unique_set2_data %>%
  select(SYMBOL, qvals.state_bool, qvals.age_bool, qvals.state_bool.age_bool, State_AdultFC, State_ChildFC,  Age_stim, Age_unstim, Sig_cat, Sig_cat2) %>%
  dplyr::rename(Adult = State_AdultFC, Child = State_ChildFC)

#make long
shared.state.long <- gather(shared.state, 
                            Age, 
                            State_FC,
                            Adult:Child, 
                            factor_key = TRUE)


shared.state.long$ordering <- ifelse(shared.state.long$Age == "Child", shared.state.long$State_FC, NA)
shared.state.long  <- shared.state.long  %>%
  arrange(ordering, SYMBOL) %>%
  mutate(SYMBOL = factor(SYMBOL, levels = unique(SYMBOL)))

#bargraph, pathways visulalisation
pdf("child.higherafterstim.fromVENN.pdf", h=8.25 , w=11.75)
graph5 <-ggplot(shared.state.long, aes(x = SYMBOL, y= State_FC, fill=Age)) + 
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values=c("#5690CC", "#EF3D39"))+
  scale_x_discrete(position = "bottom")+
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "State_FC", x= "Gene", title = "Genes are shared between baseline and higher in adults after stim")
graph5
dev.off()

###########################################################################
#adult comparison code
setwd("F:/2023/Age_bulk")
Stimhigh.adult <- read.csv("adult.up.csv")
baseline.adult <- read.csv(("baseline.adult.csv"))

# Define your sets
set1 <- baseline.adult
set2 <- Stimhigh.adult

set1 <- set1["SYMBOL"]
set2 <- set2["SYMBOL"]

#maybe here I need to make a dataframe not a list
set1_clean <- na.omit(set1)
set2_clean <- na.omit(set2)

set1 <- set1[!is.na(set1$SYMBOL), ]
set2 <- set2[!is.na(set2$SYMBOL), ]


# Create the Venn diagram
venn.plot <- venn.diagram(x = list(set1, set2),
                          category.names = c("Baseline", "Higher after stim"),
                          filename = NULL, 
                          output = TRUE,
                          area.proportional = TRUE, 
                          category.col = c("red", "blue"),
                          scaled = TRUE,
                          fontface = "bold")

# Display the Venn diagram
grid.draw(venn.plot)

# Find the intersections
intersection_1_2 <- intersect(set1, set2)

# Create a data frame to represent the intersections
intersection_df <- data.frame(
  Sets = c("Set1", "Set2", "Set1 & Set2"),
  Size = c(length(set1), length(set2),
           length(intersection_1_2)))
print(intersection_df)

#made dataframe of all shared genes
shared.genes <- data.frame(gene = intersection_1_2)

#but now what I need is to come back and have all the other information of those genes
shared.genes <- subset(Stimhigh.adult, SYMBOL %in% intersection_1_2)

##BAR GRAPHS OF THESE
shared.state <- shared.genes %>%
  select(SYMBOL, qvals.state_bool, qvals.age_bool, qvals.state_bool.age_bool, State_AdultFC, State_ChildFC,  Age_stim, Age_unstim, Sig_cat, Sig_cat2) %>%
  dplyr::rename(Adult = State_AdultFC, Child = State_ChildFC)

#make long
shared.state.long <- gather(shared.state, 
                            Age, 
                            State_FC,
                            Adult:Child, 
                            factor_key = TRUE)


shared.state.long$ordering <- ifelse(shared.state.long$Age == "Adult", shared.state.long$State_FC, NA)
shared.state.long  <- shared.state.long  %>%
  arrange(ordering, SYMBOL) %>%
  mutate(SYMBOL = factor(SYMBOL, levels = unique(SYMBOL)))

#subset.state<- shared.state.long[1:50, ]
#subset.state <- subset(shared.state.long, shared.state.long$ordering > 1)

#bargraph, pathways visulalisation
pdf("adult.sharedgenes.fromVENN.pdf", h=8.25 , w=11.75)
graph4 <-ggplot(shared.state.long, aes(x = SYMBOL, y= State_FC, fill=Age)) + 
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values=c("#5690CC", "#EF3D39"))+
  scale_x_discrete(position = "bottom")+
  scale_y_continuous(position = "right") +
  coord_flip() +
  theme_bw() + 
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(axis.title.y = element_text(size=12)) + theme (axis.title.x = element_text(size=10))+
  labs(y= "State_FC", x= "Gene", title = "Genes are shared between baseline and higher in adults after stim")
graph4
dev.off()

###############################################################################################################################################################
###############################################################################################################################################################







