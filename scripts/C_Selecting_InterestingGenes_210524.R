##From analysis performed on 15 May 2024
##Addition of a step setting set.seed(1234)

library('tidyverse')
##This will clean up the local enviroment to start the monocyte analysis from scratch.
rm(list=ls())

monos<- read.delim(file = "output_data/glmmseq_Age_Monocytes_output.txt", sep = "\t")
vd2 <- read.delim(file = "output_data/glmmseq_Age_Vd2_output.txt", sep = "\t")

dim(monos)
#15273    44
dim(vd2)
#14782    44

##Filtering based on qvals

##Filter out anything that is not significant
monos_1<- monos[monos$qvals.state_bool.age_bool<0.05|monos$qvals.state_bool<0.05|monos$qvals.age_bool<0.05,] 
dim(monos_1)
#10562    44

vd2_1<- vd2[vd2$qvals.state_bool.age_bool<0.05|vd2$qvals.state_bool<0.05|vd2$qvals.age_bool<0.05,] 
dim(vd2_1)
#8823   44

##Coding for variable SigG
## age = Sig for age but NS for state or state:age
## state = Sig for state but NS for age or state:age
## state:age = Sig for state:age
## ns = NS for age, state or state:age
# state&age = Sig for state and age but not for interaction
monos_1$SigG[monos_1$qvals.state_bool>= 0.05 & monos_1$qvals.state_bool.age_bool>=0.05 & monos_1$qvals.age_bool>=0.05 ] <- "ns"
monos_1$SigG[monos_1$qvals.state_bool>= 0.05 & monos_1$qvals.state_bool.age_bool>=0.05 & monos_1$qvals.age_bool<=0.05 ] <- "age"
monos_1$SigG[monos_1$qvals.state_bool<= 0.05 & monos_1$qvals.state_bool.age_bool>=0.05 & monos_1$qvals.age_bool>=0.05 ] <- "state"
monos_1$SigG[monos_1$qvals.state_bool<= 0.05 & monos_1$qvals.state_bool.age_bool>=0.05 & monos_1$qvals.age_bool<=0.05 ] <- "state&age"
monos_1$SigG[monos_1$qvals.state_bool.age_bool<0.05] <- "state:age"
table(monos_1$SigG, useNA = "ifany")
#age     state state:age      state&age 
# 167      7660      2271       464 

vd2_1$SigG[vd2_1$qvals.state_bool>= 0.05 & vd2_1$qvals.state_bool.age_bool>=0.05 & vd2_1$qvals.age_bool>=0.05 ] <- "ns"
vd2_1$SigG[vd2_1$qvals.state_bool>= 0.05 & vd2_1$qvals.state_bool.age_bool>=0.05 & vd2_1$qvals.age_bool<=0.05 ] <- "age"
vd2_1$SigG[vd2_1$qvals.state_bool<= 0.05 & vd2_1$qvals.state_bool.age_bool>=0.05 & vd2_1$qvals.age_bool>=0.05 ] <- "state"
vd2_1$SigG[vd2_1$qvals.state_bool<= 0.05 & vd2_1$qvals.state_bool.age_bool>=0.05 & vd2_1$qvals.age_bool<=0.05 ] <- "state&age"
vd2_1$SigG[vd2_1$qvals.state_bool.age_bool<=0.05] <- "state:age"
table(vd2_1$SigG, useNA = "ifany")
#age     state state:age      state&age 
#  57      5704      2991        71 


tapply(vd2_1$qvals.state_bool.age_bool, vd2_1$SigG, summary)
tapply(vd2_1$qvals.state_bool, vd2_1$SigG, summary)
tapply(vd2_1$qvals.age_bool, vd2_1$SigG, summary)


##Selecting only interactions genes & relevant variables
table(monos_1$SigG)
#age     state state:age state&age 
#167      7660      2271       464  
monos_1%>%
  dplyr::select("ENSEMBL", "SYMBOL", "GENENAME", "State_AdultFC", "State_ChildFC", "Age_stim", "Age_unstim",
         "qvals.state_bool", "qvals.age_bool", "qvals.state_bool.age_bool","SigG")%>%
  filter(SigG=="state:age")-> monos_filter

dim(monos_filter)
# 2271   11
table(monos_filter$SigG, useNA = "ifany")
#state:age     
#  2271       

table(vd2_1$SigG)
#age     state state:age 
#57      5704      2991 
vd2_1%>%
  dplyr::select("ENSEMBL", "SYMBOL", "GENENAME", "State_AdultFC", "State_ChildFC", "Age_stim", "Age_unstim",
         "qvals.state_bool", "qvals.age_bool", "qvals.state_bool.age_bool","SigG")%>%
  filter(SigG=="state:age") ->vd2_filter
dim(vd2_filter)
#2991   11

table(vd2_filter$SigG, useNA = "ifany")
#state:age      
#  2991     

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
#5  14   1  32  15  76  61  42  56 448 876 619  21   5

##creating groups for vd2
vd2_filter$groupsAge<-NA
vd2_filter$groupsAge[vd2_filter$Age_unstim<0 & vd2_filter$Age_stim<0 & vd2_filter$State_ChildFC<0 & vd2_filter$State_AdultFC<0] <- "A"#A
vd2_filter$groupsAge[vd2_filter$Age_unstim<0 & vd2_filter$Age_stim<0 & vd2_filter$State_ChildFC<0 & vd2_filter$State_AdultFC>0] <- "B"#B
vd2_filter$groupsAge[vd2_filter$Age_unstim<0 & vd2_filter$Age_stim<0 & vd2_filter$State_ChildFC>0 & vd2_filter$State_AdultFC<0] <- "C"#C
vd2_filter$groupsAge[vd2_filter$Age_unstim<0 & vd2_filter$Age_stim<0 & vd2_filter$State_ChildFC>0 & vd2_filter$State_AdultFC>0] <- "D"#D
vd2_filter$groupsAge[vd2_filter$Age_unstim>0 & vd2_filter$Age_stim<0 & vd2_filter$State_ChildFC>0 & vd2_filter$State_AdultFC>0] <- "E"#E
vd2_filter$groupsAge[vd2_filter$Age_unstim>0 & vd2_filter$Age_stim<0 & vd2_filter$State_ChildFC>0 & vd2_filter$State_AdultFC<0] <- "F"#F
vd2_filter$groupsAge[vd2_filter$Age_unstim>0 & vd2_filter$Age_stim<0 & vd2_filter$State_ChildFC<0 & vd2_filter$State_AdultFC<0] <- "G"#G
vd2_filter$groupsAge[vd2_filter$Age_unstim>0 & vd2_filter$Age_stim>0 & vd2_filter$State_ChildFC>0 & vd2_filter$State_AdultFC<0] <- "H"#H
vd2_filter$groupsAge[vd2_filter$Age_unstim>0 & vd2_filter$Age_stim>0 & vd2_filter$State_ChildFC<0 & vd2_filter$State_AdultFC<0] <- "I"#I
vd2_filter$groupsAge[vd2_filter$Age_unstim<0 & vd2_filter$Age_stim>0 & vd2_filter$State_ChildFC<0 & vd2_filter$State_AdultFC<0] <- "J"#J
vd2_filter$groupsAge[vd2_filter$Age_unstim<0 & vd2_filter$Age_stim>0 & vd2_filter$State_ChildFC<0 & vd2_filter$State_AdultFC>0] <- "K"#K
vd2_filter$groupsAge[vd2_filter$Age_unstim<0 & vd2_filter$Age_stim>0 & vd2_filter$State_ChildFC>0 & vd2_filter$State_AdultFC>0] <- "L"#L
vd2_filter$groupsAge[vd2_filter$Age_unstim>0 & vd2_filter$Age_stim>0 & vd2_filter$State_ChildFC>0 & vd2_filter$State_AdultFC>0] <- "M"#M
vd2_filter$groupsAge[vd2_filter$Age_unstim>0 & vd2_filter$Age_stim>0 & vd2_filter$State_ChildFC<0 & vd2_filter$State_AdultFC>0] <- "N"#N
table(vd2_filter$groupsAge)
#A    C    D    E    F    G    H    I    J    K    L    M    N 
#6    5   18 1295  689   61  205   44    6    3   21  636    2

##vd2 data doesn't have group "B"genes, for plotting purposes we can add this information before drawing the barplot.


##Creating direction groups
monos_filter$Direction_Ch_Ad<-NA
monos_filter$Direction_Ch_Ad[monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC>0] <- "UP_UP"
monos_filter$Direction_Ch_Ad[monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC<0] <- "DOWN_DOWN"
monos_filter$Direction_Ch_Ad[monos_filter$State_ChildFC<0 & monos_filter$State_AdultFC>0] <- "DOWN_UP"
monos_filter$Direction_Ch_Ad[monos_filter$State_ChildFC>0 & monos_filter$State_AdultFC<0] <- "UP_DOWN"
table(monos_filter$Direction_Ch_Ad)

#DOWN_DOWN   DOWN_UP   UP_DOWN     UP_UP 
#570       895       119       687 

vd2_filter$Direction_Ch_Ad<-NA
vd2_filter$Direction_Ch_Ad[vd2_filter$State_ChildFC>0 & vd2_filter$State_AdultFC>0] <- "UP_UP"
vd2_filter$Direction_Ch_Ad[vd2_filter$State_ChildFC<0 & vd2_filter$State_AdultFC<0] <- "DOWN_DOWN"
vd2_filter$Direction_Ch_Ad[vd2_filter$State_ChildFC<0 & vd2_filter$State_AdultFC>0] <- "DOWN_UP"
vd2_filter$Direction_Ch_Ad[vd2_filter$State_ChildFC>0 & vd2_filter$State_AdultFC<0] <- "UP_DOWN"
table(vd2_filter$Direction_Ch_Ad)

#DOWN_DOWN   DOWN_UP   UP_DOWN     UP_UP 
#117         5       899      1970 

##Creating variable higher expression pre-stim
monos_filter$HigherExp_prestim<-NA
monos_filter$HigherExp_prestim[monos_filter$Age_unstim<0] <- "Children"
monos_filter$HigherExp_prestim[monos_filter$Age_unstim>0] <- "Adults"
table(monos_filter$HigherExp_prestim)

#Adults Children 
#276     1995 

vd2_filter$HigherExp_prestim<-NA
vd2_filter$HigherExp_prestim[vd2_filter$Age_unstim<0] <- "Children"
vd2_filter$HigherExp_prestim[vd2_filter$Age_unstim>0] <- "Adults"
table(vd2_filter$HigherExp_prestim)

#Adults Children 
#2932       59 

##Creating variable higher expression post-stim
monos_filter$HigherExp_poststim<-NA
monos_filter$HigherExp_poststim[monos_filter$Age_stim<0] <- "Children"
monos_filter$HigherExp_poststim[monos_filter$Age_stim>0] <- "Adults"

table(monos_filter$HigherExp_poststim)
#Adults Children 
#2067      204 

vd2_filter$HigherExp_poststim<-NA
vd2_filter$HigherExp_poststim[vd2_filter$Age_stim<0] <- "Children"
vd2_filter$HigherExp_poststim[vd2_filter$Age_stim>0] <- "Adults"

table(vd2_filter$HigherExp_poststim)

#Adults Children 
#917     2074  

write.table(monos_filter, file = "output_data/Age_monos_SiGenelist_withInteraction_210524.txt", sep = "\t",
            row.names = TRUE)
write.table(vd2_filter, file = "output_data/Age_vd2_SiGenelist_withInteraction_210524.txt", sep = "\t",
            row.names = TRUE)

##Suggested barplot code for VD2 data missing B group

VD2_bardata <- as.data.frame(table(vd2_filter$groupsAge,vd2_filter$HigherExp_prestim, vd2_filter$HigherExp_poststim,vd2_filter$Direction_Ch_Ad))
VD2_bardata1 <- VD2_bardata[!VD2_bardata$Freq==0,]
head(VD2_bardata1)

colnames(VD2_bardata1)[1]="Group"
colnames(VD2_bardata1)[2]="HigherExp_prestim"
colnames(VD2_bardata1)[3]="HigherExp_poststim"
colnames(VD2_bardata1)[4]="Direction_Ch_Ad"
colnames(VD2_bardata1)[5]="No of DEGs"

b_row <- as.data.frame(c("x"),c(" "))
b_row$Group <- "B"
b_row$HigherExp_prestim <- NA
b_row$HigherExp_poststim<-NA
b_row$Direction_Ch_Ad <- NA
b_row$'No of DEGs' <- 0
b_row <-b_row[-c(2:5), -1]

head(b_row)
VD2_bardata1 <- rbind(VD2_bardata1, b_row)
head(VD2_bardata1)

VD2_bardata1$Group <- factor(VD2_bardata1$Group, levels=c("A", "B", "C", "D","E", "F", "G", "H", "I", "J", "K", "L", "M", "N"))

VD2_bardata2 <- VD2_bardata1 %>% arrange(Group)
head(VD2_bardata2)

##Adding coloring info
VD2_bardata2%>%
  mutate(
    GeneSelect= case_when(HigherExp_poststim=="Children"& Direction_Ch_Ad=="UP_UP"|HigherExp_poststim=="Children"&Direction_Ch_Ad=="UP_DOWN"~"Child",
                          HigherExp_poststim=="Adults"& Direction_Ch_Ad=="UP_UP"|HigherExp_poststim=="Adults"& Direction_Ch_Ad=="DOWN_UP"~"Adult",
                          TRUE~"Not selected"))-> VD2_bardata3


table(VD2_bardata3$GeneSelect)
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
summary.bar <- ggplot(VD2_bardata3, aes(x = Group, y =`No of DEGs`, fill = GeneSelect)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values=boxplot.colour)+
  geom_text(aes(label=`No of DEGs`), vjust=-0.5)+
  theme_bw() +
  labs(title = " ", x = "Group", y = "Interaction genes")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
summary.bar

ggsave(plot=summary.bar, filename = "graphs/VD2_barplot_geneSelect_210524.pdf", width = 7, height = 5, dpi = 600,scale = 1)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

