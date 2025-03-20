if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(BiocManager)

BiocManager::install("edgeR")
BiocManager::install("BiocParallel")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("myles-lewis/glmmSeq")
install.packages("statmod")
install.packages("kableExtra")
install.packages("tidyverse")


##packages from data_exploration_Rawc.R
library('edgeR')
library('BiocParallel')
library('tidyverse')
library('statmod')
###If getting this error "Error in initializePtr():
####\n  function 'chm_factor_ldetL2' not provided by package 'Matrix'\n" 
install.packages("lme4", type = "source") 

##packages from clean_glmm_Age_Rawc
library(glmmSeq)
library(edgeR)
library(tidyverse)
library(kableExtra)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)

##For PCA analysis

##reinstalling old matrix version because 1.6.4 doesnt have required function
#remotes::install_version("Matrix", version = "1.5-1")
install.packages("TMB", type = "source")
install.packages("glmmTMB", type = "source")
