# SKCM Cancer BRAF
# DT Laura Vo Ngoc
# Start: 15/06/2018

library("BSgenome.Hsapiens.UCSC.hg38")
library(deconstructSigs)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(TCGA2STAT)
library(TCGAbiolinks)
require(FirebrowseR)

#--------- WD & LOAD FILES --------------------------------------------------------------------------------------------
setwd("C:/Users/rockp/Desktop/UCL/Project/Project-Laura")

load("./Data/TCGA_SKCM_mutations.RData")             # Raw TCGA mutation data
load("./Output/TCGA_SKCM_snvs.RData")      # Selected SNVs

# Primary tumor
load("./Output/TCGA_SKCM_TP_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_SKCM_TP_weights_cut0.00.RData")     # Mutational signatures output

# Metastatic tumor
load("./Output/TCGA_SKCM_TM_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_SKCM_TM_weights_cut0.00.RData")     # Mutational signatures output

#---------  DATA OF INTEREST ------------------------------------------------------------------------------------------

mut.skcm <- getTCGA(disease="SKCM", data.type="Mutation", type="somatic")