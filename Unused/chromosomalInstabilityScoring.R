# Chromosomal instability scoring
# DT Laura Vo Ngoc
# Start: 05/06/2018

library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
require(FirebrowseR)


#--------- WD & LOAD FILES --------------------------------------------------------------------------------------------
setwd("C:/Users/rockp/Desktop/UCL/Project/Project-Laura")

# Mutational signatures and clinical data
load("./Output/TCGA_SKCM_TP_clinical.RData")    # SKCM TP
load('./Output/TCGA_SKCM_TM_clinical.RData')    # SKCM TM
load("./Output/TCGA_UVM_TP_clinical.RData")     # UVM TP 

# Stage
load("./Output/TCGA_SKCM_TP_stage.RData")       # SKCM TP
load("./Output/TCGA_SKCM_TM_stage.RData")       # SKCM TM
load("./Output/TCGA_UVM_TP_stage.RData")        # UVM TP

# Expression matrix
load("./Output/TCGA_SKCM_TP_expression_mel_matrix.RData")  # SKCM TP
load("./Output/TCGA_SKCM_TM_expression_mel_matrix.RData")  # SKCM TM
load("./Output/TCGA_UVM_TP_expression_mel_matrix.RData")   # UVM TP

# Metastatic score
load("./Output/TCGA_SKCM_TP_metastatic_score_mel.RData")  # SKCM TP
load("./Output/TCGA_SKCM_TM_metastatic_score_mel.RData")  # SKCM TM
load("./Output/TCGA_UVM_TP_metastatic_score_mel.RData")   # UVM TP


#--------- COLORS --------------------------------------------------------------------------------------------