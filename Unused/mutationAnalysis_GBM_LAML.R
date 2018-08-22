# Mutational Analysis - GBM Cancer
# DT Laura Vo Ngoc
# Start: 29/06/2018

library("BSgenome.Hsapiens.UCSC.hg38")
library(deconstructSigs)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(survival)
library(TCGAbiolinks)
require(FirebrowseR)

#--------- WD & LOAD FILES --------------------------------------------------------------------------------------------
setwd("C:/Users/rockp/Desktop/UCL/Project/Project-Laura")

# GBM data
load("./Data/TCGA_GBM_mutations.RData")             # Raw TCGA mutation data
load("./Output/TCGA_GBM_snvs.RData")                # Selected SNVs

# GBM mutsigs
load("./Output/TCGA_GBM_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_GBM_laml.weights_cut0.00.RData")     # Mutational signatures output

# Metastatic tumor
load("./Output/TCGA_GBM_TM_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_GBM_TM_laml.weights_cut0.00.RData")     # Mutational signatures output


# LAML data
load("./Data/TCGA_LAML_mutations.RData")             # Raw TCGA mutation data
load("./Output/TCGA_LAML_snvs.RData")                # Selected SNVs

# Primary tumor
load("./Output/TCGA_LAML_TP_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_LAML_TP_laml.weights_cut0.00.RData")     # Mutational signatures output

# Metastatic tumor
load("./Output/TCGA_LAML_TM_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_LAML_TM_laml.weights_cut0.00.RData")     # Mutational signatures output


order <- c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20","S21","S22","S23","S24","S25","S26","S27","S28","S29","S30","unknown")


#---------  DATA OF INTEREST ------------------------------------------------------------------------------------------

#### ~GBM ####

# Download TCGA mutations for GBM cancer:
gbm.mutations <- GDCquery_Maf("GBM", pipelines = "mutect2")
save(gbm.mutations, file="./Data/TCGA_GBM_mutations.RData")

# Select columns of interest
gbm.snvs <- data.frame(gbm.mutations[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                "Chromosome","Start_Position","End_Position",
                                "Variant_Classification","Variant_Type",
                                "Reference_Allele","Tumor_Seq_Allele1",
                                "Tumor_Seq_Allele2")])
gbm.snvs$Sample_Type <- gbm.snvs$Tumor_Sample_Barcode

# Only interested in point mutations, not small insertions/deletions:
gbm.snvs <- gbm.snvs[which((gbm.snvs$Tumor_Seq_Allele2 %in% 
                        c("A","C","G","T"))&
                       (gbm.snvs$Reference_Allele %in%
                            c("A","C","G","T"))),]

# Extract primary tumor samples ### 01 = primary tumor, 06 = metastatic
gbm.snvs$Sample_Type <- sapply(gbm.snvs$Sample_Type, function(x) strsplit(x,"-")[[1]][4])
gbm.snvs$Sample_Type <- sapply(gbm.snvs$Sample_Type, function(x) strsplit(x, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", 
                                                                  perl=TRUE)[[1]][1])

# Save
save(gbm.snvs, file="./Output/TCGA_GBM_snvs.RData")



#### ~LAML ####

# Download TCGA mutations for LAML cancer:
laml.mutations <- GDCquery_Maf("LAML", pipelines = "mutect2")
save(laml.mutations, file="./Data/TCGA_LAML_mutations.RData")

# Select columns of interest
laml.snvs <- data.frame(laml.mutations[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                    "Chromosome","Start_Position","End_Position",
                                    "Variant_Classification","Variant_Type",
                                    "Reference_Allele","Tumor_Seq_Allele1",
                                    "Tumor_Seq_Allele2")])
laml.snvs$Sample_Type <- laml.snvs$Tumor_Sample_Barcode

# Only interested in point mutations, not small insertions/deletions:
laml.snvs <- laml.snvs[which((laml.snvs$Tumor_Seq_Allele2 %in% 
                                c("A","C","G","T"))&
                               (laml.snvs$Reference_Allele %in%
                                    c("A","C","G","T"))),]

# Extract primary tumor samples ### 01 = primary tumor, 06 = metastatic
laml.snvs$Sample_Type <- sapply(laml.snvs$Sample_Type, function(x) strsplit(x,"-")[[1]][4])
laml.snvs$Sample_Type <- sapply(laml.snvs$Sample_Type, function(x) strsplit(x, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", 
                                                                          perl=TRUE)[[1]][1])

# Save
save(laml.snvs, file="./Output/TCGA_LAML_snvs.RData")

# Stores each type into separate variables
laml.snvs.tp <- laml.snvs[which(laml.snvs$Sample_Type == "01"),]           #  samples
laml.snvs.tm <- laml.snvs[which(laml.snvs$Sample_Type == "06"),]          #  samples


#--------- GBM MUTATIONAL SIGNATURES ANALYSIS -------------------------------------------------------------------------

# Convert to deconstructSigs input:
gbm.sigs.input <- mut.to.sigs.input(mut.ref = gbm.snvs, 
                                   sample.id = "Tumor_Sample_Barcode", 
                                   chr = "Chromosome", 
                                   pos = "Start_Position", 
                                   ref = "Reference_Allele", 
                                   alt = "Tumor_Seq_Allele2",
                                   bsg = BSgenome.Hsapiens.UCSC.hg38)

save(gbm.sigs.input, file="./Output/TCGA_GBM_mutsig_input.RData")


# Signatures
gbm.laml.weights <- as.data.frame(t(sapply(rownames(gbm.sigs.input), 
                                     function(x) whichSignatures(tumor.ref = gbm.sigs.input,
                                                                 signatures.ref = signatures.cosmic,
                                                                 sample.id = x, 
                                                                 contexts.needed = TRUE,
                                                                 signature.cutoff = 0.00,              # default = 0.06
                                                                 tri.counts.method = 'default')$laml.weights)), 
                            row.names=rownames(gbm.sigs.input))

gbm.laml.weight <- as.data.frame(sapply(gbm.laml.weights, function(x) unlist(x)), row.names = rownames(laml.weights))

# Change Signature.# to S#
colnames(gbm.laml.weight) <- sapply(colnames(gbm.laml.weight), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))

# Add unknown contribution column
gbm.laml.weight$unknown <- (1-rowSums(gbm.laml.weight))

save(gbm.laml.weight, file="./Output/TCGA_GBM_laml.weights_cut0.00.RData")


# Sort by S1 (descending)
gbm.sorted <- gbm.laml.weight[order(-gbm.laml.weight$S1),]


# Reformatting the laml.weights data
gbm.laml.weight <- as.matrix(gbm.laml.weight)
gbm.melted <- melt(gbm.laml.weight)
colnames(gbm.melted) <- c("PatientId", "Signature", "laml.weight")
gbm.melted$laml.weight <- sapply(gbm.melted$laml.weight, function(x) 100*x)

# Remove non-contributing signatures
gbm.sorted <- as.matrix(gbm.sorted)
gbm.trim <- gbm.sorted[,which(colSums(gbm.sorted)>0)]
gbm.melt <- melt(gbm.trim)
colnames(gbm.melt) <- c("PatientId", "Signature", "laml.weight")
gbm.melt$laml.weight <- sapply(gbm.melt$laml.weight, function(x) 100*x)


# Plot: stacked barplot signatures
colors <- c(brewer.pal(12, "Paired"), brewer.pal(11, "Set3"),"#000000")
colors2 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#BFA89E", "#826C7F", "#FF7F00",
             "#CAB2D6", "#8B786D", "#6A3D9A", "#D6DBD2", "#FFFF99", "#B15928", "#8DD3C7", "#FFFFB3", "#E75A7C", "#BEBADA",
             "#FB8072", "#40476D", "#80B1D3", "#FDB462", "#B3DE69", "#258EA6", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", 
             "#000000")
pdf("./Figures/TCGA_GBM_TP_mutsig.pdf", w=10, h=6)
ggplot(gbm.melt, aes(x=PatientId, y=laml.weight, fill=Signature)) + 
    scale_fill_manual(values=colors2) +
    geom_bar(stat="identity") +
    ggtitle("Mutational signature analysis GBM primary tumor") +
    xlab("Sample") +
    ylab("Contribution (%)") +
    coord_cartesian(ylim=c(0,100), expand=FALSE) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
dev.off()

# Plot: boxplot signatures
pdf("./Figures/TCGA_GBM_TP_mutsigs_boxplot.pdf", w=10, h=6)
ggplot(gbm.melt, aes(x=Signature, y=laml.weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis GBM primary tumor") +
    xlab("Signatures") +
    ylab("Contribution") +
    coord_cartesian(ylim=c(0,100), expand=FALSE)
dev.off()

# Plot: heatmap
gbm.melted$laml.weight <- as.numeric(gbm.melted$laml.weight)
gbm.mat <- acast(PatientId~Signature,data=gbm.melted, value="laml.weight", fun.aggregate=mean)
pdf("./Figures/TCGA_GBM_mutsigs_heatmap.pdf", w=10, h=6, onefile=FALSE)
p <- pheatmap(t(gbm.mat), show_colnames = FALSE, main="Mutational signatures GBM primary tumor")
dev.off()

#--------- ~ Mean contributions per signature ---------

gbm.means <- data.frame(Signatures=c(colnames(gbm.trim)), value=c(colMeans(gbm.trim)*100))

pdf("./Figures/TCGA_GBM_TM_mutsigs_avg.pdf", w=13, h=6)
ggplot(gbm.means, aes(x=reorder(Signatures, -value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    ggtitle("Mean mutational signature contributions GBM primary") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,70), expand=FALSE)
dev.off()

save(gbm.means, file="./Output/TCGA_GBM_TP_mutsigs_means.RData")

#--------- LAML MUTATIONAL SIGNATURES ANALYSIS -------------------------------------------------------------------------

# Convert to deconstructSigs input:
laml.sigs.input <- mut.to.sigs.input(mut.ref = laml.snvs, 
                                   sample.id = "Tumor_Sample_Barcode", 
                                   chr = "Chromosome", 
                                   pos = "Start_Position", 
                                   ref = "Reference_Allele", 
                                   alt = "Tumor_Seq_Allele2",
                                   bsg = BSgenome.Hsapiens.UCSC.hg38)

save(laml.sigs.input, file="./Output/TCGA_LAML_mutsig_input.RData")


# Signatures
laml.laml.weights <- as.data.frame(t(sapply(rownames(laml.sigs.input), 
                                     function(x) whichSignatures(tumor.ref = laml.sigs.input,
                                                                 signatures.ref = signatures.cosmic,
                                                                 sample.id = x, 
                                                                 contexts.needed = TRUE,
                                                                 signature.cutoff = 0.00,              # default = 0.06
                                                                 tri.counts.method = 'default')$laml.weights)), 
                            row.names=rownames(laml.sigs.input))

laml.laml.weight <- as.data.frame(sapply(laml.laml.weights, function(x) unlist(x)), row.names = rownames(laml.laml.weights))

# Change Signature.# to S#
colnames(laml.laml.weight) <- sapply(colnames(laml.laml.weight), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))

# Add unknown contribution column
laml.laml.weight$unknown <- (1-rowSums(laml.laml.weight))

save(laml.laml.weight, file="./Output/TCGA_LAML_laml.weights_cut0.00.RData")



# Sort by S1 (descending)
laml.sorted <- laml.weight[order(-laml.weight$S1),]


# Reformatting the laml.weights data 
laml.weight <- as.matrix(laml.weight)
laml.melted <- melt(laml.weight)
colnames(laml.melted) <- c("PatientId", "Signature", "weight")
laml.melted$weight <- sapply(laml.melted$weight, function(x) 100*x)


# Remove non-contributing signatures
laml.sorted <- as.matrix(laml.sorted)
laml.trim <- laml.sorted[,which(colSums(laml.sorted)>0)]
laml.melt <- melt(laml.trim)
colnames(laml.melt) <- c("PatientId", "Signature", "weight")
laml.melt$weight <- sapply(laml.melt$weight, function(x) 100*x)

# Plot: stacked barplot signatures
pdf("./Figures/TCGA_LAML_TB_mutsig.pdf", w=10, h=6)
ggplot(melt, aes(x=PatientId, y=laml.weight, fill=Signature)) + 
    scale_fill_manual(values=colors2) +
    geom_bar(stat="identity") +
    ggtitle("Mutational signature analysis LAML metastatic tumor ") +
    xlab("Sample") +
    ylab("Contribution (%)") +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    coord_cartesian(ylim=c(0,100), expand=FALSE)
dev.off()

# Plot: boxplot signatures
pdf("./Figures/TCGA_LAML_TB_mutsigs_boxplot.pdf", w=10, h=6)
ggplot(melt, aes(x=Signature, y=laml.weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis LAML metastatic tumor") +
    xlab("Signatures") +
    ylab("Contribution") +
    coord_cartesian(ylim=c(0,100), expand=FALSE)
dev.off()

# Plot: heatmap
melted$laml.weight <- as.numeric(melted$laml.weight)
tm.mat <- acast(PatientId~Signature, data=melted, value="laml.weight", fun.aggregate=mean)
pdf("./Figures/TCGA_LAML_TB_mutsigs_heatmap.pdf", w=10, h=6, onefile=FALSE)
p <- pheatmap(t(tm.mat), show_colnames = FALSE, main="Mutational signatures LAML metastatic tumor")
dev.off()


#--------- ~ Mean contributions per signature ---------

means <- data.frame(Signatures=c(colnames(trim)), value=c(colMeans(trim)*100))

pdf("./Figures/TCGA_LAML_TB_mutsigs_avg.pdf", w=13, h=6)
ggplot(means, aes(x=reorder(Signatures, -value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    ggtitle("Mean mutational signature contributions LAML metastatic tumor") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,70), expand=FALSE)
dev.off()

save(means, file="./Output/TCGA_LAML_TB_mutsigs_means.RData")



#--------- CLINICAL DATA ----------------------------------------------------------------------------------------------

# FirebrowseR

#### ~GBM ####
load("./Data/TCGA_GBM_clinical.RData")
load("./Output/TCGA_GBM_TP_clinical.RData")
load('./Output/TCGA_GBM_TM_clinical.RData')

# Download data:
GBM.clin = Samples.Clinical(format = "csv",
                             cohort = "GBM",
                             page_size=2000)
save(GBM.clin, file="./Data/TCGA_GBM_clinical.RData")

# Create variables with both clinical and signatures data
gbm.tp.clin <- merge(GBM.clin, tp, by="patient_id")
gbm.tm.clin <- merge(GBM.clin, tm, by="patient_id")
save(gbm.tp.clin, file="./Output/TCGA_GBM_TP_clinical.RData")
save(gbm.tm.clin, file="./Output/TCGA_GBM_TM_clinical.RData")


#### ~LAML ####
load("./Data/TCGA_LAML_clinical.RData")
load("./Output/TCGA_LAML_TP_clinical.RData")
load('./Output/TCGA_LAML_TM_clinical.RData')

# Download data:
LAML.clin = Samples.Clinical(format = "csv",
                            cohort = "LAML",
                            page_size=2000)
save(LAML.clin, file="./Data/TCGA_LAML_clinical.RData")

# Create variables with both clinical and signatures data
laml.tp.clin <- merge(LAML.clin, tp, by="patient_id")
laml.tm.clin <- merge(LAML.clin, tm, by="patient_id")
save(laml.tp.clin, file="./Output/TCGA_LAML_TP_clinical.RData")
save(laml.tm.clin, file="./Output/TCGA_LAML_TM_clinical.RData")