# Mutational Analysis - UVM Cancer
# DT Laura Vo Ngoc
# Start: 06/03/2018

library("BSgenome.Hsapiens.UCSC.hg38")
library(deconstructSigs)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(TCGAbiolinks)

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

#### WD & LOAD FILES ####
setwd("C:/Users/rockp/Desktop/UCL/Project/Project-Laura")

load("./Data/TCGA_UVM_mutations.RData")             # Raw TCGA mutation data
load("./Output/TCGA_UVM_snvs.RData")      # Selected SNVs

# Primary tumor
load("./Output/TCGA_UVM_TP_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_UVM_TP_weights_cut0.00.RData")     # Mutational signatures output


#**********************************************************************************************************************
#### DATA OF INTEREST ####

# Download TCGA mutations for SKCM cancer:
mutations <- GDCquery_Maf("UVM", pipelines = "mutect2")
save(mutations, file="./Data/TCGA_UVM_mutations.RData")

# Select columns of interest
snvs <- data.frame(mutations[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                "Chromosome","Start_Position","End_Position",
                                "Variant_Classification","Variant_Type",
                                "Reference_Allele","Tumor_Seq_Allele1",
                                "Tumor_Seq_Allele2")])
snvs$Sample_Type <- snvs$Tumor_Sample_Barcode

# only interested in point mutations, not small insertions/deletions:
snvs <- snvs[which((snvs$Tumor_Seq_Allele2 %in% 
                        c("A","C","G","T"))&
                       (snvs$Reference_Allele %in%
                            c("A","C","G","T"))),]

# Extract primary tumor samples ### 01 = primary tumor, 06 = metastatic
snvs$Sample_Type <- sapply(snvs$Sample_Type, function(x) strsplit(x,"-")[[1]][4])
snvs$Sample_Type <- sapply(snvs$Sample_Type, function(x) strsplit(x, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)[[1]][1])

# Save
save(snvs, file="./Output/TCGA_UVM_snvs.RData")

# Stores each type into separate variables
snvs_tp <- snvs[which(snvs$Sample_Type == "01"),]
snvs_met <- snvs[which(snvs$Sample_Type == "06"),]


#**********************************************************************************************************************
#### MUTATIONAL SIGNATURES ANALYSIS -- TP ####

# Convert to deconstructSigs input:
sigs_input_tp <- mut.to.sigs.input(mut.ref = snvs_tp, 
                                   sample.id = "Tumor_Sample_Barcode", 
                                   chr = "Chromosome", 
                                   pos = "Start_Position", 
                                   ref = "Reference_Allele", 
                                   alt = "Tumor_Seq_Allele2",
                                   bsg = BSgenome.Hsapiens.UCSC.hg38)

save(sigs_input_tp, file="./Output/TCGA_UVM_TP_mutsig_input.RData")


# Signatures
weights_tp <- as.data.frame(t(sapply(rownames(sigs_input_tp), 
                                     function(x) whichSignatures(tumor.ref = sigs_input_tp,
                                                                 signatures.ref = signatures.cosmic,
                                                                 sample.id = x, 
                                                                 contexts.needed = TRUE,
                                                                 signature.cutoff = 0.06,              # default = 0.06
                                                                 tri.counts.method = 'default')$weights)), 
                            row.names=rownames(sigs_input_tp))

weight_tp <- as.data.frame(sapply(weights_tp, function(x) unlist(x)), row.names = rownames(weights_tp))

# Change Signature.# to S#
colnames(weight_tp) <- sapply(colnames(weight_tp), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))

# Add unknown contribution column
weight_tp$unknown <- (1-rowSums(weight_tp))

save(weight_tp, file="./Output/TCGA_UVM_TP_weights_cut0.00.RData")



# Sort by S1 (descending)
sorted_tp <- weight_tp[order(-weight_tp$S1),]


# Reformatting the weights data 
weight_tp <- as.matrix(weight_tp)
melted_tp <- melt(weight_tp)
colnames(melted_tp) <- c("PatientId", "Signature", "weight")
melted_tp$weight <- sapply(melted_tp$weight, function(x) 100*x)


# Remove non-contributing signatures
sorted_tp <- as.matrix(sorted_tp)
trim_tp <- sorted_tp[,which(colSums(sorted_tp)>0)]
melt_tp <- melt(trim_tp)
colnames(melt_tp) <- c("PatientId", "Signature", "weight")
melt_tp$weight <- sapply(melt_tp$weight, function(x) 100*x)


# Plot: stacked barplot signatures
colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#BFA89E", "#826C7F", "#FF7F00",
             "#CAB2D6", "#8B786D", "#6A3D9A", "#D6DBD2", "#FFFF99", "#B15928", "#8DD3C7", "#FFFFB3", "#E75A7C", "#BEBADA",
             "#FB8072", "#40476D", "#80B1D3", "#FDB462", "#B3DE69", "#258EA6", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", 
             "#000000")
pdf("./Figures/TCGA_UVM_TP_mutsig.pdf", w=10, h=6)
ggplot(melt_tp, aes(x=PatientId, y=weight, fill=Signature)) + 
    scale_fill_manual(values=colors) +
    geom_bar(stat="identity") +
    ggtitle("Mutational signature analysis UVM primary tumor ") +
    xlab("Sample") +
    ylab("Contribution (%)") +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    coord_cartesian(ylim=c(0,100), expand=FALSE)
dev.off()

# Plot: boxplot signatures
pdf("./Figures/TCGA_UVM_TP_mutsigs_boxplot.pdf", w=10, h=6)
ggplot(melt_tp, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis UVM primary tumor") +
    xlab("Signatures") +
    ylab("Contribution") +
    coord_cartesian(ylim=c(0,100), expand=FALSE)
dev.off()

# Plot: heatmap
melted_tp$weight <- as.numeric(melted_tp$weight)
tp_mat <- acast(PatientId~Signature,data=melted_tp, value="weight", fun.aggregate=mean)
pdf("./Figures/TCGA_UVM_TP_mutsigs_heatmap.pdf", w=10, h=6, onefile=FALSE)
p_tp <- pheatmap(t(tp_mat), show_colnames = FALSE, main="Mutational signatures UVM primary tumor")
dev.off()


#### ~ Mean contributions per signature ####
means <- data.frame(Signatures=c(colnames(trim_tp)), value=c(colMeans(trim_tp)*100))

pdf("./Figures/TCGA_UVM_TP_mutsigs_avg.pdf", w=12, h=6)
ggplot(means, aes(x=reorder(Signatures, -value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    ggtitle("Mean mutational signature contributions UVM primary tumor") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,70), expand=FALSE)
dev.off()

save(means, file="./Output/TCGA_UVM_TP_mutsigs_means.RData")

#**********************************************************************************************************************
#### CLINICAL DATA ####

load("./Data/TCGA_UVM_clinical.RData")
load("./Output/TCGA_UVM_TP_clinical.RData")
load("./Output/TCGA_UVM_TP_stage.RData")

# Mutsigs data frame with patient id column
tp_uvm <- as.data.frame(weight_tp[,which(colSums(weight_tp)>0)])
tp_uvm$patient_id <- sapply(rownames(tp_uvm), function(x) as.factor(tolower(strsplit(x,"-")[[1]][3])))


# Download data:
UVM.clin = Samples.Clinical(format = "csv",
                           cohort = "UVM",
                           page_size=2000)
save(UVM.clin, file="./Data/TCGA_UVM_clinical.RData")

# Create variables with both clinical and signatures data
clin <- merge(UVM.clin, tp_uvm, by="patient_id")

save(clin, file="./Output/TCGA_UVM_TP_clinical.RData")


#### ~ Stage ####

# Get just id, cancer stage, and sigs
stage <- clin[, c(1:2, 15, 97:127)]
stage$stage <- sapply(stage$clinical_stage, function(x) toupper(strsplit(x,"\\ ")[[1]][2]))

# Loop to remove a, b, c specifications of cancer stage
for (i in 1:nrow(stage)) {
    if (grepl("a", stage$stage[i], ignore.case = TRUE) | grepl("b", stage$stage[i], ignore.case = TRUE) | 
        grepl("c", stage$stage[i], ignore.case = TRUE)) {
        stage$stage[i] <- substr(stage$stage[i], 1, nchar(stage$stage[i])-1)
    }
    else {
        stage$stage[i] <- stage$stage[i]
    }
}

# Add column for early/late stage
stage$earlyLate <- stage$stage
for (i in 1:nrow(stage)) {
    if (is.na(stage$earlyLate[i])){
        stage$earlyLate[i] <- NA
    }
    else if (stage$earlyLate[i] == "I" | stage$earlyLate[i] == "II") {
        stage$earlyLate[i] <- "early"
    }
    else if (stage$earlyLate[i] == "III" | stage$earlyLate[i] == "IV") {
        stage$earlyLate[i] <- "late"
    }
}

save(stage, file="./Output/TCGA_UVM_TP_stage.RData")


## Overview
stage.melt <- melt(stage[,c(2, 4:35)])
colnames(stage.melt) <- c("PatientId", "Stage", "Signatures", "weight")
stage.melt$weight <- sapply(stage.melt$weight, function(x) 100*x)

early.late.melt <- melt(stage[,c(2,4:34,36)])
colnames(early.late.melt) <- c("PatientId", "earlyLate", "Signatures", "weight")
early.late.melt$weight <- sapply(early.late.melt$weight, function(x) 100*x)


# Barplot by stage
colors7 <- c("#1F78B4", "#B2DF8A", "#33A02C", "#000000")
pdf("./Figures/TCGA_UVM_TP_mutsigs_stage_barplot.pdf", w=13, h=6)
ggplot(stage.melt, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors7) +
    geom_bar(stat="identity", aes(fill=Stage), position="dodge") +
    ggtitle("Mutational signature contributions UVM TP by cancer stage") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian( ylim=c(0,100), expand = FALSE )
dev.off()


# Mean sig barplot per stage
# Extract mean sig per cancer stage
stagex <- stage[,c(4:35)]
rownames(stagex) <- stage[,2]
stage.mean <- sapply(stagex[, 1:31], function(x) tapply(x, stagex[, 32], mean))
stage.mean <- melt(stage.mean)
colnames(stage.mean) <- c("Stage", "Signatures", "weight")
stage.mean$weight <- sapply(stage.mean$weight, function(x) 100*x)

# Plot barplot
pdf("./Figures/TCGA_UVM_TP_mutsigs_stage_avgbarplot.pdf", w=13, h=6)
ggplot(stage.mean, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors7) +
    geom_bar(stat="identity", aes(fill=Stage), position="dodge") +
    ggtitle("Mean mutational signature contributions UVM TP by cancer stage") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,50), expand = FALSE)
dev.off()

# Mean sig barplot per progression
stagep <- stage[,c(4:34, 36)]
rownames(stagep) <- stage[,2]
stagep.mean <- sapply(stagep[, 1:31], function(x) tapply(x, stagep[, 32], mean))
stagep.mean <- melt(stagep.mean)
colnames(stagep.mean) <- c("Progression", "Signatures", "weight")
stagep.mean$weight <- sapply(stagep.mean$weight, function(x) 100*x) 

# Plot barplot
pdf("./Figures/TCGA_UVM_TP_mutsigs_prog_avgbarplot.pdf", w=13, h=6)
ggplot(stagep.mean, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors4) +
    geom_bar(stat="identity", aes(fill=Progression), position="dodge") +
    ggtitle("Mean mutational signature contributions UVM TP by disease progression") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,50), expand = FALSE)
dev.off()


# Loop over each signature, by stage
# Data frame for t-test output
prog.ttest <- data.frame(Signature=colnames(stage)[4:34], t=NA, df=NA, pval = NA)
for (i in 4:34) {
    sig <- noquote(colnames(stage[i]))
    var <- na.omit(data.frame(Progression=stage$earlyLate, weight=stage[,i]))
    title <- paste0("UVM TP ", sig, " contribution by disease progression")
    file <- paste0("./Figures/TCGA_UVM_TP_prog_", sig, ".pdf")
    
    early <- var[which(var$Progression == "early"),]
    late <- var[which(var$Progression == "late"),]
    ttest <- t.test(early$weight, late$weight)
    
    prog.ttest[which(prog.ttest$Signature==sig),2] <- ttest$statistic
    prog.ttest[which(prog.ttest$Signature==sig),3] <- ttest$parameter
    prog.ttest[which(prog.ttest$Signature==sig),4] <- ttest$p.value
    #prog.ttest[which(prog.ttest$Signature==sig),5] <- ttest$conf.int
    #prog.ttest[which(prog.ttest$Signature==sig),6] <- ttest$estimate
    
    pdf(file, w=7, h=5)
    plot <- ggplot(var, aes(x=Progression, y=weight*100)) +
        #geom_bar(stat="identity", fill="blue") +
        geom_boxplot() +
        ggtitle(title) +
        xlab("Disease progression") +
        ylab("Contribution (%)")
    #coord_cartesian(ylim=c(0,100), expand = FALSE)
    print(plot) 
    dev.off()
}
prog.ttest$padj <- p.adjust(prog.ttest$pval, method="BH")

write.csv(prog.ttest, file="./Output/TCGA_UVM_TP_prog_eachsigttest.csv")

#--------- ~~T/N/M Stage TP -------
tnm <- clin[, c(1:2, 13, 14, 16, 97:127)]
tnm$t <- tnm$clinical_t
tnm$n <- tnm$clinical_n
tnm$m <- tnm$clinical_m


# Loop to remove a, b, c specifications of tnm stage
for (i in 1:nrow(tnm)) {
    if (grepl("a", tnm$t[i], ignore.case = TRUE) | grepl("b", tnm$t[i], ignore.case = TRUE) | 
        grepl("c", tnm$t[i], ignore.case = TRUE) | grepl("d", tnm$t[i], ignore.case = TRUE) | 
        grepl("e", tnm$t[i], ignore.case = TRUE)) {
        tnm$t[i] <- substr(tnm$t[i], 1, nchar(tnm$t[i])-1)
    }
    else {
        tnm$t[i] <- tnm$t[i]
    }
}
for (i in 1:nrow(tnm)) {
    if (grepl("a", tnm$n[i], ignore.case = TRUE) | grepl("b", tnm$n[i], ignore.case = TRUE) | 
        grepl("c", tnm$n[i], ignore.case = TRUE)) {
        tnm$n[i] <- substr(tnm$n[i], 1, nchar(tnm$n[i])-1)
    }
    else {
        tnm$n[i] <- tnm$n[i]
    }
}
for (i in 1:nrow(tnm)) {
    if (grepl("a", tnm$m[i], ignore.case = TRUE) | grepl("b", tnm$m[i], ignore.case = TRUE) | 
        grepl("c", tnm$m[i], ignore.case = TRUE)) {
        tnm$m[i] <- substr(tnm$m[i], 1, nchar(tnm$m[i])-1)
    }
    else {
        tnm$m[i] <- tnm$m[i]
    }
}

#### ~ Therapies ####
# Get list of annotated therapies
therapies.uvm <- array()
for (i in 1:ncol(UVM.clin)) {
    if (grepl("therapy", colnames(UVM.clin[i])) | grepl("treatment", colnames(UVM.clin[i])) | grepl("drug", colnames(UVM.clin)[i])) {
        therapies.uvm <- append(therapies.uvm, colnames(UVM.clin[i]))
    }
}
therapies.uvm

# Summary
table(UVM.clin$therapy_type, useNA="always")
table(UVM.clin$radiation_therapy, useNA="always")
table(UVM.clin$drug_name, useNA="always")