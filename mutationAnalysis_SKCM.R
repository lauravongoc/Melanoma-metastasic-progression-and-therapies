# Mutational Analysis - SKCM Cancer
# DT Laura Vo Ngoc
# Start: 18/02/2018

library("BSgenome.Hsapiens.UCSC.hg38")
library(deconstructSigs)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(TCGAbiolinks)
require(FirebrowseR)

#--------- WD & LOAD FILES --------------------------------------------------------------------------------------------
setwd("C:/Users/rockp/Desktop/UCL/Project/Project-Laura")

load("./Data/TCGA_SKCM_mutations.RData")             # Raw TCGA mutation data
load("./Output/TCGA_SKCM_snvs.RData")      # Selected SNVs

# Primary tumor
load("./Output/TCGA_SKCM_PT_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_SKCM_PT_weights_cut0.00.RData")     # Mutational signatures output

# Metastatic tumor
load("./Output/TCGA_SKCM_TM_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_SKCM_TM_weights_cut0.00.RData")     # Mutational signatures output




#---------  DATA OF INTEREST ------------------------------------------------------------------------------------------

# Download TCGA mutations for SKCM cancer:
mutations <- GDCquery_Maf("SKCM", pipelines = "mutect2")
save(mutations, file="./Data/TCGA_SKCM_mutations.RData")

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
snvs$Sample_Type <- sapply(snvs$Sample_Type, function(x) strsplit(x, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", 
                                                                  perl=TRUE)[[1]][1])

# Save
save(snvs, file="./Output/TCGA_SKCM_snvs.RData")

# Stores each type into separate variables
snvs_pt <- snvs[which(snvs$Sample_Type == "01"),]           # 61990 samples
snvs_met <- snvs[which(snvs$Sample_Type == "06"),]          # 326634 samples




#--------- MUTATIONAL SIGNATURES ANALYSIS: PT -------------------------------------------------------------------------

# Convert to deconstructSigs input:
sigs_input_pt <- mut.to.sigs.input(mut.ref = snvs_pt, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

save(sigs_input_pt, file="./Output/TCGA_SKCM_PT_mutsig_input.RData")


# Signatures
weights_pt <- as.data.frame(t(sapply(rownames(sigs_input_pt), 
                                    function(x) whichSignatures(tumor.ref = sigs_input_pt,
                                                                signatures.ref = signatures.cosmic,
                                                                sample.id = x, 
                                                                contexts.needed = TRUE,
                                                                signature.cutoff = 0.00,              # default = 0.06
                                                                tri.counts.method = 'exome')$weights)), 
                           row.names=rownames(sigs_input_pt))

weight_pt <- as.data.frame(sapply(weights_pt, function(x) unlist(x)), row.names = rownames(weights_pt))

# Change Signature.# to S#
colnames(weight_pt) <- sapply(colnames(weight_pt), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))

# Add unknown contribution column
weight_pt$unknown <- (1-rowSums(weight_pt))

save(weight_pt, file="./Output/TCGA_SKCM_PT_weights_cut0.00.RData")



# Sort by S1 (descending)
sorted_pt <- weight_pt[order(-weight_pt$S1),]


# Reformatting the weights data 
weight_pt <- as.matrix(weight_pt)
melted_pt <- melt(weight_pt)
colnames(melted_pt) <- c("PatientId", "Signature", "weight")
melted_pt$weight <- sapply(melted_pt$weight, function(x) 100*x)


# Remove non-contributing signatures
sorted_pt <- as.matrix(sorted_pt)
trim_pt <- sorted_pt[,which(colSums(sorted_pt)>0)]
melt_pt <- melt(trim_pt)
colnames(melt_pt) <- c("PatientId", "Signature", "weight")
melt_pt$weight <- sapply(melt_pt$weight, function(x) 100*x)


# Plot: stacked barplot signatures
colors <- c(brewer.pal(12, "Paired"), brewer.pal(11, "Set3"),"#000000")
pdf("./Figures/TCGA_SKCM_PT_mutsig.pdf", w=10, h=6)
ggplot(melt_pt, aes(x=PatientId, y=weight, fill=Signature)) + 
    scale_fill_manual(values=colors) +
    geom_bar(stat="identity") +
    ggtitle("Mutational signature analysis SKCM primary tumor ") +
    xlab("Sample") +
    ylab("Contribution (%)") +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
dev.off()
    
# Plot: boxplot signatures
pdf("./Figures/TCGA_SKCM_PT_mutsigs_boxplot.pdf", w=10, h=6)
ggplot(melt_pt, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM primary tumor") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

# Plot: heatmap
melted_pt$weight <- as.numeric(melted_pt$weight)
pt_mat <- acast(PatientId~Signature,data=melted_pt, value="weight", fun.aggregate=mean)
pdf("./Figures/TCGA_SKCM_PT_mutsigs_heatmap.pdf", w=10, h=6, onefile=FALSE)
p_pt <- pheatmap(t(pt_mat), show_colnames = FALSE, main="Mutational signatures SKCM primary tumor")
dev.off()


#--------- ~ S1 and S7 analysis: PT ---------
sp1_7 <- as.matrix(weight_pt[,c(1, 7)])
melt_sp1_7 <- melt(sp1_7)
colnames(melt_sp1_7) <- c("PatientId", "Signature", "weight")
melt_sp1_7$weight <- sapply(melt_sp1_7$weight, function(x) 100*x)

# Plot: boxplot signatures
pdf("./Figures/TCGA_SKCM_PT_mutsigs1_7_boxplot.pdf", w=4, h=5)
ggplot(melt_sp1_7, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("SKCM primary tumor S1 vs S7") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()


#--------- ~ Mean contributions per signature ---------

means_pt <- data.frame(Signatures=c(colnames(trim_pt)), value=c(colMeans(trim_pt)*100))

pdf("./Figures/TCGA_SKCM_PT_mutsigs_avg.pdf", w=10, h=6)
ggplot(means_pt, aes(x=reorder(Signatures, -value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    ggtitle("Mean mutational signature contributions SKCM primary tumor") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5))
dev.off()

save(means_pt, file="./Output/TCGA_SKCM_PT_mutsigs_means.RData")



#--------- MUTATIONAL SIGNATURES ANALYSIS: TM -------------------------------------------------------------------------

# Convert to deconstructSigs input:
sigs_input_tm <- mut.to.sigs.input(mut.ref = snvs_met, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

save(sigs_input_tm, file="./Output/TCGA_SKCM_TM_mutsig_input.RData")


# Signatures
weights_tm <- as.data.frame(t(sapply(rownames(sigs_input_tm), 
                                     function(x) whichSignatures(tumor.ref = sigs_input_tm,
                                                                 signatures.ref = signatures.cosmic,
                                                                 sample.id = x, 
                                                                 contexts.needed = TRUE,
                                                                 signature.cutoff = 0.00,              # default = 0.06
                                                                 tri.counts.method = 'exome')$weights)), 
                            row.names=rownames(sigs_input_tm))

weight_tm <- as.data.frame(sapply(weights_tm, function(x) unlist(x)), row.names = rownames(weights_tm))

# Change Signature.# to S#
colnames(weight_tm) <- sapply(colnames(weight_tm), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))

# Add unknown contribution column
weight_tm$unknown <- (1-rowSums(weight_tm))

save(weight_tm, file="./Output/TCGA_SKCM_TM_weights_cut0.00.RData")



# Sort by S1 (descending)
sorted_tm <- weight_tm[order(-weight_tm$S1),]


# Reformatting the weights data 
weight_tm <- as.matrix(weight_tm)
melted_tm <- melt(weight_tm)
colnames(melted_tm) <- c("PatientId", "Signature", "weight")
melted_tm$weight <- sapply(melted_tm$weight, function(x) 100*x)


# Remove non-contributing signatures
sorted_tm <- as.matrix(sorted_tm)
trim_tm <- sorted_tm[,which(colSums(sorted_tm)>0)]
melt_tm <- melt(trim_tm)
colnames(melt_tm) <- c("PatientId", "Signature", "weight")
melt_tm$weight <- sapply(melt_tm$weight, function(x) 100*x)


colors2 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#BFA89E", "#826C7F", "#FF7F00",
             "#CAB2D6", "#8B786D", "#6A3D9A", "#D6DBD2", "#FFFF99", "#B15928", "#8DD3C7", "#FFFFB3", "#E75A7C", "#BEBADA",
             "#FB8072", "#40476D", "#80B1D3", "#FDB462", "#B3DE69", "#258EA6", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", 
             "#000000")
# Plot: stacked barplot signatures
pdf("./Figures/TCGA_SKCM_TM_mutsig.pdf", w=10, h=6)
ggplot(melt_tm, aes(x=PatientId, y=weight, fill=Signature)) + 
    scale_fill_manual(values=colors2) +
    geom_bar(stat="identity") +
    ggtitle("Mutational signature analysis SKCM metastatic tumor ") +
    xlab("Sample") +
    ylab("Contribution (%)") +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
dev.off()

# Plot: boxplot signatures
pdf("./Figures/TCGA_SKCM_TM_mutsigs_boxplot.pdf", w=10, h=6)
ggplot(melt_tm, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM metastatic tumor") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

# Plot: heatmap
melted_tm$weight <- as.numeric(melted_tm$weight)
tm_mat <- acast(PatientId~Signature, data=melted_tm, value="weight", fun.aggregate=mean)
pdf("./Figures/TCGA_SKCM_TM_mutsigs_heatmap.pdf", w=10, h=6, onefile=FALSE)
p_tm <- pheatmap(t(tm_mat), show_colnames = FALSE, main="Mutational signatures SKCM metastatic tumor")
dev.off()


#--------- ~ S1 and S7 analysis: TM ---------
sm1_7 <- as.matrix(weight_tm[,c(1, 7)])
melt_sm1_7 <- melt(sm1_7)
colnames(melt_sm1_7) <- c("PatientId", "Signature", "weight")
melt_sm1_7$weight <- sapply(melt_sm1_7$weight, function(x) 100*x)

# Plot: boxplot signatures
pdf("./Figures/TCGA_SKCM_TM_mutsigs1_7_boxplot.pdf", w=4, h=5)
ggplot(melt_sm1_7, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("SKCM metastatic tumor S1 vs S7") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()


#--------- ~ Mean contributions per signature ---------

means_tm <- data.frame(Signatures=c(colnames(trim_tm)), value=c(colMeans(trim_tm)*100))

pdf("./Figures/TCGA_SKCM_TM_mutsigs_avg.pdf", w=13, h=6)
ggplot(means_tm, aes(x=reorder(Signatures, -value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    ggtitle("Mean mutational signature contributions SKCM metastatic tumor") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5))
dev.off()

save(means_tm, file="./Output/TCGA_SKCM_TM_mutsigs_means.RData")



#--------- PRIMARY VS. METASTATIC -------------------------------------------------------------------------------------

load("./Output/TCGA_SKCM_PT_weights_trim_cut0.00.RData")
load("./Output/TCGA_SKCM_TM_weights_trim_cut0.00.RData")

pt <- weight_pt[,which(colSums(weight_pt)>0)]
tm <- weight_tm[,which(colSums(weight_tm)>0)]
t.test(pt$S1, pt$S7)
t.test(tm$S1, tm$S7)

# Boxplot S1 in both samples
pt_s1 <- data.frame(Type="PT", value=pt$S1*100)
tm_s1 <- data.frame(Type="TM", value=tm$S1*100)
S1 <- rbind(pt_s1, tm_s1)

colors3 <- brewer.pal(3, "Paired")
pdf("./Figures/TCGA_SKCM_S1_mutsigs_boxplot.pdf", w=4, h=5)
ggplot(S1, aes(x=Type, y=value, fill=Type)) +
    geom_boxplot() +
    scale_fill_manual(values=colors3) +
    ggtitle("SKCM S1") +
    xlab("Sample type") +
    ylab("Contribution (%)")
dev.off()

# Boxplot S7 in both samples
pt_s7 <- data.frame(Type="PT", value=pt$S7*100)
tm_s7 <- data.frame(Type="TM", value=tm$S7*100)
S7 <- rbind(pt_s7, tm_s7)

colors3 <- brewer.pal(3, "Paired")
pdf("./Figures/TCGA_SKCM_S7_mutsigs_boxplot.pdf", w=4, h=5)
ggplot(S7, aes(x=Type, y=value, fill=Type)) +
    geom_boxplot() +
    scale_fill_manual(values=colors3) +
    ggtitle("SKCM S7") +
    xlab("Sample type") +
    ylab("Contribution (%)")
dev.off()

t.test(pt$S1, tm$S1)
t.test(pt$S7, tm$S7)


pt$patient_id <- sapply(rownames(pt), function(x) as.factor(tolower(strsplit(x,"-")[[1]][3])))
tm$patient_id <- sapply(rownames(tm), function(x) as.factor(tolower(strsplit(x,"-")[[1]][3])))

matched <- merge(pt, tm, by="patient_id")       #none


save(pt, file="./Output/TCGA_SKCM_PT_weights_trim_cut0.00.RData")
save(tm, file="./Output/TCGA_SKCM_TM_weights_trim_cut0.00.RData")




#--------- DOMINANT SIGNATURE -------------------------------------------------------------------------------------

# Extract dominant signature per patient in PT
pt.max <- pt[,-25]
pt.max$max <- sapply(1:nrow(pt.max), function(x) max(pt.max[x, 1:24]))
pt.max$maxSig <- sapply(1:nrow(pt.max) , function(x) noquote(colnames(pt.max[,which(pt.max[x,]==pt.max$max[x])])[1]))

# Extract dominant signature per patient in TM
tm.max <- tm[,-32]
tm.max$max <- sapply(1:nrow(tm.max), function(x) max(tm.max[x, 1:31]))
tm.max$maxSig <- sapply(1:nrow(tm.max) , function(x) noquote(colnames(tm.max[,which(tm.max[x,]==tm.max$max[x])])[1]))

# Generate data frame of dominant signatures per cohorts
pt.domsigs <- as.data.frame(table(pt.max$maxSig))
tm.domsigs <- as.data.frame(table(tm.max$maxSig))
domsig <- merge(pt.domsigs, tm.domsigs, by="Var1", all=TRUE)
domsig[is.na(domsig)] <- 0
colnames(domsig) <- c("Signature", "SKCM_PT", "SKCM_TM")
domsig <- domsig[match(c("S1","S3","S6","S7","S15","S16","S18","S20","S21","S26"),domsig$Signature),]

# Generates contingency tables and runs Fisher's exact test --> outputs into data frame
fisher <- data.frame(Signature=domsig$Signature, OR=NA, pval=NA, padj=NA)
for (i in 1:nrow(domsig)) {
    p <- domsig[i,2]
    m <- domsig[i,3]
    sump <- colSums(domsig[2])-p
    summ <- colSums(domsig[3])-m
    mat <- array(c(p, sump, m, summ), c(2,2))
    fisher[i,2] <- fisher.test(mat)$estimate
    fisher[i,3] <- fisher.test(mat)$p.value
}
fisher$padj <- p.adjust(fisher$pval, method="BH")  
write.csv(fisher, file="./Output/TCGA_SKCM_fisher.csv")   
    
    
    
    
#--------- SKCM_PT, SKCM_TM, UVM --------------------------------------------------------------------------------------

load("./Output/TCGA_SKCM_PT_mutsigs_means.RData")   # means_pt
load("./Output/TCGA_SKCM_TM_mutsigs_means.RData")   # means_tm
load("./Output/TCGA_UVM_PT_mutsigs_means.RData")    # means

colnames(means_pt)[2] <- "SKCM_PT"
colnames(means_tm)[2] <- "SKCM_TM"
colnames(means)[2] <- "UVM_PT"

merged <- merge(means_pt, means_tm, by="Signatures", all=TRUE)
merged <- merge(merged, means, by="Signatures", all=TRUE)
merged[is.na(merged)] <- 0

# Ordered by signature number
order <- c(colnames(weight_pt))
merged <- merged[match(order,merged$Signatures),]

melt_merged <- melt(merged)
colnames(melt_merged) <- c("Signatures", "Sample", "value")


colors3 <- brewer.pal(3, "Paired")
pdf("./Figures/TCGA_SKCM_UVM_mutsigs_avg.pdf", w=13, h=6)
ggplot(melt_merged, aes(x=Signatures, y=value)) + 
    scale_fill_manual(values=colors3) +
    geom_bar(stat="identity", aes(fill=Sample), position="dodge") +
    ggtitle("Mean mutational signature contributions SKCM and UVM") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5))
dev.off()




#--------- CLINICAL DATA ----------------------------------------------------------------------------------------------

load("./Data/TCGA_SKCM_clinical.RData")
load("./Output/TCGA_SKCM_PT_clinical.RData")
load('./Output/TCGA_SKCM_TM_clinical.RData')

# Download data:
SKCM.clin = Samples.Clinical(format = "csv",
                           cohort = "SKCM",
                           page_size=2000)
save(SKCM.clin, file="./Data/TCGA_SKCM_clinical.RData")

# Create variables with both clinical and signatures data
pt.clin <- merge(SKCM.clin, pt, by="patient_id")
tm.clin <- merge(SKCM.clin, tm, by="patient_id")
save(pt.clin, file="./Output/TCGA_SKCM_PT_clinical.RData")
save(tm.clin, file="./Output/TCGA_SKCM_TM_clinical.RData")

#--------- ~ Stage PT ---------

load("./Output/TCGA_SKCM_PT_stage.RData")

# Get just id, cancer stage, and sigs
pt.stage <- pt.clin[, c(1:2, 47, 71:94)]
pt.stage$stage <- sapply(pt.stage$pathologic_stage, function(x) toupper(strsplit(x,"\\ ")[[1]][2]))
pt.stage$stageNum <- pt.stage$stage

# Loop to remove a, b, c specifications of cancer stage
for (i in 1:nrow(pt.stage)) {
    if (grepl("a", pt.stage$stage[i], ignore.case = TRUE) | grepl("b", pt.stage$stage[i], ignore.case = TRUE) | 
        grepl("c", pt.stage$stage[i], ignore.case = TRUE)) {
        pt.stage$stage[i] <- substr(pt.stage$stage[i], 1, nchar(pt.stage$stage[i])-1)
    }
    else {
        pt.stage$stage[i] <- pt.stage$stage[i]
    }
}

# Add column for early/late stage
pt.stage$earlyLate <- pt.stage$stage
for (i in 1:nrow(pt.stage)) {
    if (is.na(pt.stage$earlyLate[i])){
        pt.stage$earlyLate[i] <- NA
    }
    else if (pt.stage$earlyLate[i] == "I" | pt.stage$earlyLate[i] == "II") {
        pt.stage$earlyLate[i] <- "early"
    }
    else if (pt.stage$earlyLate[i] == "III" | pt.stage$earlyLate[i] == "IV") {
        pt.stage$earlyLate[i] <- "late"
    }
}

# Loop for numeric version of stage
for (i in 1:nrow(pt.stage)){
    if (is.na(pt.stage$stageNum[i])){
        pt.stage$stageNum[i] <- 0
    }
    else if (pt.stage$stageNum[i] == "I") {
        pt.stage$stageNum[i] <- 1
    }
    else if (pt.stage$stageNum[i] == "II") {
        pt.stage$stageNum[i] <- 2
    }
    else if (pt.stage$stageNum[i] == "III") {
        pt.stage$stageNum[i] <- 3
    }
    else if (pt.stage$stageNum[i] == "IV") {
        pt.stage$stageNum[i] <- 4
    }
}
pt.stage$stage2 <- pt.stage$stage
pt.stage <- pt.stage[,-28]



# Remove patient with pathological stage = i/ii nos (1/104)
pt.stage <- pt.stage[-which(pt.stage$stage == "NOS"),]

save(pt.stage, file="./Output/TCGA_SKCM_PT_stage.RData")


## Overview
pt.stage.melt <- melt(pt.stage)[,c(2,4:6)]
colnames(pt.stage.melt) <- c("PatientId", "Stage", "Signatures", "weight")
pt.stage.melt$weight <- sapply(pt.stage.melt$weight, function(x) 100*x)

pt.early.late.melt <- melt(pt.stage[,c(2,4:27,29)])
colnames(pt.early.late.melt) <- c("PatientId", "earlyLate", "Signatures", "weight")
pt.early.late.melt$weight <- sapply(pt.early.late.melt$weight, function(x) 100*x)


# Barplot by stage
colors4 <- c(brewer.pal(4, "Paired"),"#000000")
pdf("./Figures/TCGA_SKCM_PT_mutsigs_stage_barplot.pdf", w=13, h=6)
ggplot(pt.stage.melt, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors4) +
    geom_bar(stat="identity", aes(fill=Stage), position="dodge") +
    ggtitle("Mutational signature contributions SKCM PT by cancer stage") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian( ylim=c(0,100), expand = FALSE )
dev.off()


# Mean sig barplot per stage
# Extract mean sig per cancer stage
pt.stagex <- pt.stage[,c(4:28)]
rownames(pt.stagex) <- pt.stage[,2]
pt.stage.mean <- sapply(pt.stagex[, 1:24], function(x) tapply(x, pt.stagex[, 25], mean))
pt.stage.mean <- melt(pt.stage.mean)
colnames(pt.stage.mean) <- c("Stage", "Signatures", "weight")
pt.stage.mean$weight <- sapply(pt.stage.mean$weight, function(x) 100*x)

# Plot barplot
pdf("./Figures/TCGA_SKCM_PT_mutsigs_stage_avgbarplot.pdf", w=13, h=6)
ggplot(pt.stage.mean, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors4) +
    geom_bar(stat="identity", aes(fill=Stage), position="dodge") +
    ggtitle("Mean mutational signature contributions SKCM PT by cancer stage") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,60), expand = FALSE)
dev.off()

# Mean sig barplot per progression
pt.stagep <- pt.stage[,c(4:27, 29)]
rownames(pt.stagep) <- pt.stage[,2]
pt.stagep.mean <- sapply(pt.stagep[, 1:24], function(x) tapply(x, pt.stagep[, 25], mean))
pt.stagep.mean <- melt(pt.stagep.mean)
colnames(pt.stagep.mean) <- c("Progression", "Signatures", "weight")
pt.stagep.mean$weight <- sapply(pt.stagep.mean$weight, function(x) 100*x) 

# Plot barplot
pdf("./Figures/TCGA_SKCM_PT_mutsigs_prog_avgbarplot.pdf", w=13, h=6)
ggplot(pt.stagep.mean, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors4) +
    geom_bar(stat="identity", aes(fill=Progression), position="dodge") +
    ggtitle("Mean mutational signature contributions SKCM PT by disease progression") +
    xlab("Signaturtes") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,60), expand = FALSE)
dev.off()


# Loop over each signature, by stage
# Data frame for t-test output
pt.prog.ttest <- data.frame(Signature=colnames(pt.stage)[4:27], t=NA, df=NA, pval = NA)
for (i in 4:27) {
    sig <- noquote(colnames(pt.stage[i]))
    var <- na.omit(data.frame(Progression=pt.stage$earlyLate, weight=pt.stage[,i]))
    title <- paste0("SKCM PT ", sig, " contribution by disease progression")
    file <- paste0("./Figures/TCGA_SKCM_PT_prog_", sig, ".pdf")
    
    early <- var[which(var$Progression == "early"),]
    late <- var[which(var$Progression == "late"),]
    ttest <- t.test(early$weight, late$weight)
    
    pt.prog.ttest[which(pt.prog.ttest$Signature==sig),2] <- ttest$statistic
    pt.prog.ttest[which(pt.prog.ttest$Signature==sig),3] <- ttest$parameter
    pt.prog.ttest[which(pt.prog.ttest$Signature==sig),4] <- ttest$p.value
    #pt.prog.ttest[which(pt.prog.ttest$Signature==sig),5] <- ttest$conf.int
    #pt.prog.ttest[which(pt.prog.ttest$Signature==sig),6] <- ttest$estimate
    
    pdf(file, w=7, h=5)
    plot <- ggplot(var, aes(x=Progression, y=weight*100)) +
         #geom_bar(stat="identity", fill="blue") +
         geom_boxplot() +
         ggtitle(title) +
         xlab("Disease progression") +
         ylab("Contribution (%)")
    print(plot)
    dev.off()
}
pt.prog.ttest$padj <- p.adjust(pt.prog.ttest$pval, method="BH")

write.csv(pt.prog.ttest, file="./Output/TCGA_SKCM_PT_prog_eachsigttest.csv")



# Group by stage (4/104 = NA)
pt.stage1 <- pt.stage[which(pt.stage$stage == "I"),c(2,4:27)]
rownames(pt.stage1) <- pt.stage1$tcga_participant_barcode
pt.stage1 <- pt.stage1[,-1]

pt.stage2 <- pt.stage[which(pt.stage$stage == "II"),c(2,4:27)]
rownames(pt.stage2) <- pt.stage2$tcga_participant_barcode
pt.stage2 <- pt.stage2[,-1]

pt.stage3 <- pt.stage[which(pt.stage$stage == "III"),c(2,4:27)]
rownames(pt.stage3) <- pt.stage3$tcga_participant_barcode
pt.stage3 <- pt.stage3[,-1]

pt.stage4 <- pt.stage[which(pt.stage$stage == "IV"),c(2,4:27)]
rownames(pt.stage4) <- pt.stage4$tcga_participant_barcode
pt.stage4 <- pt.stage4[,-1]


## Stage I
pt.stage1 <- as.matrix(pt.stage1)
pt.stage1.trim <- pt.stage1[,which(colSums(pt.stage1)>0)]
pt.stage1.melt <- melt(pt.stage1.trim)
colnames(pt.stage1.melt) <- c("PatientId", "Signature", "weight")
pt.stage1.melt$weight <- sapply(pt.stage1.melt$weight, function(x) 100*x)

#Boxplot
pdf("./Figures/TCGA_SKCM_PT_mutsigs_stage1_boxplot.pdf", w=10, h=6)
ggplot(pt.stage1.melt, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM PT Stage I") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

#Anova
pt.stage1.anova <- aov(weight~Signature,  data=pt.stage1.melt)
summary(pt.stage1.anova)


## Stage II
pt.stage2 <- as.matrix(pt.stage2)
pt.stage2.trim <- pt.stage2[,which(colSums(pt.stage2)>0)]
pt.stage2.melt <- melt(pt.stage2.trim)
colnames(pt.stage2.melt) <- c("PatientId", "Signature", "weight")
pt.stage2.melt$weight <- sapply(pt.stage2.melt$weight, function(x) 100*x)

#Boxplot
pdf("./Figures/TCGA_SKCM_PT_mutsigs_stage2_boxplot.pdf", w=10, h=6)
ggplot(pt.stage2.melt, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM PT Stage II") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

#Anova
pt.stage2.anova <- aov(weight~Signature,  data=pt.stage2.melt)
summary(pt.stage2.anova)


## Stage III
pt.stage3 <- as.matrix(pt.stage3)
pt.stage3.trim <- pt.stage3[,which(colSums(pt.stage3)>0)]
pt.stage3.melt <- melt(pt.stage3.trim)
colnames(pt.stage3.melt) <- c("PatientId", "Signature", "weight")
pt.stage3.melt$weight <- sapply(pt.stage3.melt$weight, function(x) 100*x)

pdf("./Figures/TCGA_SKCM_PT_mutsigs_stage3_boxplot.pdf", w=10, h=6)
ggplot(pt.stage3.melt, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM PT Stage III") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

#Anova
pt.stage3.anova <- aov(weight~Signature,  data=pt.stage3.melt)
summary(pt.stage3.anova)


## Stage IV
pt.stage4 <- as.matrix(pt.stage4)
pt.stage4.trim <- pt.stage4[,which(colSums(pt.stage4)>0)]
pt.stage4.melt <- melt(pt.stage4.trim)
colnames(pt.stage4.melt) <- c("PatientId", "Signature", "weight")
pt.stage4.melt$weight <- sapply(pt.stage4.melt$weight, function(x) 100*x)

#Boxplot
pdf("./Figures/TCGA_SKCM_PT_mutsigs_stage4_boxplot.pdf", w=10, h=6)
ggplot(pt.stage4.melt, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM PT Stage IV") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

#Anova
pt.stage4.anova <- aov(weight~Signature,  data=pt.stage4.melt)
summary(pt.stage4.anova)



#--------- ~ Stage TM ---------

load("./Output/TCGA_SKCM_TM_stage.RData")

# Get just id, cancer stage, and sigs
tm.stage <- tm.clin[, c(1:2, 47, 71:101)]
tm.stage$stage <- sapply(tm.stage$pathologic_stage, function(x) toupper(strsplit(x,"\\ ")[[1]][2]))

# Loop to remove a, b, c specifications of cancer stage
for (i in 1:nrow(tm.stage)) {
    if (grepl("a", tm.stage$stage[i], ignore.case = TRUE) | grepl("b", tm.stage$stage[i], ignore.case = TRUE) | 
        grepl("c", tm.stage$stage[i], ignore.case = TRUE)) {
        tm.stage$stage[i] <- substr(tm.stage$stage[i], 1, nchar(tm.stage$stage[i])-1)
    }
    else {
        tm.stage$stage[i] <- tm.stage$stage[i]
    }
}

# Add column for early/late stage
tm.stage$earlyLate <- tm.stage$stage
for (i in 1:nrow(tm.stage)) {
    if (is.na(tm.stage$earlyLate[i])){
        tm.stage$earlyLate[i] <- NA
    }
    else if (tm.stage$earlyLate[i] == 0 | tm.stage$earlyLate[i] == "I" | tm.stage$earlyLate[i] == "II") {
        tm.stage$earlyLate[i] <- "early"
    }
    else if (tm.stage$earlyLate[i] == "III" | tm.stage$earlyLate[i] == "IV") {
        tm.stage$earlyLate[i] <- "late"
    }
}


# Remove patients with pathological stage = i/ii nos (13/104)
tm.stage <- tm.stage[-which(tm.stage$stage == "NOS"),]

save(tm.stage, file="./Output/TCGA_SKCM_TM_stage.RData")


## Overview
tm.stage.melt <- melt(tm.stage[,c(2, 4:35)])
colnames(tm.stage.melt) <- c("PatientId", "Stage", "Signatures", "weight")
tm.stage.melt$weight <- sapply(tm.stage.melt$weight, function(x) 100*x)

tm.early.late.melt <- melt(tm.stage[,c(2,4:34,36)])
colnames(tm.early.late.melt) <- c("PatientId", "earlyLate", "Signatures", "weight")
tm.early.late.melt$weight <- sapply(tm.early.late.melt$weight, function(x) 100*x)


# Barplot by stage
colors6 <- c("#fb9a99", brewer.pal(4, "Paired"),"#000000")
pdf("./Figures/TCGA_SKCM_TM_mutsigs_stage_barplot.pdf", w=13, h=6)
ggplot(tm.stage.melt, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors6) +
    geom_bar(stat="identity", aes(fill=Stage), position="dodge") +
    ggtitle("Mutational signature contributions SKCM TM by cancer stage") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian( ylim=c(0,100), expand = FALSE )
dev.off()


# Mean sig barplot per stage
# Extract mean sig per cancer stage
tm.stagex <- tm.stage[,c(4:35)]
rownames(tm.stagex) <- tm.stage[,2]
tm.stage.mean <- sapply(tm.stagex[, 1:31], function(x) tapply(x, tm.stagex[, 32], mean))
tm.stage.mean <- melt(tm.stage.mean)
colnames(tm.stage.mean) <- c("Stage", "Signatures", "weight")
tm.stage.mean$weight <- sapply(tm.stage.mean$weight, function(x) 100*x)

# Plot barplot
pdf("./Figures/TCGA_SKCM_TM_mutsigs_stage_avgbarplot.pdf", w=13, h=6)
ggplot(tm.stage.mean, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors6) +
    geom_bar(stat="identity", aes(fill=Stage), position="dodge") +
    ggtitle("Mean mutational signature contributions SKCM TM by cancer stage") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,60), expand = FALSE)
dev.off()

# Mean sig barplot per progression
tm.stagep <- tm.stage[,c(4:34, 36)]
rownames(tm.stagep) <- tm.stage[,2]
tm.stagep.mean <- sapply(tm.stagep[, 1:31], function(x) tapply(x, tm.stagep[, 32], mean))
tm.stagep.mean <- melt(tm.stagep.mean)
colnames(tm.stagep.mean) <- c("Progression", "Signatures", "weight")
tm.stagep.mean$weight <- sapply(tm.stagep.mean$weight, function(x) 100*x) 

# Plot barplot
pdf("./Figures/TCGA_SKCM_TM_mutsigs_prog_avgbarplot.pdf", w=13, h=6)
ggplot(tm.stagep.mean, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors4) +
    geom_bar(stat="identity", aes(fill=Progression), position="dodge") +
    ggtitle("Mean mutational signature contributions SKCM TM by disease progression") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,60), expand = FALSE)
dev.off()


# Loop over each signature, by stage
# Data frame for t-test output
tm.prog.ttest <- data.frame(Signature=colnames(tm.stage)[4:34], t=NA, df=NA, pval = NA)
for (i in 4:34) {
    sig <- noquote(colnames(tm.stage[i]))
    var <- na.omit(data.frame(Progression=tm.stage$earlyLate, weight=tm.stage[,i]))
    title <- paste0("SKCM TM ", sig, " contribution by disease progression")
    file <- paste0("./Figures/TCGA_SKCM_TM_prog_", sig, ".pdf")
    
    early <- var[which(var$Progression == "early"),]
    late <- var[which(var$Progression == "late"),]
    ttest <- t.test(early$weight, late$weight)
    
    tm.prog.ttest[which(tm.prog.ttest$Signature==sig),2] <- ttest$statistic
    tm.prog.ttest[which(tm.prog.ttest$Signature==sig),3] <- ttest$parameter
    tm.prog.ttest[which(tm.prog.ttest$Signature==sig),4] <- ttest$p.value
    #tm.prog.ttest[which(tm.prog.ttest$Signature==sig),5] <- ttest$conf.int
    #tm.prog.ttest[which(tm.prog.ttest$Signature==sig),6] <- ttest$estimate
    
    pdf(file, w=7, h=5)
    plot <- ggplot(var, aes(x=Progression, y=weight*100)) +
        #geom_bar(stat="identity", fill="blue") +
        geom_boxplot() +
        ggtitle(title) +
        xlab("Disease progression") +
        ylab("Contribution (%)")
    print(plot) 
    dev.off()
}
tm.prog.ttest$padj <- p.adjust(tm.prog.ttest$pval, method="BH")

write.csv(tm.prog.ttest, file="./Output/TCGA_SKCM_TM_prog_eachsigttest.csv")



#--------- ~ Therapies ---------
# Get list of annotated therapies
therapies <- array()
for (i in 1:ncol(SKCM.clin)) {
    if (grepl("therapy", colnames(SKCM.clin[i])) | grepl("treatment", colnames(SKCM.clin[i])) | grepl("drug", colnames(SKCM.clin[i]))) {
        therapies <- append(therapies, colnames(SKCM.clin[i]))
    }
}
therapies

# Summary
table(pt.clin$prior_systemic_therapy_type, useNA="always")
table(pt.clin$radiation_therapy, useNA="always")
table(pt.clin$history_of_neoadjuvant_treatment, useNA="always")

table(tm.clin$prior_systemic_therapy_type, useNA="always")
table(tm.clin$radiation_therapy, useNA="always").
table(tm.clin$history_of_neoadjuvant_treatment, useNA="always")

# Therapy subgroups



### TO ASK
# Therapy annotations:  interferon 90 day prior excision admin indicator -- immunotherapy before surgery
#                       history of neoadjuvant treatment (either no or NA)
#                               = chemo after surgery




# cancer stage barplot: use avg --> complementary to boxplot sigs analysis

