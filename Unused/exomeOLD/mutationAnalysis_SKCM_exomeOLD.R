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
load("./Output/TCGA_SKCM_TP_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_SKCM_TP_weights_cut0.00.RData")     # Mutational signatures output

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
snvs.tp <- snvs[which(snvs$Sample_Type == "01"),]           # 61990 samples
snvs.tm <- snvs[which(snvs$Sample_Type == "06"),]          # 326634 samples




#--------- MUTATIONAL SIGNATURES ANALYSIS: TP -------------------------------------------------------------------------

# Convert to deconstructSigs input:
sigs.input.tp <- mut.to.sigs.input(mut.ref = snvs.tp, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

save(sigs.input.tp, file="./Output/TCGA_SKCM_TP_mutsig_input.RData")


# Signatures
weights.tp <- as.data.frame(t(sapply(rownames(sigs.input.tp), 
                                    function(x) whichSignatures(tumor.ref = sigs.input.tp,
                                                                signatures.ref = signatures.cosmic,
                                                                sample.id = x, 
                                                                contexts.needed = TRUE,
                                                                signature.cutoff = 0.00,              # default = 0.06
                                                                tri.counts.method = 'exome')$weights)), 
                           row.names=rownames(sigs.input.tp))

weight.tp <- as.data.frame(sapply(weights.tp, function(x) unlist(x)), row.names = rownames(weights.tp))

# Change Signature.# to S#
colnames(weight.tp) <- sapply(colnames(weight.tp), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))

# Add unknown contribution column
weight.tp$unknown <- (1-rowSums(weight.tp))

save(weight.tp, file="./Output/TCGA_SKCM_TP_weights_cut0.00.RData")



# Sort by S1 (descending)
sorted.tp <- weight.tp[order(-weight.tp$S1),]


# Reformatting the weights data 
weight.tp <- as.matrix(weight.tp)
melted.tp <- melt(weight.tp)
colnames(melted.tp) <- c("PatientId", "Signature", "weight")
melted.tp$weight <- sapply(melted.tp$weight, function(x) 100*x)


# Remove non-contributing signatures
sorted.tp <- as.matrix(sorted.tp)
trim.tp <- sorted.tp[,which(colSums(sorted.tp)>0)]
melt.tp <- melt(trim.tp)
colnames(melt.tp) <- c("PatientId", "Signature", "weight")
melt.tp$weight <- sapply(melt.tp$weight, function(x) 100*x)


# Plot: stacked barplot signatures
colors <- c(brewer.pal(12, "Paired"), brewer.pal(11, "Set3"),"#000000")
pdf("./Figures/TCGA_SKCM_TP_mutsig.pdf", w=10, h=6)
ggplot(melt.tp, aes(x=PatientId, y=weight, fill=Signature)) + 
    scale_fill_manual(values=colors) +
    geom_bar(stat="identity") +
    ggtitle("Mutational signature analysis SKCM primary tumor ") +
    xlab("Sample") +
    ylab("Contribution (%)") +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
dev.off()
    
# Plot: boxplot signatures
pdf("./Figures/TCGA_SKCM_TP_mutsigs_boxplot.pdf", w=10, h=6)
ggplot(melt.tp, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM primary tumor") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

# Plot: heatmap
melted.tp$weight <- as.numeric(melted.tp$weight)
tp.mat <- acast(PatientId~Signature,data=melted.tp, value="weight", fun.aggregate=mean)
pdf("./Figures/TCGA_SKCM_TP_mutsigs_heatmap.pdf", w=10, h=6, onefile=FALSE)
p.tp <- pheatmap(t(tp.mat), show_colnames = FALSE, main="Mutational signatures SKCM primary tumor")
dev.off()


#--------- ~ S1 and S7 analysis: TP ---------
sp1_7 <- as.matrix(weight.tp[,c(1, 7)])
melt_sp1_7 <- melt(sp1_7)
colnames(melt_sp1_7) <- c("PatientId", "Signature", "weight")
melt_sp1_7$weight <- sapply(melt_sp1_7$weight, function(x) 100*x)

# Plot: boxplot signatures
pdf("./Figures/TCGA_SKCM_TP_mutsigs1_7_boxplot.pdf", w=4, h=5)
ggplot(melt_sp1_7, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("SKCM primary tumor S1 vs S7") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()


#--------- ~ Mean contributions per signature ---------

means.tp <- data.frame(Signatures=c(colnames(trim.tp)), value=c(colMeans(trim.tp)*100))

pdf("./Figures/TCGA_SKCM_TP_mutsigs_avg.pdf", w=10, h=6)
ggplot(means.tp, aes(x=reorder(Signatures, -value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    ggtitle("Mean mutational signature contributions SKCM primary tumor") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5))
dev.off()

save(means.tp, file="./Output/TCGA_SKCM_TP_mutsigs_means.RData")



#--------- MUTATIONAL SIGNATURES ANALYSIS: TM -------------------------------------------------------------------------

# Convert to deconstructSigs input:
sigs_input.tm <- mut.to.sigs.input(mut.ref = snvs_met, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

save(sigs_input.tm, file="./Output/TCGA_SKCM_TM_mutsig_input.RData")


# Signatures
weights.tm <- as.data.frame(t(sapply(rownames(sigs_input.tm), 
                                     function(x) whichSignatures(tumor.ref = sigs_input.tm,
                                                                 signatures.ref = signatures.cosmic,
                                                                 sample.id = x, 
                                                                 contexts.needed = TRUE,
                                                                 signature.cutoff = 0.00,              # default = 0.06
                                                                 tri.counts.method = 'exome')$weights)), 
                            row.names=rownames(sigs_input.tm))

weight.tm <- as.data.frame(sapply(weights.tm, function(x) unlist(x)), row.names = rownames(weights.tm))

# Change Signature.# to S#
colnames(weight.tm) <- sapply(colnames(weight.tm), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))

# Add unknown contribution column
weight.tm$unknown <- (1-rowSums(weight.tm))

save(weight.tm, file="./Output/TCGA_SKCM_TM_weights_cut0.00.RData")



# Sort by S1 (descending)
sorted.tm <- weight.tm[order(-weight.tm$S1),]


# Reformatting the weights data 
weight.tm <- as.matrix(weight.tm)
melted.tm <- melt(weight.tm)
colnames(melted.tm) <- c("PatientId", "Signature", "weight")
melted.tm$weight <- sapply(melted.tm$weight, function(x) 100*x)


# Remove non-contributing signatures
sorted.tm <- as.matrix(sorted.tm)
trim.tm <- sorted.tm[,which(colSums(sorted.tm)>0)]
melt.tm <- melt(trim.tm)
colnames(melt.tm) <- c("PatientId", "Signature", "weight")
melt.tm$weight <- sapply(melt.tm$weight, function(x) 100*x)


colors2 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#BFA89E", "#826C7F", "#FF7F00",
             "#CAB2D6", "#8B786D", "#6A3D9A", "#D6DBD2", "#FFFF99", "#B15928", "#8DD3C7", "#FFFFB3", "#E75A7C", "#BEBADA",
             "#FB8072", "#40476D", "#80B1D3", "#FDB462", "#B3DE69", "#258EA6", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", 
             "#000000")
# Plot: stacked barplot signatures
pdf("./Figures/TCGA_SKCM_TM_mutsig.pdf", w=10, h=6)
ggplot(melt.tm, aes(x=PatientId, y=weight, fill=Signature)) + 
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
ggplot(melt.tm, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM metastatic tumor") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

# Plot: heatmap
melted.tm$weight <- as.numeric(melted.tm$weight)
tm.mat <- acast(PatientId~Signature, data=melted.tm, value="weight", fun.aggregate=mean)
pdf("./Figures/TCGA_SKCM_TM_mutsigs_heatmap.pdf", w=10, h=6, onefile=FALSE)
p.tm <- pheatmap(t(tm.mat), show_colnames = FALSE, main="Mutational signatures SKCM metastatic tumor")
dev.off()


#--------- ~ S1 and S7 analysis: TM ---------
sm1_7 <- as.matrix(weight.tm[,c(1, 7)])
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

means.tm <- data.frame(Signatures=c(colnames(trim.tm)), value=c(colMeans(trim.tm)*100))

pdf("./Figures/TCGA_SKCM_TM_mutsigs_avg.pdf", w=13, h=6)
ggplot(means.tm, aes(x=reorder(Signatures, -value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    ggtitle("Mean mutational signature contributions SKCM metastatic tumor") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5))
dev.off()

save(means.tm, file="./Output/TCGA_SKCM_TM_mutsigs_means.RData")



#--------- PRIMARY VS. METASTATIC -------------------------------------------------------------------------------------

load("./Output/TCGA_SKCM_TP_weights_trim_cut0.00.RData")
load("./Output/TCGA_SKCM_TM_weights_trim_cut0.00.RData")

tp <- weight.tp[,which(colSums(weight.tp)>0)]
tm <- weight.tm[,which(colSums(weight.tm)>0)]
t.test(tp$S1, tp$S7)
t.test(tm$S1, tm$S7)

# Boxplot S1 in both samples
tp.s1 <- data.frame(Type="TP", value=tp$S1*100)
tm.s1 <- data.frame(Type="TM", value=tm$S1*100)
S1 <- rbind(tp.s1, tm.s1)

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
tp.s7 <- data.frame(Type="TP", value=tp$S7*100)
tm.s7 <- data.frame(Type="TM", value=tm$S7*100)
S7 <- rbind(tp.s7, tm.s7)

colors3 <- brewer.pal(3, "Paired")
pdf("./Figures/TCGA_SKCM_S7_mutsigs_boxplot.pdf", w=4, h=5)
ggplot(S7, aes(x=Type, y=value, fill=Type)) +
    geom_boxplot() +
    scale_fill_manual(values=colors3) +
    ggtitle("SKCM S7") +
    xlab("Sample type") +
    ylab("Contribution (%)")
dev.off()

t.test(tp$S1, tm$S1)
t.test(tp$S7, tm$S7)


tp$patient_id <- sapply(rownames(tp), function(x) as.factor(tolower(strsplit(x,"-")[[1]][3])))
tm$patient_id <- sapply(rownames(tm), function(x) as.factor(tolower(strsplit(x,"-")[[1]][3])))

matched <- merge(tp, tm, by="patient_id")       #none


save(tp, file="./Output/TCGA_SKCM_TP_weights_trim_cut0.00.RData")
save(tm, file="./Output/TCGA_SKCM_TM_weights_trim_cut0.00.RData")




#--------- DOMINANT SIGNATURE -------------------------------------------------------------------------------------

# Extract dominant signature per patient in TP
tp.max <- tp[,-25]
tp.max$max <- sapply(1:nrow(tp.max), function(x) max(tp.max[x, 1:24]))
tp.max$maxSig <- sapply(1:nrow(tp.max) , function(x) noquote(colnames(tp.max[,which(tp.max[x,]==tp.max$max[x])])[1]))

# Extract dominant signature per patient in TM
tm.max <- tm[,-32]
tm.max$max <- sapply(1:nrow(tm.max), function(x) max(tm.max[x, 1:31]))
tm.max$maxSig <- sapply(1:nrow(tm.max) , function(x) noquote(colnames(tm.max[,which(tm.max[x,]==tm.max$max[x])])[1]))

# Generate data frame of dominant signatures per cohorts
tp.domsigs <- as.data.frame(table(tp.max$maxSig))
tm.domsigs <- as.data.frame(table(tm.max$maxSig))
domsig <- merge(tp.domsigs, tm.domsigs, by="Var1", all=TRUE)
domsig[is.na(domsig)] <- 0
colnames(domsig) <- c("Signature", "SKCM_TP", "SKCM_TM")
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
    
    
    
    
#--------- SKCM_TP, SKCM_TM, UVM --------------------------------------------------------------------------------------

load("./Output/TCGA_SKCM_TP_mutsigs_means.RData")   # means.tp
load("./Output/TCGA_SKCM_TM_mutsigs_means.RData")   # means.tm
load("./Output/TCGA_UVM_TP_mutsigs_means.RData")    # means

colnames(means.tp)[2] <- "SKCM_TP"
colnames(means.tm)[2] <- "SKCM_TM"
colnames(means)[2] <- "UVM_TP"

merged <- merge(means.tp, means.tm, by="Signatures", all=TRUE)
merged <- merge(merged, means, by="Signatures", all=TRUE)
merged[is.na(merged)] <- 0

# Ordered by signature number
order <- c(colnames(weight.tp))
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
load("./Output/TCGA_SKCM_TP_clinical.RData")
load('./Output/TCGA_SKCM_TM_clinical.RData')

# Download data:
SKCM.clin = Samples.Clinical(format = "csv",
                           cohort = "SKCM",
                           page_size=2000)
save(SKCM.clin, file="./Data/TCGA_SKCM_clinical.RData")

# Create variables with both clinical and signatures data
tp.clin <- merge(SKCM.clin, tp, by="patient_id")
tm.clin <- merge(SKCM.clin, tm, by="patient_id")
save(tp.clin, file="./Output/TCGA_SKCM_TP_clinical.RData")
save(tm.clin, file="./Output/TCGA_SKCM_TM_clinical.RData")

#--------- ~ Stage TP ---------

load("./Output/TCGA_SKCM_TP_stage.RData")

# Get just id, cancer stage, and sigs
tp.stage <- tp.clin[, c(1:2, 47, 71:94)]
tp.stage$stage <- sapply(tp.stage$pathologic_stage, function(x) toupper(strsplit(x,"\\ ")[[1]][2]))
tp.stage$stageNum <- tp.stage$stage

# Loop to remove a, b, c specifications of cancer stage
for (i in 1:nrow(tp.stage)) {
    if (grepl("a", tp.stage$stage[i], ignore.case = TRUE) | grepl("b", tp.stage$stage[i], ignore.case = TRUE) | 
        grepl("c", tp.stage$stage[i], ignore.case = TRUE)) {
        tp.stage$stage[i] <- substr(tp.stage$stage[i], 1, nchar(tp.stage$stage[i])-1)
    }
    else {
        tp.stage$stage[i] <- tp.stage$stage[i]
    }
}

# Add column for early/late stage
tp.stage$earlyLate <- tp.stage$stage
for (i in 1:nrow(tp.stage)) {
    if (is.na(tp.stage$earlyLate[i])){
        tp.stage$earlyLate[i] <- NA
    }
    else if (tp.stage$earlyLate[i] == "I" | tp.stage$earlyLate[i] == "II") {
        tp.stage$earlyLate[i] <- "early"
    }
    else if (tp.stage$earlyLate[i] == "III" | tp.stage$earlyLate[i] == "IV") {
        tp.stage$earlyLate[i] <- "late"
    }
}

# Loop for numeric version of stage
for (i in 1:nrow(tp.stage)){
    if (is.na(tp.stage$stageNum[i])){
        tp.stage$stageNum[i] <- 0
    }
    else if (tp.stage$stageNum[i] == "I") {
        tp.stage$stageNum[i] <- 1
    }
    else if (tp.stage$stageNum[i] == "II") {
        tp.stage$stageNum[i] <- 2
    }
    else if (tp.stage$stageNum[i] == "III") {
        tp.stage$stageNum[i] <- 3
    }
    else if (tp.stage$stageNum[i] == "IV") {
        tp.stage$stageNum[i] <- 4
    }
}
tp.stage$stage2 <- tp.stage$stage
tp.stage <- tp.stage[,-28]



# Remove patient with pathological stage = i/ii nos (1/104)
tp.stage <- tp.stage[-which(tp.stage$stage == "NOS"),]

save(tp.stage, file="./Output/TCGA_SKCM_TP_stage.RData")


## Overview
tp.stage.melt <- melt(tp.stage)[,c(2,4:6)]
colnames(tp.stage.melt) <- c("PatientId", "Stage", "Signatures", "weight")
tp.stage.melt$weight <- sapply(tp.stage.melt$weight, function(x) 100*x)

tp.early.late.melt <- melt(tp.stage[,c(2,4:27,29)])
colnames(tp.early.late.melt) <- c("PatientId", "earlyLate", "Signatures", "weight")
tp.early.late.melt$weight <- sapply(tp.early.late.melt$weight, function(x) 100*x)


# Barplot by stage
colors4 <- c(brewer.pal(4, "Paired"),"#000000")
pdf("./Figures/TCGA_SKCM_TP_mutsigs_stage_barplot.pdf", w=13, h=6)
ggplot(tp.stage.melt, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors4) +
    geom_bar(stat="identity", aes(fill=Stage), position="dodge") +
    ggtitle("Mutational signature contributions SKCM TP by cancer stage") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian( ylim=c(0,100), expand = FALSE )
dev.off()


# Mean sig barplot per stage
# Extract mean sig per cancer stage
tp.stagex <- tp.stage[,c(4:28)]
rownames(tp.stagex) <- tp.stage[,2]
tp.stage.mean <- sapply(tp.stagex[, 1:24], function(x) tapply(x, tp.stagex[, 25], mean))
tp.stage.mean <- melt(tp.stage.mean)
colnames(tp.stage.mean) <- c("Stage", "Signatures", "weight")
tp.stage.mean$weight <- sapply(tp.stage.mean$weight, function(x) 100*x)

# Plot barplot
pdf("./Figures/TCGA_SKCM_TP_mutsigs_stage_avgbarplot.pdf", w=13, h=6)
ggplot(tp.stage.mean, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors4) +
    geom_bar(stat="identity", aes(fill=Stage), position="dodge") +
    ggtitle("Mean mutational signature contributions SKCM TP by cancer stage") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,60), expand = FALSE)
dev.off()

# Mean sig barplot per progression
tp.stagep <- tp.stage[,c(4:27, 29)]
rownames(tp.stagep) <- tp.stage[,2]
tp.stagep.mean <- sapply(tp.stagep[, 1:24], function(x) tapply(x, tp.stagep[, 25], mean))
tp.stagep.mean <- melt(tp.stagep.mean)
colnames(tp.stagep.mean) <- c("Progression", "Signatures", "weight")
tp.stagep.mean$weight <- sapply(tp.stagep.mean$weight, function(x) 100*x) 

# Plot barplot
pdf("./Figures/TCGA_SKCM_TP_mutsigs_prog_avgbarplot.pdf", w=13, h=6)
ggplot(tp.stagep.mean, aes(x=Signatures, y=weight)) + 
    scale_fill_manual(values=colors4) +
    geom_bar(stat="identity", aes(fill=Progression), position="dodge") +
    ggtitle("Mean mutational signature contributions SKCM TP by disease progression") +
    xlab("Signaturtes") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,60), expand = FALSE)
dev.off()


# Loop over each signature, by stage
# Data frame for t-test output
tp.prog.ttest <- data.frame(Signature=colnames(tp.stage)[4:27], t=NA, df=NA, pval = NA)
for (i in 4:27) {
    sig <- noquote(colnames(tp.stage[i]))
    var <- na.omit(data.frame(Progression=tp.stage$earlyLate, weight=tp.stage[,i]))
    title <- paste0("SKCM TP ", sig, " contribution by disease progression")
    file <- paste0("./Figures/TCGA_SKCM_TP_prog_", sig, ".pdf")
    
    early <- var[which(var$Progression == "early"),]
    late <- var[which(var$Progression == "late"),]
    ttest <- t.test(early$weight, late$weight)
    
    tp.prog.ttest[which(tp.prog.ttest$Signature==sig),2] <- ttest$statistic
    tp.prog.ttest[which(tp.prog.ttest$Signature==sig),3] <- ttest$parameter
    tp.prog.ttest[which(tp.prog.ttest$Signature==sig),4] <- ttest$p.value
    #tp.prog.ttest[which(tp.prog.ttest$Signature==sig),5] <- ttest$conf.int
    #tp.prog.ttest[which(tp.prog.ttest$Signature==sig),6] <- ttest$estimate
    
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
tp.prog.ttest$padj <- p.adjust(tp.prog.ttest$pval, method="BH")

write.csv(tp.prog.ttest, file="./Output/TCGA_SKCM_TP_prog_eachsigttest.csv")



# Group by stage (4/104 = NA)
tp.stage1 <- tp.stage[which(tp.stage$stage == "I"),c(2,4:27)]
rownames(tp.stage1) <- tp.stage1$tcga_participant_barcode
tp.stage1 <- tp.stage1[,-1]

tp.stage2 <- tp.stage[which(tp.stage$stage == "II"),c(2,4:27)]
rownames(tp.stage2) <- tp.stage2$tcga_participant_barcode
tp.stage2 <- tp.stage2[,-1]

tp.stage3 <- tp.stage[which(tp.stage$stage == "III"),c(2,4:27)]
rownames(tp.stage3) <- tp.stage3$tcga_participant_barcode
tp.stage3 <- tp.stage3[,-1]

tp.stage4 <- tp.stage[which(tp.stage$stage == "IV"),c(2,4:27)]
rownames(tp.stage4) <- tp.stage4$tcga_participant_barcode
tp.stage4 <- tp.stage4[,-1]


## Stage I
tp.stage1 <- as.matrix(tp.stage1)
tp.stage1.trim <- tp.stage1[,which(colSums(tp.stage1)>0)]
tp.stage1.melt <- melt(tp.stage1.trim)
colnames(tp.stage1.melt) <- c("PatientId", "Signature", "weight")
tp.stage1.melt$weight <- sapply(tp.stage1.melt$weight, function(x) 100*x)

#Boxplot
pdf("./Figures/TCGA_SKCM_TP_mutsigs_stage1_boxplot.pdf", w=10, h=6)
ggplot(tp.stage1.melt, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM TP Stage I") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

#Anova
tp.stage1.anova <- aov(weight~Signature,  data=tp.stage1.melt)
summary(tp.stage1.anova)


## Stage II
tp.stage2 <- as.matrix(tp.stage2)
tp.stage2.trim <- tp.stage2[,which(colSums(tp.stage2)>0)]
tp.stage2.melt <- melt(tp.stage2.trim)
colnames(tp.stage2.melt) <- c("PatientId", "Signature", "weight")
tp.stage2.melt$weight <- sapply(tp.stage2.melt$weight, function(x) 100*x)

#Boxplot
pdf("./Figures/TCGA_SKCM_TP_mutsigs_stage2_boxplot.pdf", w=10, h=6)
ggplot(tp.stage2.melt, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM TP Stage II") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

#Anova
tp.stage2.anova <- aov(weight~Signature,  data=tp.stage2.melt)
summary(tp.stage2.anova)


## Stage III
tp.stage3 <- as.matrix(tp.stage3)
tp.stage3.trim <- tp.stage3[,which(colSums(tp.stage3)>0)]
tp.stage3.melt <- melt(tp.stage3.trim)
colnames(tp.stage3.melt) <- c("PatientId", "Signature", "weight")
tp.stage3.melt$weight <- sapply(tp.stage3.melt$weight, function(x) 100*x)

pdf("./Figures/TCGA_SKCM_TP_mutsigs_stage3_boxplot.pdf", w=10, h=6)
ggplot(tp.stage3.melt, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM TP Stage III") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

#Anova
tp.stage3.anova <- aov(weight~Signature,  data=tp.stage3.melt)
summary(tp.stage3.anova)


## Stage IV
tp.stage4 <- as.matrix(tp.stage4)
tp.stage4.trim <- tp.stage4[,which(colSums(tp.stage4)>0)]
tp.stage4.melt <- melt(tp.stage4.trim)
colnames(tp.stage4.melt) <- c("PatientId", "Signature", "weight")
tp.stage4.melt$weight <- sapply(tp.stage4.melt$weight, function(x) 100*x)

#Boxplot
pdf("./Figures/TCGA_SKCM_TP_mutsigs_stage4_boxplot.pdf", w=10, h=6)
ggplot(tp.stage4.melt, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis SKCM TP Stage IV") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

#Anova
tp.stage4.anova <- aov(weight~Signature,  data=tp.stage4.melt)
summary(tp.stage4.anova)



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
table(tp.clin$prior_systemic_therapy_type, useNA="always")
table(tp.clin$radiation_therapy, useNA="always")
table(tp.clin$history_of_neoadjuvant_treatment, useNA="always")

table(tm.clin$prior_systemic_therapy_type, useNA="always")
table(tm.clin$radiation_therapy, useNA="always")
table(tm.clin$history_of_neoadjuvant_treatment, useNA="always")

# Therapy subgroups



### TO ASK
# Therapy annotations:  interferon 90 day prior excision admin indicator -- immunotherapy before surgery
#                       history of neoadjuvant treatment (either no or NA)
#                               = chemo after surgery




# cancer stage barplot: use avg --> complementary to boxplot sigs analysis

