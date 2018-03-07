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

#### WD & LOAD FILES ####
setwd("C:/Users/rockp/Desktop/UCL/Project/Project-Laura")

load("./Data/TCGA_SKCM_mutations.RData")             # Raw TCGA mutation data
load("./Output/TCGA_SKCM_snvs.RData")      # Selected SNVs

# Primary tumor
load("./Output/TCGA_SKCM_PT_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_SKCM_PT_weights_cut0.00.RData")     # Mutational signatures output

# Metastatic tumor
load("./Output/TCGA_SKCM_TM_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_SKCM_TM_weights_cut0.00.RData")     # Mutational signatures output


#**********************************************************************************************************************
#### DATA OF INTEREST ####

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
snvs$Sample_Type <- sapply(snvs$Sample_Type, function(x) strsplit(x, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)[[1]][1])

# Save
save(snvs, file="./Output/TCGA_SKCM_snvs.RData")

# Stores each type into separate variables
snvs_pt <- snvs[which(snvs$Sample_Type == "01"),]           # 61990 samples
snvs_met <- snvs[which(snvs$Sample_Type == "06"),]          # 326634 samples


#**********************************************************************************************************************
#### MUTATIONAL SIGNATURES ANALYSIS -- PT ####

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


#### S1 and S7 analysis -- PT ####
sp1_7 <- weight_pt[,c(1, 7)]
sp1_7 <- as.matrix(sp1_7)
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


#### Mean contributions per signature ####

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

#**********************************************************************************************************************
#### MUTATIONAL SIGNATURES ANALYSIS -- TM ####

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


#### S1 and S7 analysis -- TM ####
sm1_7 <- weight_tm[,c(1, 7)]
sm1_7 <- as.matrix(sm1_7)
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


#### Mean contributions per signature ####

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

#**********************************************************************************************************************
#### PRIMARY VS. METASTATIC ####

pt <- weight_pt[,which(colSums(weight_pt)>0)]
tm <- weight_tm[,which(colSums(weight_tm)>0)]
t.test(pt$S1, pt$S7)
t.test(tm$S1, tm$S7)
t.test(pt$S1, tm$S1)
t.test(pt$S7, tm$S7)



pt$participant <- sapply(rownames(pt), function(x) paste0(strsplit(x,"-")[[1]][2], "-", strsplit(x,"-")[[1]][3]))
tm$participant <- sapply(rownames(tm), function(x) paste0(strsplit(x,"-")[[1]][2], "-", strsplit(x,"-")[[1]][3]))
matched <- pt[intersect(pt$participant, tm$participant),]       #none

save(pt, file="./Output/TCGA_SKCM_PT_weights_trim_cut0.00.RData")
save(tm, file="./Output/TCGA_SKCM_TM_weights_trim_cut0.00.RData")




#**********************************************************************************************************************
#### SKCM_PT, SKCM_TM, UVM ####

load("./Output/TCGA_SKCM_PT_mutsigs_means.RData")   # means_pt
load("./Output/TCGA_SKCM_TM_mutsigs_means.RData")   # means_tm
load("./Output/TCGA_UVM_PT_mutsigs_means.RData")    # means

colnames(means_pt)[2] <- "SKCM_PT"
colnames(means_tm)[2] <- "SKCM_TM"
colnames(means)[2] <- "UVM_PT"

merged <- merge(means_pt, means_tm, by="Signatures", all=TRUE)
merged <- merge(merged, means, by="Signatures", all=TRUE)
merged[is.na(merged)] <- 0
merged <- merged[order(-weight_pt$S1),]

melt_merged <- melt(merged)
colnames(melt_merged) <- c("Signatures", "Sample", "value")

factors = c(colnames(weight_pt))

colors <- c(brewer.pal(12, "Paired"), brewer.pal(11, "Set3"),"#000000")

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