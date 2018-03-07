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
load("./Output/TCGA_UVM_PT_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_UVM_PT_weights_cut0.00.RData")     # Mutational signatures output


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
snvs_pt <- snvs[which(snvs$Sample_Type == "01"),]
snvs_met <- snvs[which(snvs$Sample_Type == "06"),]


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

save(sigs_input_pt, file="./Output/TCGA_UVM_PT_mutsig_input.RData")


# Signatures
weights_pt <- as.data.frame(t(sapply(rownames(sigs_input_pt), 
                                     function(x) whichSignatures(tumor.ref = sigs_input_pt,
                                                                 signatures.ref = signatures.cosmic,
                                                                 sample.id = x, 
                                                                 contexts.needed = TRUE,
                                                                 signature.cutoff = 0.06,              # default = 0.06
                                                                 tri.counts.method = 'exome')$weights)), 
                            row.names=rownames(sigs_input_pt))

weight_pt <- as.data.frame(sapply(weights_pt, function(x) unlist(x)), row.names = rownames(weights_pt))

# Change Signature.# to S#
colnames(weight_pt) <- sapply(colnames(weight_pt), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))

# Add unknown contribution column
weight_pt$unknown <- (1-rowSums(weight_pt))

save(weight_pt, file="./Output/TCGA_UVM_PT_weights_cut0.00.RData")



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
colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#BFA89E", "#826C7F", "#FF7F00", 
             "#CAB2D6", "#6A3D9A", "#D6DBD2", "#FFFF99", "#B15928", "#8DD3C7", "#FFFFB3", "#E75A7C", "#BEBADA", 
             "#FB8072", "#40476D", "#80B1D3", "#FDB462", "#B3DE69", "#258EA6", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", 
             "#000000")
pdf("./Figures/TCGA_UVM_PT_mutsig.pdf", w=10, h=6)
ggplot(melt_pt, aes(x=PatientId, y=weight, fill=Signature)) + 
    scale_fill_manual(values=colors) +
    geom_bar(stat="identity") +
    ggtitle("Mutational signature analysis UVM primary tumor ") +
    xlab("Sample") +
    ylab("Contribution (%)") +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
dev.off()

# Plot: boxplot signatures
pdf("./Figures/TCGA_UVM_PT_mutsigs_boxplot.pdf", w=10, h=6)
ggplot(melt_pt, aes(x=Signature, y=weight)) +
    geom_boxplot() +
    ggtitle("Mutational signature analysis UVM primary tumor") +
    xlab("Signatures") +
    ylab("Contribution")
dev.off()

# Plot: heatmap
melted_pt$weight <- as.numeric(melted_pt$weight)
pt_mat <- acast(PatientId~Signature,data=melted_pt, value="weight", fun.aggregate=mean)
pdf("./Figures/TCGA_UVM_PT_mutsigs_heatmap.pdf", w=10, h=6, onefile=FALSE)
p_pt <- pheatmap(t(pt_mat), show_colnames = FALSE, main="Mutational signatures UVM primary tumor")
dev.off()

#trim_pt <- as.data.frame(trim_pt)


#### Mean contributions per signature ####
means <- data.frame(Signatures=c(colnames(trim_pt)), value=c(colMeans(trim_pt)*100))

pdf("./Figures/TCGA_UVM_PT_mutsigs_avg.pdf", w=12, h=6)
ggplot(means, aes(x=reorder(Signatures, -value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    ggtitle("Mean mutational signature contributions UVM primary tumor") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5))
dev.off()

save(means, file="./Output/TCGA_UVM_PT_mutsigs_means.RData")