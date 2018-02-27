# Mutational Analysis - SKCM Cancer
# DT Laura Vo Ngoc
# Start: 18/02/2018

library("BSgenome.Hsapiens.UCSC.hg38")
library(deconstructSigs)
library(ggplot2)
library(reshape2)
library(TCGAbiolinks)


#### WD & LOAD FILES ####
setwd("C:/Users/rockp/Desktop/UCL/Project/Project-Laura")

load("./Data/TCGAmutations_SKCM.RData")             # Raw TCGA mutation data
load("./Output/TCGAmutations_SKCM_snvs.RData")      # Selected SNVs
# Primary tumor
load("./Output/mutsig_input_SKCM_PT.RData")            # Mutation analysis input
load("./Output/TCGAmutsigs_cut0.00_SKCM.RData")     # Mutational signatures output
# Metastatic tumor
load("./Output/mutsig_input_SKCM_TM.RData")            # Mutation analysis input
load("./Output/TCGAmutsigs_cut0.00_SKCM_TM.RData")     # Mutational signatures output

#### DATA OF INTEREST ####
# Download TCGA mutations for SKCM cancer:
mutations <- GDCquery_Maf("SKCM", pipelines = "mutect2")
save(mutations, file="./Data/TCGAmutations_SKCM.RData")

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
save(snvs, file="./Output/TCGAmutations_SKCM_snvs.RData")

# Stores each type into separate variables
snvs_pt <- snvs[which(snvs$Sample_Type == "01"),]           # 61990 samples
snvs_met <- snvs[which(snvs$Sample_Type == "06"),]          # 326634 samples


#### MUTATIONAL SIGNATURES ANALYSIS -- PT ####
# Convert to deconstructSigs input:
sigs_input_pt <- mut.to.sigs.input(mut.ref = snvs_pt, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

save(sigs_input_pt, file="./Output/mutsig_input_SKCM_PT.RData")


# Signatures for all samples
#mut_sigs_pt <- as.data.frame(t(sapply(1:nrow(sigs_input_pt), 
 #                                  function(x) whichSignatures(tumor.ref = sigs_input_pt,
  #                                       signatures.ref = signatures.cosmic,
   #                                      sample.id = rownames(sigs_input_pt[x,]), 
    #                                     contexts.needed = TRUE,
     #                                    signature.cutoff = 0.00,              # default = 0.06
      #                                   tri.counts.method = 'exome'))), 
       #                   row.names=rownames(sigs_input_pt))

weights_pt <- as.data.frame(t(sapply(1:nrow(sigs_input_pt), 
                                    function(x) whichSignatures(tumor.ref = sigs_input_pt,
                                                                signatures.ref = signatures.cosmic,
                                                                sample.id = rownames(sigs_input_pt[x,]), 
                                                                contexts.needed = TRUE,
                                                                signature.cutoff = 0.00,              # default = 0.06
                                                                tri.counts.method = 'exome')$weights)), 
                           row.names=rownames(sigs_input_pt))

weights <- sapply(weights_pt, function(x) unlist(x))
weights <- weights[,-31]

rownames(weights) <- rownames(weights_pt)

save(mut_sigs_pt, file="./Output/TCGAmutsigs_cut0.00_SKCM_PT.RData")
save(weights_pt, file="./Output/TCGAweights_cut0.00_SKCM_PT.RData")


# Reformatting the weights data 
melted <- melt(weights)
colnames(melted) <- c("PatientId", "Signature", "weight")
melted$weight <- sapply(melted$weight, function(x) 100*x)


# Plot: stacked barplot signatures
pdf("./Figures/TCGAmutsig_SCKM_PT.pdf", w=10, h=6)
ggplot(melted, aes(x=PatientId, y=weight, fill=Signature)) + 
    geom_bar(stat="identity") +
    ggtitle("Mutational signature analysis SCKM primary tumor ") +
    xlab("Sample") +
    ylab("Contribution (%)") +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) #+
    #geom_text(label=paste0(melted$weight,"%"))
dev.off()
    





#### MUTATIONAL SIGNATURES ANALYSIS -- TM ####
# Convert to deconstructSigs input:
sigs_input_tm <- mut.to.sigs.input(mut.ref = snvs_met, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

save(sigs_input_tm, file="./Output/mutsig_input_SKCM_TM.RData")


# Signatures for all samples

## sapply(rownames(), function(x) sample.id = x)

mut_sigs_tm <- as.data.frame(t(sapply(1:nrow(sigs_input_tm), 
                                   function(x) whichSignatures(tumor.ref = sigs_input_tm,
                                                               signatures.ref = signatures.cosmic,
                                                               sample.id = rownames(sigs_input_tm[x,]), 
                                                               contexts.needed = TRUE,
                                                               signature.cutoff = 0.00,              # default = 0.06
                                                               tri.counts.method = 'exome'))), 
                          row.names=rownames(sigs_input_tm))

save(mut_sigs_tm, file="./Output/TCGAmutsigs_cut0.00_SKCM_TM.RData")



# Pie chart test (single sample)
pdf("./Figure/mutsig_SKCM_sample3_tm_pie.pdf")
makePie(mut_sigs_tm[3,], sub="TCGA SKCM TM")
dev.off()



