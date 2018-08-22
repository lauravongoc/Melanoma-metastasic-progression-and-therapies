# Excision repair
# DT Laura Vo Ngoc
# Start: 05/08/2018

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

# Clinical and mutsigs data
load("./Output/TCGA_SKCM_TP_clinical.RData")    # SKCM TP
load("./Output/TCGA_SKCM_TM_clinical.RData")    # SKCM TM
load("./Output/TCGA_UVM_TP_clinical.RData")     # UVM TP 

# Stage
load("./Output/TCGA_SKCM_TP_stage.RData")       # SKCM TP
load("./Output/TCGA_SKCM_TM_stage.RData")       # SKCM TM

# Mutsigs, metscores, and stage
load("./Output/TCGA_SKCM_TP_metastatic_score_stage.RData")      # SKCM TP
load("./Output/TCGA_SKCM_TM_metastatic_score_stage.RData")      # SKCM TM
load("./Output/TCGA_UVM_TP_metastatic_score_stage.RData")       # UVM TP

# Survival
load("./Output/TCGA_SKCM_TP_survival.RData")        # SKCM TP
load("./Output/TCGA_SKCM_TM_survival.RData")        # SKCM TM
load("./Output/TCGA_UVM_TP_survival.RData")         # UVM

# Genelists
load("./Output/NER_genelist.RData")     # NER genes
load("./Output/MMR_genelist.RData")     # MMR genes

# Expression data NER/MMR
load("./Data/TCGA_SKCM_TP_NER_expression.RData")    # SKCM TP NER
load("./Data/TCGA_SKCM_TM_NER_expression.RData")    # SKCM TM NER
load("./Data/TCGA_UVM_NER_expression.RData")        # UVM NER

load("./Data/TCGA_SKCM_TP_MMR_expression.RData")    # SKCM TP MMR
load("./Data/TCGA_SKCM_TM_MMR_expression.RData")    # SKCM TM MMR
load("./Data/TCGA_UVM_MMR_expression.RData")        # UVM MMR


# Total mutsigs
load("./Output/TCGA_SKCM_TP_total_mutsigs.RData")       # SKCM TP
load("./Output/TCGA_SKCM_TM_total_mutsigs.RData")       # SKCM TM
load("./Output/TCGA_UVM_TP_total_mutsigs.RData")        # UVM TP

# Mutation files
load("./Data/TCGA_UVM_ALLmut.RData")    # UVM

#--------- COLORS -----------------------------------------------------------------------------------------------------
colors <- c(brewer.pal(12, "Paired"), brewer.pal(11, "Set3"),"#000000")
colors2 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#BFA89E", "#826C7F", "#FF7F00",
             "#CAB2D6", "#8B786D", "#6A3D9A", "#D6DBD2", "#FFFF99", "#B15928", "#8DD3C7", "#FFFFB3", "#E75A7C", "#BEBADA",
             "#FB8072", "#40476D", "#80B1D3", "#FDB462", "#B3DE69", "#258EA6", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", 
             "#000000")
colors3 <- brewer.pal(3, "Paired")
colors4 <- c(brewer.pal(4, "Paired"),"#000000")
colors6 <- c("#fb9a99", brewer.pal(4, "Paired"),"#000000")
colors7 <- c("#1F78B4", "#B2DF8A", "#33A02C", "#000000")


#--------- IMPORT GENES -------------------------------------------------------------------------------------------------
# Read table
ER.genes <- read.table(file="./ExcisionRepairGenes.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)

# Select only genes in NER or MMR pathways
ER.genes <- ER.genes[which(ER.genes$Pathway.1=="NER" | ER.genes$Pathway.1=="MMR"),]
#ER.genes1 <- c("XPA", "RPA1", "RPA2", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5")      # 8 genes

ner.genes <- ER.genes$Gene.ID[which(ER.genes$Pathway.1=="NER")]     # 70 genes
mmr.genes <- ER.genes$Gene.ID[which(ER.genes$Pathway.1=="MMR")]     # 27 genes

save(ner.genes, file="./Output/NER_genelist.RData")
save(mmr.genes, file="./Output/MMR_genelist.RData")

# Split NER genelist into 2
ner.genes1 <- ner.genes[1:35]
ner.genes2 <- ner.genes[36:70]


#--------- EXPRESSION DATA ----------------------------------------------------------------------------------------------
# SKCM TP
tp.skcm.expr.ner1 <- Samples.mRNASeq(format = "csv",
                                    gene = ner.genes1,
                                    cohort = "SKCM",
                                    sample_type = "TP",
                                    protocol = "RSEM",
                                    page_size=6000
)
tp.skcm.expr.ner2 <- Samples.mRNASeq(format = "csv",
                                    gene = ner.genes2,
                                    cohort = "SKCM",
                                    sample_type = "TP",
                                    protocol = "RSEM",
                                    page_size=6000
)
tp.skcm.expr.mmr <- Samples.mRNASeq(format = "csv",
                                    gene = mmr.genes,
                                    cohort = "SKCM",
                                    sample_type = "TP",
                                    protocol = "RSEM",
                                    page_size=6000
)
tp.skcm.expr.ner <- rbind(tp.skcm.expr.ner1, tp.skcm.expr.ner2)
save(tp.skcm.expr.ner, file="./Data/TCGA_SKCM_TP_NER_expression.RData")
save(tp.skcm.expr.mmr, file="./Data/TCGA_SKCM_TP_MMR_expression.RData")


# SKCM TM
tm.skcm.expr.ner1 <- Samples.mRNASeq(format = "csv",
                                    gene = ner.genes1,
                                    cohort = "SKCM",
                                    sample_type = "TM",
                                    protocol = "RSEM",
                                    page_size=6000
)
tm.skcm.expr.ner2 <- Samples.mRNASeq(format = "csv",
                                    gene = ner.genes2,
                                    cohort = "SKCM",
                                    sample_type = "TM",
                                    protocol = "RSEM",
                                    page_size=6000
)
tm.skcm.expr.mmr <- Samples.mRNASeq(format = "csv",
                                    gene = mmr.genes,
                                    cohort = "SKCM",
                                    sample_type = "TM",
                                    protocol = "RSEM",
                                    page_size=6000
)
tm.skcm.expr.ner <- rbind(tm.skcm.expr.ner1, tm.skcm.expr.ner2)
save(tm.skcm.expr.ner, file="./Data/TCGA_SKCM_TM_NER_expression.RData")
save(tm.skcm.expr.mmr, file="./Data/TCGA_SKCM_TM_MMR_expression.RData")

# UVM TP
uvm.expr.ner1 <- Samples.mRNASeq(format = "csv",
                                     gene = ner.genes1,
                                     cohort = "UVM",
                                     sample_type = "TP",
                                     protocol = "RSEM",
                                     page_size=6000
)
uvm.expr.ner2 <- Samples.mRNASeq(format = "csv",
                                     gene = ner.genes2,
                                     cohort = "UVM",
                                     sample_type = "TP",
                                     protocol = "RSEM",
                                     page_size=6000
)
uvm.expr.mmr <- Samples.mRNASeq(format = "csv",
                                    gene = mmr.genes,
                                    cohort = "UVM",
                                    sample_type = "TP",
                                    protocol = "RSEM",
                                    page_size=6000
)
uvm.expr.ner <- rbind(uvm.expr.ner1, uvm.expr.ner2)
save(uvm.expr.ner, file="./Data/TCGA_UVM_NER_expression.RData")
save(uvm.expr.mmr, file="./Data/TCGA_UVM_MMR_expression.RData")


#--------- DISTRIBUTION ----------------------------------------------------------------------------------------------
tp.skcm.expr.ner$expression_log2 <- as.numeric(tp.skcm.expr.ner$expression_log2)
tp.skcm.expr.mmr$expression_log2 <- as.numeric(tp.skcm.expr.mmr$expression_log2)
tm.skcm.expr.ner$expression_log2 <- as.numeric(tm.skcm.expr.ner$expression_log2)
tm.skcm.expr.mmr$expression_log2 <- as.numeric(tm.skcm.expr.mmr$expression_log2)
uvm.expr.ner$expression_log2 <- as.numeric(uvm.expr.ner$expression_log2)
uvm.expr.mmr$expression_log2 <- as.numeric(uvm.expr.mmr$expression_log2)


pdf("./Figures/TCGA_SKCM_UVM_metscore_hist_small_cutoff.pdf", w=8, h=6)
hist(tm.skcm.expr.ner$expression_log2, 
     col="#1F78B4B3",  
     main="NER genes expression distribution", 
     xlab="NER genes expression",
     xlim=c(-2.2,15.3),
     ylim=c(0,7000))
hist(tp.skcm.expr.ner$expression_log2, 
     col="#fb9a99B3", 
     main="Metastatic score distribution SKCM TP", 
     xlab="Metastatic score",
     xlim=c(-2.2,15.3),
     ylim=c(0,7000),
     add=TRUE)
hist(uvm.expr.ner$expression_log2,
     col="#B2DF8AB3", 
     main="Metastatic score distribution UVM TP", 
     xlab="Metastatic score",
     #breaks=14,
     xlim=c(-2.24,15.3),
     ylim=c(0,7000),
     add=TRUE)
legend("topright", c("SKCM TP", "SKCM TM", "UVM TP"), col=c("#fb9a99B3","#1F78B4B3", "#B2DF8AB3"), lwd=10, bty="n")
abline(v=1, col="black", lty=2, lwd=2)
dev.off()




#--------- METSCORE CORRELATION ----------------------------------------------------------------------------------------------

# Expression and metscore variables
tp.skcm.ner.metscore <- merge(tp.skcm.expr.ner, tp.skcm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)
tp.skcm.mmr.metscore <- merge(tp.skcm.expr.mmr, tp.skcm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)

tm.skcm.ner.metscore <- merge(tm.skcm.expr.ner, tm.skcm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)
tm.skcm.mmr.metscore <- merge(tm.skcm.expr.mmr, tm.skcm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)

uvm.ner.metscore <- merge(uvm.expr.ner, uvm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)
uvm.mmr.metscore <- merge(uvm.expr.mmr, uvm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)

# Loop running correlation for each gene

# NER genes
genes.ner <- unique(tp.skcm.expr.ner$gene)
ner.cor <- data.frame(gene=NA, cor_skcm_tp=NA, pval=NA, padj=NA, sign=NA, 
                      cor_skcm_tm=NA, pval_tm=NA, padj_tm=NA,sign_tm=NA, 
                      cor_uvm=NA, pval_uvm=NA, padj_uvm=NA, sign_uvm=NA)

for (i in 1:length(genes.ner)) {
    gen <- genes.ner[i]
    data <- tp.skcm.ner.metscore[which(tp.skcm.ner.metscore$gene==gen),]
    c <- cor.test(data$expression_log2, data$metscore)
    data1 <- tm.skcm.ner.metscore[which(tm.skcm.ner.metscore$gene==gen),]
    c1 <- cor.test(data1$expression_log2, data1$metscore)
    data2 <- uvm.ner.metscore[which(uvm.ner.metscore$gene==gen),]
    c2 <- cor.test(data2$expression_log2, data2$metscore)
    
    ner.cor[i,1] <- gen
    ner.cor[i,2] <- c$estimate[[1]]
    ner.cor[i,3] <- c$p.value[[1]]
    
    ner.cor[i,6] <- c1$estimate[[1]]
    ner.cor[i,7] <- c1$p.value[[1]]
    
    ner.cor[i,10] <- c2$estimate[[1]]
    ner.cor[i,11] <- c2$p.value[[1]]
}

ner.cor$padj <- p.adjust(ner.cor$pval, method="BH")
ner.cor$sign <- ifelse(ner.cor$padj < 0.05, "*","")
ner.cor$padj_tm <- p.adjust(ner.cor$pval_tm, method="BH")
ner.cor$sign_tm <- ifelse(ner.cor$padj_tm < 0.05, "*","")
ner.cor$padj_uvm <- p.adjust(ner.cor$pval_uvm, method="BH")
ner.cor$sign_uvm <- ifelse(ner.cor$padj_uvm < 0.05, "*","")

write.csv(ner.cor, file="./Output/TCGA_SKCM_UVM_NER_metscore_cor.csv")

# MMR genes
genes.mmr <- unique(tp.skcm.expr.mmr$gene)
mmr.cor <- data.frame(gene=NA, cor_skcm_tp=NA, pval=NA, padj=NA, sign=NA, 
                      cor_skcm_tm=NA, pval_tm=NA, padj_tm=NA, sign_tm=NA, 
                      cor_uvm=NA, pval_uvm=NA, padj_uvm=NA, sign_uvm=NA)

for (i in 1:length(genes.mmr)) {
    gen <- genes.mmr[i]
    data <- tp.skcm.mmr.metscore[which(tp.skcm.mmr.metscore$gene==gen),]
    c <- cor.test(data$expression_log2, data$metscore)
    data1 <- tm.skcm.mmr.metscore[which(tm.skcm.mmr.metscore$gene==gen),]
    c1 <- cor.test(data1$expression_log2, data1$metscore)
    data2 <- uvm.mmr.metscore[which(uvm.mmr.metscore$gene==gen),]
    c2 <- cor.test(data2$expression_log2, data2$metscore)
    
    mmr.cor[i,1] <- gen
    mmr.cor[i,2] <- c$estimate[[1]]
    mmr.cor[i,3] <- c$p.value[[1]]
    
    mmr.cor[i,6] <- c1$estimate[[1]]
    mmr.cor[i,7] <- c1$p.value[[1]]

    mmr.cor[i,10] <- c2$estimate[[1]]
    mmr.cor[i,11] <- c2$p.value[[1]]
}
mmr.cor$padj <- p.adjust(mmr.cor$pval, method="BH")
mmr.cor$sign <- ifelse(mmr.cor$padj < 0.05, "*","")
mmr.cor$padj_tm <- p.adjust(mmr.cor$pval_tm, method="BH")
mmr.cor$sign_tm <- ifelse(mmr.cor$padj_tm < 0.05, "*","")
mmr.cor$padj_uvm <- p.adjust(mmr.cor$pval_uvm, method="BH")
mmr.cor$sign_uvm <- ifelse(mmr.cor$padj_uvm < 0.05, "*","")

write.csv(mmr.cor, file="./Output/TCGA_SKCM_UVM_MMR_metscore_cor.csv")

#--------- TOTAL MUTSIGS CORRELATION ----------------------------------------------------------------------------------------------

# Expression and metscore variables
tp.skcm.ner.mutsig <- merge(tp.skcm.expr.ner, tp.skcm.totalmuts, by="tcga_participant_barcode", all.x=TRUE)
tp.skcm.mmr.mutsig <- merge(tp.skcm.expr.mmr, tp.skcm.totalmuts, by="tcga_participant_barcode", all.x=TRUE)

tm.skcm.ner.mutsig <- merge(tm.skcm.expr.ner, tm.skcm.totalmuts, by="tcga_participant_barcode", all.x=TRUE)
tm.skcm.mmr.mutsig <- merge(tm.skcm.expr.mmr, tm.skcm.totalmuts, by="tcga_participant_barcode", all.x=TRUE)

uvm.ner.mutsig <- merge(uvm.expr.ner, uvm.totalmuts, by="tcga_participant_barcode", all.x=TRUE)
uvm.mmr.mutsig <- merge(uvm.expr.mmr, uvm.totalmuts, by="tcga_participant_barcode", all.x=TRUE)


# ---- ~ Loops running correlation for each gene ----
# NER genes
genes.ner <- unique(tp.skcm.expr.ner$gene)
tp.ner.cor.sigs <- data.frame(gene=genes.ner)

for (j in c(12, 15, 19, 29, 30, 32, 37)) {
    sig <- colnames(tp.skcm.ner.mutsig)[j]
    
    df <- data.frame(gene=NA, cor=NA, pval=NA, padj=NA, sign=NA)
    colnames(df) <- c("gene", paste0(sig, "_cor"), paste0(sig, "_p"), paste0(sig, "_padj"), paste0(sig,"_sign"))
    
    for (i in 1:length(genes.ner)) {
        gen <- genes.ner[i]
        
        df[i,1] <- gen
        
        data <- tp.skcm.ner.mutsig[which(tp.skcm.ner.mutsig$gene==gen),]
        
        c <- cor.test(data$expression_log2, data[,j])
        
        df[i,2] <- c$estimate[[1]]
        df[i,3] <- c$p.value[[1]]
    }
    df[,4] <- p.adjust(df[,3], method="BH")
    df[,5] <- ifelse(df[,4]<0.05, "*", "")
    
    tp.ner.cor.sigs <- merge(tp.ner.cor.sigs, df, by="gene")
}

write.csv(tp.ner.cor.sigs, file="./Output/TCGA_SKCM_TP_NER_totmutsigs_cor.csv")


tm.ner.cor.sigs <- data.frame(gene=genes.ner)

for (j in c(12, 15, 19, 29, 30, 32, 37)) {
    sig <- colnames(tm.skcm.ner.mutsig)[j]
    
    df <- data.frame(gene=NA, cor=NA, pval=NA, padj=NA, sign=NA)
    colnames(df) <- c("gene", paste0(sig, "_cor"), paste0(sig, "_p"), paste0(sig, "_padj"), paste0(sig,"_sign"))
    
    for (i in 1:length(genes.ner)) {
        gen <- genes.ner[i]
        
        df[i,1] <- gen
        
        data <- tm.skcm.ner.mutsig[which(tm.skcm.ner.mutsig$gene==gen),]
        
        c <- cor.test(data$expression_log2, data[,j])
        
        df[i,2] <- c$estimate[[1]]
        df[i,3] <- c$p.value[[1]]
    }
    df[,4] <- p.adjust(df[,3], method="BH")
    df[,5] <- ifelse(df[,4]<0.05, "*", "")
    
    tm.ner.cor.sigs <- merge(tm.ner.cor.sigs, df, by="gene")
}

write.csv(tm.ner.cor.sigs, file="./Output/TCGA_SKCM_TM_NER_totmutsigs_cor.csv")

uvm.ner.cor.sigs <- data.frame(gene=genes.ner)

for (j in c(12, 15, 19, 29, 30, 32, 37)) {
    sig <- colnames(uvm.ner.mutsig)[j]
    
    df <- data.frame(gene=NA, cor=NA, pval=NA, padj=NA, sign=NA)
    colnames(df) <- c("gene", paste0(sig, "_cor"), paste0(sig, "_p"), paste0(sig, "_padj"), paste0(sig,"_sign"))
    
    for (i in 1:length(genes.ner)) {
        gen <- genes.ner[i]
        
        df[i,1] <- gen
        
        data <- uvm.ner.mutsig[which(uvm.ner.mutsig$gene==gen),]
        
        c <- cor.test(data$expression_log2, data[,j])
        
        df[i,2] <- c$estimate[[1]]
        df[i,3] <- c$p.value[[1]]
    }
    df[,4] <- p.adjust(df[,3], method="BH")
    df[,5] <- ifelse(df[,4]<0.05, "*", "")
    
    uvm.ner.cor.sigs <- merge(uvm.ner.cor.sigs, df, by="gene")
}

write.csv(uvm.ner.cor.sigs, file="./Output/TCGA_UVM_TP_NER_topmutsigs_cor.csv")

# MMR genes
genes.mmr <- unique(tp.skcm.expr.mmr$gene)
tp.mmr.cor.sigs <- data.frame(gene=genes.mmr)

for (j in c(14, 23, 28, 33)) {
    sig <- colnames(tp.skcm.mmr.mutsig)[j]
    
    df <- data.frame(gene=NA, cor=NA, pval=NA, padj=NA, sign=NA)
    colnames(df) <- c("gene", paste0(sig, "_cor"), paste0(sig, "_p"), paste0(sig, "_padj"), paste0(sig,"_sign"))
    
    for (i in 1:length(genes.mmr)) {
        gen <- genes.mmr[i]
        
        df[i,1] <- gen
        
        data <- tp.skcm.mmr.mutsig[which(tp.skcm.mmr.mutsig$gene==gen),]
        
        c <- cor.test(data$expression_log2, data[,j])
        
        df[i,2] <- c$estimate[[1]]
        df[i,3] <- c$p.value[[1]]
    }
    df[,4] <- p.adjust(df[,3], method="BH")
    df[,5] <- ifelse(df[,4]<0.05, "*", "")
    
    tp.mmr.cor.sigs <- merge(tp.mmr.cor.sigs, df, by="gene")
}

write.csv(tp.mmr.cor.sigs, file="./Output/TCGA_SKCM_TP_MMR_totmutsigs_cor.csv")

tm.mmr.cor.sigs <- data.frame(gene=genes.mmr)

for (j in c(14, 23, 28, 33)) {
    sig <- colnames(tm.skcm.mmr.mutsig)[j]
    
    df <- data.frame(gene=NA, cor=NA, pval=NA, padj=NA, sign=NA)
    colnames(df) <- c("gene", paste0(sig, "_cor"), paste0(sig, "_p"), paste0(sig, "_padj"), paste0(sig,"_sign"))
    
    for (i in 1:length(genes.mmr)) {
        gen <- genes.mmr[i]
        
        df[i,1] <- gen
        
        data <- tm.skcm.mmr.mutsig[which(tm.skcm.mmr.mutsig$gene==gen),]
        
        c <- cor.test(data$expression_log2, data[,j])
        
        df[i,2] <- c$estimate[[1]]
        df[i,3] <- c$p.value[[1]]
    }
    df[,4] <- p.adjust(df[,3], method="BH")
    df[,5] <- ifelse(df[,4]<0.05, "*", "")
    
    tm.mmr.cor.sigs <- merge(tm.mmr.cor.sigs, df, by="gene")
}

write.csv(tm.mmr.cor.sigs, file="./Output/TCGA_SKCM_TM_MMR_totmutsigs_cor.csv")

uvm.mmr.cor.sigs <- data.frame(gene=genes.mmr)

for (j in c(14, 23, 28, 33)) {
    sig <- colnames(uvm.mmr.mutsig)[j]
    
    df <- data.frame(gene=NA, cor=NA, pval=NA, padj=NA, sign=NA)
    colnames(df) <- c("gene", paste0(sig, "_cor"), paste0(sig, "_p"), paste0(sig, "_padj"), paste0(sig,"_sign"))
    
    for (i in 1:length(genes.mmr)) {
        gen <- genes.mmr[i]
        
        df[i,1] <- gen
        
        data <- uvm.mmr.mutsig[which(uvm.mmr.mutsig$gene==gen),]
        
        c <- cor.test(data$expression_log2, data[,j])
        
        df[i,2] <- c$estimate[[1]]
        df[i,3] <- c$p.value[[1]]
    }
    df[,4] <- p.adjust(df[,3], method="BH")
    df[,5] <- ifelse(df[,4]<0.05, "*", "")
    
    uvm.mmr.cor.sigs <- merge(uvm.mmr.cor.sigs, df, by="gene")
}

write.csv(uvm.mmr.cor.sigs, file="./Output/TCGA_UVM_TP_MMR_topmutsigs_cor.csv")


#--------- COMBINE MUTSIGS ----------------------------------------------------------------------------------------------

tp.skcm.ner.mutsig$NER <- tp.skcm.ner.mutsig$S4 + tp.skcm.ner.mutsig$S7 + tp.skcm.ner.mutsig$S11 + tp.skcm.ner.mutsig$S22 + tp.skcm.ner.mutsig$S24 + tp.skcm.ner.mutsig$S29
tp.skcm.mmr.mutsig$MMR <- tp.skcm.mmr.mutsig$S6 + tp.skcm.mmr.mutsig$S15 + tp.skcm.mmr.mutsig$S20 + tp.skcm.mmr.mutsig$S26

tm.skcm.ner.mutsig$NER <- tm.skcm.ner.mutsig$S4 + tm.skcm.ner.mutsig$S7 + tm.skcm.ner.mutsig$S11 + tm.skcm.ner.mutsig$S22 + tm.skcm.ner.mutsig$S24 + tm.skcm.ner.mutsig$S29
tm.skcm.mmr.mutsig$MMR <- tm.skcm.mmr.mutsig$S6 + tm.skcm.mmr.mutsig$S15 + tm.skcm.mmr.mutsig$S20 + tm.skcm.mmr.mutsig$S26

uvm.ner.mutsig$NER <- uvm.ner.mutsig$S4 + uvm.ner.mutsig$S7 + uvm.ner.mutsig$S11 + uvm.ner.mutsig$S22 + uvm.ner.mutsig$S24 + uvm.ner.mutsig$S29
uvm.mmr.mutsig$MMR <- uvm.mmr.mutsig$S6 + uvm.mmr.mutsig$S15 + uvm.mmr.mutsig$S20 + uvm.mmr.mutsig$S26


# ---- ~ NER/MMR genes vs mutsigs ----

# NER genes
genes.ner <- unique(tp.skcm.expr.ner$gene)
ner.cor.comb <- data.frame(gene=NA, cor_skcm_tp=NA, pval=NA, padj=NA, sign=NA, 
                           cor_skcm_tm=NA, pval_tm=NA, padj_tm=NA,sign_tm=NA, 
                           cor_uvm=NA, pval_uvm=NA, padj_uvm=NA, sign_uvm=NA)

for (i in 1:length(genes.ner)) {
    gen <- genes.ner[i]
    data <- tp.skcm.ner.mutsig[which(tp.skcm.ner.mutsig$gene==gen),]
    c <- cor.test(data$expression_log2, data$NER)
    data1 <- tm.skcm.ner.mutsig[which(tm.skcm.ner.mutsig$gene==gen),]
    c1 <- cor.test(data1$expression_log2, data1$NER)
    data2 <- uvm.ner.mutsig[which(uvm.ner.mutsig$gene==gen),]
    c2 <- cor.test(data2$expression_log2, data2$NER)
    
    ner.cor.comb[i,1] <- gen
    ner.cor.comb[i,2] <- c$estimate[[1]]
    ner.cor.comb[i,3] <- c$p.value[[1]]
    
    ner.cor.comb[i,6] <- c1$estimate[[1]]
    ner.cor.comb[i,7] <- c1$p.value[[1]]
    
    ner.cor.comb[i,10] <- c2$estimate[[1]]
    ner.cor.comb[i,11] <- c2$p.value[[1]]
}

ner.cor.comb$padj <- p.adjust(ner.cor.comb$pval, method="BH")
ner.cor.comb$sign <- ifelse(ner.cor.comb$padj < 0.05, "*","")
ner.cor.comb$padj_tm <- p.adjust(ner.cor.comb$pval_tm, method="BH")
ner.cor.comb$sign_tm <- ifelse(ner.cor.comb$padj_tm < 0.05, "*","")
ner.cor.comb$padj_uvm <- p.adjust(ner.cor.comb$pval_uvm, method="BH")
ner.cor.comb$sign_uvm <- ifelse(ner.cor.comb$padj_uvm < 0.05, "*","")

write.csv(ner.cor.comb, file="./Output/TCGA_SKCM_UVM_NER_mutsigs_combo.csv")

# MMR genes
genes.mmr <- unique(tp.skcm.expr.mmr$gene)
mmr.cor.comb <- data.frame(gene=NA, cor_skcm_tp=NA, pval=NA, padj=NA, sign=NA, 
                           cor_skcm_tm=NA, pval_tm=NA, padj_tm=NA,sign_tm=NA, 
                           cor_uvm=NA, pval_uvm=NA, padj_uvm=NA, sign_uvm=NA)

for (i in 1:length(genes.mmr)) {
    gen <- genes.mmr[i]
    data <- tp.skcm.mmr.mutsig[which(tp.skcm.mmr.mutsig$gene==gen),]
    c <- cor.test(data$expression_log2, data$MMR)
    data1 <- tm.skcm.mmr.mutsig[which(tm.skcm.mmr.mutsig$gene==gen),]
    c1 <- cor.test(data1$expression_log2, data1$MMR)
    data2 <- uvm.mmr.mutsig[which(uvm.mmr.mutsig$gene==gen),]
    c2 <- cor.test(data2$expression_log2, data2$MMR)
   
     mmr.cor.comb[i,1] <- gen
    mmr.cor.comb[i,2] <- c$estimate[[1]]
    mmr.cor.comb[i,3] <- c$p.value[[1]]
    
    mmr.cor.comb[i,6] <- c1$estimate[[1]]
    mmr.cor.comb[i,7] <- c1$p.value[[1]]
    
    mmr.cor.comb[i,10] <- c2$estimate[[1]]
    mmr.cor.comb[i,11] <- c2$p.value[[1]]
}

mmr.cor.comb$padj <- p.adjust(mmr.cor.comb$pval, method="BH")
mmr.cor.comb$sign <- ifelse(mmr.cor.comb$padj < 0.05, "*","")
mmr.cor.comb$padj_tm <- p.adjust(mmr.cor.comb$pval_tm, method="BH")
mmr.cor.comb$sign_tm <- ifelse(mmr.cor.comb$padj_tm < 0.05, "*","")
mmr.cor.comb$padj_uvm <- p.adjust(mmr.cor.comb$pval_uvm, method="BH")
mmr.cor.comb$sign_uvm <- ifelse(mmr.cor.comb$padj_uvm < 0.05, "*","")

write.csv(mmr.cor.comb, file="./Output/TCGA_SKCM_UVM_MMR_mutsigs_combo.csv")

# ---- ~ Metscore vs mutsigs ----

tp.skcm.ner <- merge(tp.skcm.ner.mutsig, tp.skcm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)
tp.skcm.mmr <- merge(tp.skcm.mmr.mutsig, tp.skcm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)

tm.skcm.ner <- merge(tm.skcm.ner.mutsig, tm.skcm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)
tm.skcm.mmr <- merge(tm.skcm.mmr.mutsig, tm.skcm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)

uvm.ner <- merge(uvm.ner.mutsig, uvm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)
uvm.mmr <- merge(uvm.mmr.mutsig, uvm.metscore[c(1,3)], by="tcga_participant_barcode", all.x=TRUE)

tp.skcm.ner$genelist <- "NER"
tp.skcm.mmr$genelist <- "MMR"

colnames(tp.skcm.ner)[41] <- "comb_sigs"
colnames(tp.skcm.mmr)[41] <- "comb_sigs"

tp.skcm.ner.mmr <- rbind(tp.skcm.ner[,c(1:3,9:43)], tp.skcm.mmr[,c(1:3,9:43)])

tm.skcm.ner$genelist <- "NER"
tm.skcm.mmr$genelist <- "MMR"

colnames(tm.skcm.ner)[41] <- "comb_sigs"
colnames(tm.skcm.mmr)[41] <- "comb_sigs"

tm.skcm.ner.mmr <- rbind(tm.skcm.ner[,c(1:3,9:43)], tm.skcm.mmr[,c(1:3,9:43)])

uvm.ner$genelist <- "NER"
uvm.mmr$genelist <- "MMR"

colnames(uvm.ner)[41] <- "comb_sigs"

uvm.ner.mmr <- rbind(uvm.ner[,c(1:3,9:43)], uvm.mmr[,c(1:3,9:43)])


# Correlations
cor.test(tp.skcm.ner$metscore, tp.skcm.ner$NER)
cor.test(tp.skcm.mmr$metscore, tp.skcm.mmr$MMR)
cor.test(tm.skcm.ner$metscore, tm.skcm.ner$NER)
cor.test(tm.skcm.mmr$metscore, tm.skcm.mmr$MMR)
cor.test(uvm.ner$metscore, uvm.ner$NER)
cor.test(uvm.mmr$metscore, uvm.mmr$MMR)

# Plots

# SKCM TP
quantile(wt.s7$S7)
test <- cor.test(wt.s7$S7, wt.s7$metscore)
pval <- round(test$p.value, digits=4)
cor <- round(test$estimate[[1]], digits=7)

 ggplot(tp.skcm.ner, aes(x=NER, y=metscore)) +
    geom_point(shape=16, size=5) +
    geom_smooth(method=lm) +
    #facet_wrap(~genelist) +
    ggtitle("BRAF WT") +
    xlab("Total mutations contribution") +
    ylab("Metastatic score") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,2000,250)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,1976), ylim=c(0.6,1.6), expand=FALSE) +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13))
    #annotate("text", x=1500, y=c(1.55, 1.48), label = c(paste0("Correlation: ", cor), paste0("p-value: ", pval)))

 
 #--------- CETN2 -- SKCM TP ----------------------------------------------------------------------------------------------
 
cetn2.sig <- tp.skcm.ner.mutsig[which(tp.skcm.ner.mutsig$gene=="CETN2"),]     # 103 samples
cetn2.metscore <- tp.skcm.ner.metscore[which(tp.skcm.ner.metscore$gene=="CETN2"),]     # 103 samples
cetn2 <- merge(cetn2.sig, cetn2.metscore[,c(1,9)], by="tcga_participant_barcode")

s11.q <- quantile(cetn2.sig$S11)

cetn2 <- cetn2[which(cetn2$S11>s11.q[2][[1]] & cetn2$S11 < s11.q[4][[1]]),]

quantile(cetn2$metscore)
test <- cor.test(cetn2$S7, cetn2$metscore)
pval <- round(test$p.value, digits=4)
cor <- round(test$estimate[[1]], digits=7)
 
pdf("./Figures/TCGA_SKCM_TP_BRAF_WT_S7_metscore_scatter.pdf", w=8, h=6)
ggplot(uvm.mmr.mutsig[which(uvm.mmr.mutsig$gene=="RFC3"),], aes(x=expression_log2, y=S6)) +
    geom_point(shape=16, size=5, color="#1F78B4B3") + 
    geom_smooth(color="black", method=lm) + 
    #ggtitle("BRAF WT") +
    xlab("log2(RFC3 expresssion)") +
    ylab("Total mutation contribution S6") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    #scale_x_continuous(breaks=seq(0,2000,250)) +
    #scale_y_continuous(breaks=seq(0,10,0.2)) +
    #coord_cartesian(xlim=c(0,1976), ylim=c(0.6,1.6), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
           axis.title = element_text(size=13))
    #annotate("text", x=1500, y=c(1.55, 1.48), label = c(paste0("Correlation: ", cor), paste0("p-value: ", pval)))
dev.off()
 

 #--------- SURVIVAL ----------------------------------------------------------------------------------------------
 tp.surv.repair <- merge(tp.surv[,c(1:105)], tp.skcm.ner.mmr[,c(1:3,36,38)], by="tcga_participant_barcode")
 tp.surv.ner <- tp.surv.repair[which(tp.surv.repair$genelist=="NER"),]
 tp.surv.mmr <- tp.surv.repair[which(tp.surv.repair$genelist=="MMR"),]
 
 summary(coxph(Surv(days_to_death, vital_status)~stage+S4+S6+S7+S11+S15+S20+S22+S23+S26,
               data=tp.surv.repair))
 
 summary(coxph(Surv(days_to_death, vital_status)~stage+S4+S6+S7+S11+S15+S20+S22+S23+S26,
               data=tp.surv.ner))
 
 # Define cutoff
 tp.surv.ner$comb_sigs_hilow <- tp.surv.ner$comb_sigs
 ner_q <-quantile(tp.surv.ner$comb_sigs)
 
 for (i in 1:nrow(tp.surv.ner)) {
     if (tp.surv.ner$comb_sigs[i] < ner_q[2][[1]]) {
         tp.surv.ner$comb_sigs_hilow[i] <- "low"
     }
     else if (tp.surv.ner$comb_sigs[i] >= ner_q[2][[1]] & tp.surv.ner$comb_sigs[i] <= ner_q[4][[1]]) {
         tp.surv.ner$comb_sigs_hilow[i] <- "med"
     }
     else if (tp.surv.ner$comb_sigs[i] > ner_q[4][[1]]) {
         tp.surv.ner$comb_sigs_hilow[i] <- "high"
     }
     
 }
 tp.surv.ner$comb_sigs_hilow <- factor(tp.surv.ner$comb_sigs_hilow, levels=c("low", "med", "high"))
 
 # Just high and low
 tp.surv.ner$comb_sigs_hilow2 <- tp.surv.ner$comb_sigs
 
 for (i in 1:nrow(tp.surv.ner)) {
     if (tp.surv.ner$comb_sigs[i] < ner_q[2][[1]]) {
         tp.surv.ner$comb_sigs_hilow[i] <- "low"
     }
     else if (tp.surv.ner$comb_sigs[i] >= ner_q[2][[1]] & tp.surv.ner$comb_sigs[i] <= ner_q[4][[1]]) {
         tp.surv.ner$comb_sigs_hilow[i] <- "med"
     }
     else if (tp.surv.ner$comb_sigs[i] > ner_q[4][[1]]) {
         tp.surv.ner$comb_sigs_hilow[i] <- "high"
     }
     
 }
 tp.surv.ner$comb_sigs_hilow <- factor(tp.surv.ner$comb_sigs_hilow, levels=c("low", "med", "high"))
 
 #tp.surv.mmr <- tp.surv.mmr[which(tp.surv.mmr$comb_sigs!=0),] # Remove 150 samples
 tp.surv.mmr$comb_sigs_hilow <- tp.surv.mmr$comb_sigs
 mmr_q <-quantile(tp.surv.mmr$comb_sigs)
 
 for (i in 1:nrow(tp.surv.mmr)) {
     if (tp.surv.mmr$comb_sigs[i] < mmr_q[2][[1]]) {
         tp.surv.mmr$comb_sigs_hilow[i] <- "low"
     }
     else if (tp.surv.mmr$comb_sigs[i] >= mmr_q[2][[1]] & tp.surv.mmr$comb_sigs[i] <= mmr_q[4][[1]]) {
         tp.surv.mmr$comb_sigs_hilow[i] <- "med"
     }
     else if (tp.surv.mmr$comb_sigs[i] > mmr_q[4][[1]]) {
         tp.surv.mmr$comb_sigs_hilow[i] <- "high"
     }
     
 }
 
 # Survival by NER 
 
 tp.ner.trim <- tp.surv.ner[!duplicated(tp.surv.ner$tcga_participant_barcode),]
 tp.ner.fitKM <- survfit(Surv(days_to_death, vital_status)~comb_sigs_hilow,
                     data=tp.ner.trim)
 
 tp.ner.fitKM <- survfit(Surv(days_to_death, vital_status)~comb_sigs_hilow,
                         data=tp.ner.trim)
 summary(tp.fitKM)
 
 tp.ner.sigs.cox <- coxph(Surv(days_to_death, vital_status)~comb_sigs_hilow,
                          data=tp.surv.ner)
 tp.ner.sigs.cox <- coxph(Surv(days_to_death, vital_status)~comb_sigs_hilow,
                          data=tp.ner.trim)
 ner.cox <- summary(tp.ner.sigs.cox)
 
 ner.hr <- round(ner.cox$coefficients[2], digits=2)
 ner.CI1 <- round(ner.cox$conf.int[3], digits=2)
 ner.CI2 <- round(ner.cox$conf.int[4], digits=2)
 ner.pval <- round(ner.cox$sctest[3], digits=4)
 
 pdf("./Figures/TCGA_SKCM_TP_NER_survival_KM_annot.pdf", w=7, h=6)
 plot(tp.ner.fitKM, col=c("red","orange", "blue"),
      mark.time = TRUE,
      xlab="Survival time (days)", ylab="% Overall survival",
      lwd = 3,
      cex.axis = 1.3,
      cex.lab = 1.5)
 legend("topright", 
        #legend.title = "Combined signature\ncontribution",
        legend = c("High (n=1560)", "Medium (n=3250)", "Low (n=1560)"), 
        fill=c("red","orange","blue"), 
        bty="n",
        cex = 1.5)
 #y.intersp=1, x.intersp=2, text.width=0.8)
 legend("bottomleft", 
        legend = c("HR: ", paste0(ner.hr, ", 95% CI: ", ner.CI1, "-", ner.CI2), paste0("p = ", ner.pval)),  
        bty = "n",
        y.intersp=0.8,x.intersp=0.5,text.width=0.1,
        cex = 1.5)
 dev.off()
 
 # Survival by NER 
 tp.mmr.fitKM <- survfit(Surv(days_to_death, vital_status)~comb_sigs_hilow,
                         data=tp.surv.mmr)
 summary(tp.fitKM)
 
 tp.mmr.sigs.cox <- coxph(Surv(days_to_death, vital_status)~comb_sigs_hilow,
                          data=tp.surv.mmr)
 mmr.cox <- summary(tp.mmr.sigs.cox)
 
 mmr.hr <- round(mmr.cox$coefficients[2], digits=2)
 mmr.CI1 <- round(mmr.cox$conf.int[3], digits=2)
 mmr.CI2 <- round(mmr.cox$conf.int[4], digits=2)
 mmr.pval <- round(mmr.cox$sctest[3], digits=4)
 

 
 pdf("./Figures/TCGA_SKCM_TP_MMR_survival_KM_annot.pdf", w=7, h=6)
 plot(tp.mmr.fitKM, col=c("red","orange", "blue"),
      mark.time = TRUE,
      xlab="Survival time (days)", ylab="% Overall survival",
      lwd = 3,
      cex.axis = 1.3,
      cex.lab = 1.5)
 legend("topright", 
        #legend.title=,
        legend = c("High (n=624)", "Medium (n=1924)", "Low (n=0)"), 
        fill=c("red","orange","blue"), 
        bty="n",
        cex=1.5)
 legend("bottomleft", 
        legend = c("HR: ", paste0(mmr.hr, ", 95% CI: ", mmr.CI1, "-", mmr.CI2), paste0("p = ", mmr.pval)), 
        bty = "n",
        y.intersp=0.8,x.intersp=0.5,text.width=0.1,
        cex=1.5)
 dev.off()
 
 
 
#--------- AMPLIFICATION ----------------------------------------------------------------------------------------------

# Read in data
skcm.cna <- read.table(file="./ner_mmr_cna.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
skcm.cna$SAMPLE_ID <- sapply(skcm.cna$SAMPLE_ID, function(x) substr(x, 1, 12))
colnames(skcm.cna)[2] <- "tcga_participant_barcode"
save(skcm.cna, file="./Data/TCGA_SKCM_NER_copynumber.RData")

# 287 samples, TM

amp <- c("tcga_participant_barcode","POLD3", "POLD4", "PARP1", "POLH", "MSH5", "MSH4", "PMS2", "POLE", "EXO1", "WRN")
ner.sigs <- c("S4", "S7", "S11", "S22", "S24", "S29")
mmr.sigs <- c("S6", "S15", "S20", "S26")

skcm.amp <- skcm.cna[,which(colnames(skcm.cna) %in% amp)]
for (i in 2:11) {
    skcm.amp[,i] <- ifelse(skcm.amp[,i]==2, "Amplified", "Not amplified")
}
save(skcm.amp, file="./Output/TCGA_SKCM_NER_amplification.RData")

skcm.amp.mutsigs.metscore <- merge(skcm.amp, tm.skcm.ner.mmr[,c(1,4:34,37:38)], by="tcga_participant_barcode", all.x=TRUE, all.y=FALSE)
skcm.amp.mutsigs.metscore <- skcm.amp.mutsigs.metscore[!duplicated(skcm.amp.mutsigs.metscore),]


# Metscore
for (i in 2:11) {
    amp <- 0
    notamp <- 0
    
    gene <- colnames(skcm.amp.mutsigs.metscore)[i]
    temp <- noquote(gene)
    
    len <- length(table(skcm.amp.mutsigs.metscore[,i]))
    if (len==1) {
        amp <- 0
        notamp <- table(skcm.amp.mutsigs.metscore[,i])[[1]]
    }
    if (len==2) {
        amp <- table(skcm.amp.mutsigs.metscore[,i])[[1]]
        notamp <- table(skcm.amp.mutsigs.metscore[,i])[[2]]
    }
    
    if (amp>0 & notamp>0) {
        print(paste0(i, " ", gene, " enter\n"))
        p <- round(wilcox.test(skcm.amp.mutsigs.metscore$metscore[which(skcm.amp.mutsigs.metscore[,i]=="Amplified")], 
                    skcm.amp.mutsigs.metscore$metscore[which(skcm.amp.mutsigs.metscore[,i]=="Not amplified")])$p.value, digits=4)
        
        if (p<0.05) {
            print(ggplot(skcm.amp.mutsigs.metscore, aes_string(x=temp, y="metscore", fill=temp)) +
                geom_boxplot() +
                scale_fill_manual(values=c(colors3[2], colors3[1])) +
                ggtitle("") +
                xlab(gene) +
                ylab("Metastatic score") +
                #scale_y_continuous(breaks=seq(0,200,25)) +
                #coord_cartesian(ylim=c(0,140), expand=FALSE) +
                theme_bw() +
                scale_x_discrete(labels=c(paste0("Amplified\nn=",amp), paste0("Not amplified\nn=", notamp))) +
                theme(legend.position="none",
                      axis.text = element_text(size=12),
                      axis.title = element_text(size=13)) +
                annotate("text", x=0.82, y=1.8, label=paste0("p=",p))
            )
        }
    }
}

# NER mutsigs

for (j in 1:length(ner.sigs)) {
    
    sig <- ner.sigs[j]
    
    for (i in 2:11) {
        h <- match(sig, colnames(skcm.amp.mutsigs.metscore))
        
        amp <- 0
        notamp <- 0
        
        gene <- colnames(skcm.amp.mutsigs.metscore)[i]
        temp <- noquote(gene)
        
        len <- length(table(skcm.amp.mutsigs.metscore[,i]))
        if (len==1) {
            amp <- 0
            notamp <- table(skcm.amp.mutsigs.metscore[,i])[[1]]
        }
        if (len==2) {
            amp <- table(skcm.amp.mutsigs.metscore[,i])[[1]]
            notamp <- table(skcm.amp.mutsigs.metscore[,i])[[2]]
        }
        
        if (amp>0 & notamp>0) {
            print(paste0(i, " ", gene, " enter\n"))
            p <- round(wilcox.test(skcm.amp.mutsigs.metscore[which(skcm.amp.mutsigs.metscore[,i]=="Amplified"),h], 
                                   skcm.amp.mutsigs.metscore[which(skcm.amp.mutsigs.metscore[,i]=="Not amplified"),h])$p.value, digits=4)
            
            if (p < 0.05) {
                print(ggplot(skcm.amp.mutsigs.metscore, aes_string(x=temp, y=sig, fill=temp)) +
                      geom_boxplot() +
                      scale_fill_manual(values=c(colors3[2], colors3[1])) +
                      ggtitle("") +
                      xlab(gene) +
                      ylab(sig) +
                      #scale_y_continuous(breaks=seq(0,200,25)) +
                      #coord_cartesian(ylim=c(0,140), expand=FALSE) +
                      theme_bw() +
                      scale_x_discrete(labels=c(paste0("Amplified\nn=",amp), paste0("Not amplified\nn=", notamp))) +
                      theme(legend.position="none",
                            axis.text = element_text(size=12),
                            axis.title = element_text(size=13)) +
                      annotate("text", x=0.82, y=1.8, label=paste0("p=",p))
                )
            }
        }
    }
    
}


# MMR mutsigs

for (j in 1:length(mmr.sigs)) {
    
    sig <- mmr.sigs[j]
    
    for (i in 2:11) {
        h <- match(sig, colnames(skcm.amp.mutsigs.metscore))
        
        amp <- 0
        notamp <- 0
        
        gene <- colnames(skcm.amp.mutsigs.metscore)[i]
        temp <- noquote(gene)
        
        len <- length(table(skcm.amp.mutsigs.metscore[,i]))
        if (len==1) {
            amp <- 0
            notamp <- table(skcm.amp.mutsigs.metscore[,i])[[1]]
        }
        if (len==2) {
            amp <- table(skcm.amp.mutsigs.metscore[,i])[[1]]
            notamp <- table(skcm.amp.mutsigs.metscore[,i])[[2]]
        }
        
        if (amp>0 & notamp>0) {
            print(paste0(i, " ", gene, " enter\n"))
            p <- round(wilcox.test(skcm.amp.mutsigs.metscore[which(skcm.amp.mutsigs.metscore[,i]=="Amplified"),h], 
                                   skcm.amp.mutsigs.metscore[which(skcm.amp.mutsigs.metscore[,i]=="Not amplified"),h])$p.value, digits=4)
            
            if (p < 0.05) {
                print(ggplot(skcm.amp.mutsigs.metscore, aes_string(x=temp, y=sig, fill=temp)) +
                          geom_boxplot() +
                          scale_fill_manual(values=c(colors3[2], colors3[1])) +
                          ggtitle("") +
                          xlab(gene) +
                          ylab(sig) +
                          #scale_y_continuous(breaks=seq(0,200,25)) +
                          #coord_cartesian(ylim=c(0,140), expand=FALSE) +
                          theme_bw() +
                          scale_x_discrete(labels=c(paste0("Amplified\nn=",amp), paste0("Not amplified\nn=", notamp))) +
                          theme(legend.position="none",
                                axis.text = element_text(size=12),
                                axis.title = element_text(size=13)) +
                          annotate("text", x=0.82, y=1.8, label=paste0("p=",p))
                )
            }
        }
    }
    
}