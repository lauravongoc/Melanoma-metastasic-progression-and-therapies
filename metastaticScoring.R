# Metastatic Scoring
# DT Laura Vo Ngoc
# Start: 28/05/2018

library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
require(FirebrowseR)


#--------- WD & LOAD FILES --------------------------------------------------------------------------------------------
setwd("C:/Users/rockp/Desktop/UCL/Project/Project-Laura")

load("./Data/TCGA_SKCM_mutations.RData")            # Raw TCGA mutation data
load("./Output/TCGA_SKCM_snvs.RData")               # Selected SNVs

# Primary tumor
load("./Output/TCGA_SKCM_TP_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_SKCM_TP_weights_cut0.00.RData")     # Mutational signatures output

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

# Gene lists
load("./Output/TCGA_SKCM_metastatic_genelist.RData")            # SKCM metastasis genelist
load("./Output/TCGA_SKCM_metastatic_genelist_upreg.RData")      # SKCM metastasis upreg genelist
load("./Output/TCGA_SKCM_metastatic_genelist_downreg.RData")    # SKCM metastasis downreg genelist

# Metastatic score
load("./Output/TCGA_SKCM_TP_metastatic_score_mel.RData")  # SKCM TP
load("./Output/TCGA_SKCM_TM_metastatic_score_mel.RData")  # SKCM TM
load("./Output/TCGA_UVM_TP_metastatic_score_mel.RData")   # UVM TP

# Mutsigs, metscores, and stage
load("./Output/TCGA_SKCM_TP_metastatic_score_stage.RData")  # SKCM TP
load("./Output/TCGA_SKCM_TM_metastatic_score_stage.RData")  # SKCM TM
load("./Output/TCGA_UVM_TP_metastatic_score_stage.RData")     # UVM TP

# Total mutsig mutations
load("./Output/TCGA_SKCM_TP_totalmutsigs.RData")
load("./Output/TCGA_SKCM_TM_totalmutsigs.RData")

#--------- COLORS --------------------------------------------------------------------------------------------
colors <- c(brewer.pal(12, "Paired"), brewer.pal(11, "Set3"),"#000000")
colors2 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#BFA89E", "#826C7F", "#FF7F00",
             "#CAB2D6", "#8B786D", "#6A3D9A", "#D6DBD2", "#FFFF99", "#B15928", "#8DD3C7", "#FFFFB3", "#E75A7C", "#BEBADA",
             "#FB8072", "#40476D", "#80B1D3", "#FDB462", "#B3DE69", "#258EA6", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", 
             "#000000")
colors3 <- brewer.pal(3, "Paired")
colors4 <- c(brewer.pal(4, "Paired"),"#000000")
colors6 <- c("#fb9a99", brewer.pal(4, "Paired"),"#000000")
colors7 <- c("#1F78B4", "#B2DF8A", "#33A02C", "#000000")


#--------- METASTATIS GENES -------------------------------------------------------------------------------------------

# Download expression data for metastasis genes
    # cBioPortal
genes_cbp <- c("MMP1","MMP2","MMP3","MMP7","MMP9","MMP10","MMP11","MMP12","MMP13","MMP14","MMP15","MMP16","MMP17",
               "MMP19","MMP21","MMP23B","MMP24","MMP25","MMP26","MMP27","MMP28","ITGB3","ITGAV","PTK2","CDH1","SPARC",
               "WFDC2")

    # (Rikers,2008) + table OneNote
genes1 <- c("NEDD9","TWIST1","HIF1A","EPAS1", "C7", "AKT3", "BCL2A1", "BIRC5", "BUB1", "C1orf90", "CDC45L", "CDC6", 
           "CDK2", "CSAG2", "DUSP4", "DUSP6", "GDF15", "GPR19", "GYPC", "HEY1", "HOXA10", "HOXB6", "HOXB7",
           "HOXB9", "KIFC1", "MAGEA1", "MAGEA2", "MAGEA3", "MAGEA5", "MAGEA6", "MAGEA12", "MMP14", "PEG10", 
           "RASGRF1", "RGS20", "RRM2", "SLC16A4", "SOX5", "SPP1", "SPRED1", "TRIM51", "TYMS")

    # (Qiu, 2015)
genes2 <- c("AMDHD2", "BRCA1", "CDK1", "CCNB1", "CCNB2", "CCNE1", "CCNE2", "CDC20", "CDK1", "CHEK1", 
           "COL11A1", "DNA2", "EIF2AK3", "FPGT", "GALK1", "GFPT1", "HMGA2", "MAD2L1", "MAD2L1BP", "MAGEA3", 
           "MAGEA6", "MCM6", "MCM7", "MSH2", "MSH6", "PCNA", "PGM3", "POLA2", "POLE2", "POLR2K", "POLR3D", 
           "PRIM1", "PRIM2", "RFC2", "RFC3", "RFC4", "SEC23A", "SEC23B", "SEC24A", "SEC24B", "SEC24D", "SMC2", 
           "SMC4", "SPC25", "UAP1", "UBE2C", "UBE2D", "UBE2S", "UBE2W")

    # UVM genes OneNote
genes_uvm <- c("HTR2B","ECM1","RAB31","CDH1","PRAME","FXR1","LTA4H","EIF1B","ID2","ROBO1","LMCD1","SATB1","MTUS1")

tp.skcm.expr.cbp <- Samples.mRNASeq(format = "csv",
                           gene = genes_cbp,
                           cohort = "SKCM",
                           sample_type = "TP",
                           protocol = "RSEM",
                           page_size=6000
)
tp.skcm.expr1 <- Samples.mRNASeq(format = "csv",
                                    gene = genes1,
                                    cohort = "SKCM",
                                    sample_type = "TP",
                                    protocol = "RSEM",
                                    page_size=6000
)
tp.skcm.expr2 <- Samples.mRNASeq(format = "csv",
                                    gene = genes2,
                                    cohort = "SKCM",
                                    sample_type = "TP",
                                    protocol = "RSEM",
                                    page_size=6000
)
tp.skcm.expr <- rbind(tp.skcm.expr.cbp, tp.skcm.expr1, tp.skcm.expr2)
tp.skcm.expr <- tp.skcm.expr[!duplicated(tp.skcm.expr),]
save(tp.skcm.expr, file="./Data/TCGA_SKCM_TP_expression_mel.RData")

tm.skcm.expr.cbp <- Samples.mRNASeq(format = "csv",
                                gene = genes_cbp,
                                cohort = "SKCM",
                                sample_type = "TM",
                                protocol = "RSEM",
                                page_size=6000
)
tm.skcm.expr1 <- Samples.mRNASeq(format = "csv",
                                gene = genes1,
                                cohort = "SKCM",
                                sample_type = "TM",
                                protocol = "RSEM",
                                page_size=6000
)
tm.skcm.expr2 <- Samples.mRNASeq(format = "csv",
                                gene = genes2,
                                cohort = "SKCM",
                                sample_type = "TM",
                                protocol = "RSEM",
                                page_size=6000
)
tm.skcm.expr <- rbind(tm.skcm.expr.cbp, tm.skcm.expr1, tm.skcm.expr2)
tm.skcm.expr <- tm.skcm.expr[!duplicated(tm.skcm.expr),]
save(tm.skcm.expr, file="./Data/TCGA_SKCM_TM_expression_mel.RData")

uvm.expr <- Samples.mRNASeq(format = "csv",
                                gene = genes_uvm,
                                cohort = "UVM",
                                sample_type = "TP",
                                protocol = "RSEM",
                                page_size=6000
)
save(uvm.expr, file="./Data/TCGA_UVM_TP_expression_mel.RData")

# Generate matrix of participant and gene expr
tp.skcm.expr$expression_log2 <- as.numeric(tp.skcm.expr$expression_log2)
mat.tp.skcm.expr <- dcast(tcga_participant_barcode~gene,data=tp.skcm.expr,
                value.var="expression_log2", fun.aggregate=mean)
rownames(mat.tp.skcm.expr) <- mat.tp.skcm.expr$tcga_participant_barcode
mat.tp.skcm.expr <- mat.tp.skcm.expr[,-1]

tm.skcm.expr$expression_log2 <- as.numeric(tm.skcm.expr$expression_log2)
mat.tm.skcm.expr <- dcast(tcga_participant_barcode~gene,data=tm.skcm.expr,
                         value.var="expression_log2", fun.aggregate=mean)
rownames(mat.tm.skcm.expr) <- mat.tm.skcm.expr$tcga_participant_barcode
mat.tm.skcm.expr <- mat.tm.skcm.expr[,-1]

uvm.expr$expression_log2 <- as.numeric(uvm.expr$expression_log2)
mat.uvm.expr <- dcast(tcga_participant_barcode~gene,data=uvm.expr,
                         value.var="expression_log2", fun.aggregate=mean)
rownames(mat.uvm.expr) <- mat.uvm.expr$tcga_participant_barcode
mat.uvm.expr <- mat.uvm.expr[,-1]


mat.tp.skcm.expr[is.na(mat.tp.skcm.expr)] <- 0 
mat.tm.skcm.expr[is.na(mat.tm.skcm.expr)] <- 0 
mat.uvm.expr[is.na(mat.uvm.expr)] <- 0

save(mat.tp.skcm.expr, file="./Output/TCGA_SKCM_TP_expression_mel_matrix.RData")
save(mat.tm.skcm.expr, file="./Output/TCGA_SKCM_TM_expression_mel_matrix.RData")
save(mat.uvm.expr, file="./Output/TCGA_UVM_TP_expression_mel_matrix.RData")

# Heatmap overview of each cohort
tp.skcm <- as.data.frame(mat.tp.skcm.expr)
tp.skcm[is.na(tp.skcm)] <- 0                        # Replace NA with 0
pdf("./Figures/TCGA_SKCM_TP_expression_mel_heatmap.pdf", w=10, h=6, onefile=FALSE)
tp.skcm.map <-pheatmap(t(tp.skcm), show_colnames = FALSE, main="Metastatis genes expression SKCM TP")
dev.off()

tm.skcm <- as.data.frame(mat.tm.skcm.expr)
tm.skcm[is.na(tm.skcm)] <- 0                        # Replace NA with 0
pdf("./Figures/TCGA_SKCM_TM_expression_mel_heatmap.pdf", w=10, h=6, onefile=FALSE)
tm.skcm.map <-pheatmap(t(tm.skcm), show_colnames = FALSE, main="Metastatis genes expression SKCM TM")
dev.off()

uvm <- as.data.frame(mat.uvm.expr)
uvm[is.na(uvm)] <- 0                                # Replace NA with 0
pdf("./Figures/TCGA_UVM_TP_expression_mel_heatmap.pdf", w=10, h=6, onefile=FALSE)
uvm.map <-pheatmap(t(uvm), show_colnames = FALSE, main="Metastatis genes expression UVM TP")
dev.off()


#--------- GENE SELECTION (LM): SKCM TP VS TM ----------------------------------------------------------------------

mat.tp.skcm.expr$sample_type <- "TP"
mat.tm.skcm.expr$sample_type <- "TM"
mat.skcm.expr <- rbind(mat.tp.skcm.expr, mat.tm.skcm.expr)

pval.skcm.expr <- data.frame(Gene=noquote(colnames(mat.skcm.expr)[1:110]), pval=NA, sign=NA)
for (i in 1:110) {
    mod <- lm(mat.skcm.expr[,i] ~ mat.skcm.expr$sample_type)
    pval.skcm.expr[i,2] <- coef(summary(mod))[2,4]
}
pval.skcm.expr$sign <- ifelse(pval.skcm.expr$pval<0.05, 1, 0)

genelist <- pval.skcm.expr[which(pval.skcm.expr$sign==1),1]
save(genelist, file="./Output/TCGA_SKCM_metastatic_genelist.RData")


#--------- GENE SELECTION (TTEST): SKCM TP VS TM --------------------------------------------------------------------

# Get genes with significant difference in means
# Data frame for t-test output
skcm.expr.ttest <- data.frame(Gene=colnames(mat.tp.skcm.expr), t=NA, df=NA, pval = NA)
for (i in 1:110) {
    gene <- noquote(colnames(mat.tp.skcm.expr[i]))
    
    tp <- na.omit(mat.tp.skcm.expr[i])
    tm <- na.omit(mat.tm.skcm.expr[i])
    ttest <- t.test(mat.tp.skcm.expr[,i], mat.tm.skcm.expr[,i])
    
    skcm.expr.ttest[which(skcm.expr.ttest$Gene==gene),2] <- ttest$statistic
    skcm.expr.ttest[which(skcm.expr.ttest$Gene==gene),3] <- ttest$parameter
    skcm.expr.ttest[which(skcm.expr.ttest$Gene==gene),4] <- ttest$p.value
}
skcm.expr.ttest$sign <- NA
skcm.expr.ttest$sign <- ifelse(skcm.expr.ttest$pval<0.05, 1, 0)         # put 1 if significant
write.csv(skcm.expr.ttest, file="./Output/TCGA_SKCM_expression_eachgenettest.csv")

genelist <- skcm.expr.ttest[which(skcm.expr.ttest$sign==1), 1]           # 63 genes significantly different

# Generate variables containing only data for these genes
tp.skcm.genelist.expr <- mat.tp.skcm.expr[,which(colnames(mat.tp.skcm.expr)%in%genelist)]
tm.skcm.genelist.expr <- mat.tm.skcm.expr[,which(colnames(mat.tm.skcm.expr)%in%genelist)]

# Heatmap overview of each cohort
tp.skcm <- as.data.frame(tp.skcm.genelist.expr)
tp.skcm[is.na(tp.skcm)] <- 0                        # Replace NA with 0
pdf("./Figures/TCGA_SKCM_TP_expression_mel_heatmap.pdf", w=10, h=8, onefile=FALSE)
tp.skcm.map <-pheatmap(t(tp.skcm), show_colnames = FALSE, main="Metastatis genes expression SKCM TP")
dev.off()

tm.skcm <- as.data.frame(tm.skcm.genelist.expr)
tm.skcm[is.na(tm.skcm)] <- 0                        # Replace NA with 0
pdf("./Figures/TCGA_SKCM_TM_expression_mel_heatmap.pdf", w=10, h=8, onefile=FALSE)
tm.skcm.map <-pheatmap(t(tm.skcm), show_colnames = FALSE, main="Metastatis genes expression SKCM TM")
dev.off()


#--------- MEAN GENE EXPRESSION (genelist) ---------------------------------------------------------------------------

tp.skcm.expr.mean <- data.frame(gene=colnames(tp.skcm.genelist.expr), SKCM_TP=colMeans(tp.skcm.genelist.expr, na.rm=TRUE))
tm.skcm.expr.mean <- data.frame(gene=colnames(tm.skcm.genelist.expr), SKCM_TM=colMeans(tm.skcm.genelist.expr, na.rm=TRUE))
#uvm.expr.mean <- data.frame(gene=colnames(mat.uvm.expr), UVM=colMeans(mat.uvm.expr, na.rm=TRUE))

skcm.expr.mean <- merge(tp.skcm.expr.mean, tm.skcm.expr.mean, by="gene", all=TRUE)      # merge SKCM cohorts
#merged <- merge(merged, uvm.expr.mean, by="gene", all=TRUE)                     # add UVM 
skcm.expr.mean[is.na(skcm.expr.mean)] <- 0

order.updown <- c(as.vector(genelist.up), as.vector(genelist.down))
skcm.expr.mean$gene <- factor(skcm.expr.mean$gene, levels=order.updown)
skcm.expr.mean <- skcm.expr.mean[match(order.updown, skcm.expr.mean$gene),]

melt.skcm.expr.mean <- melt(skcm.expr.mean)


pdf("./Figures/TCGA_SKCM_expression_avg.pdf", w=13, h=6)
ggplot(melt.skcm.expr.mean, aes(x=Gene, y=mean_val)) + 
    scale_fill_manual(values=colors3) +
    geom_bar(stat="identity", aes(fill=Cohort), width=0.8, position = position_dodge(width = 0.8)) +
    ggtitle("Mean gene expression SKCM TP vs. TM") +
    xlab("Gene") +
    ylab("Mean gene expression (log2)") +
    scale_y_continuous(breaks=seq(0,100,1)) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1), legend.position=c(0.96,0.95)) +
    coord_cartesian(ylim=c(-0.1,15), expand=FALSE)
dev.off()


skcm.expr.mean$upreg <- NA
skcm.expr.mean$upreg <- ifelse(skcm.expr.mean$SKCM_TM > skcm.expr.mean$SKCM_TP, 1, 0)
write.csv(skcm.expr.mean, file="./Output/TCGA_SKCM_expression_mean.csv")

genelist.up <- skcm.expr.mean[which(skcm.expr.mean$upreg==1), 1]        # 42 genes > in TM
genelist.down <- skcm.expr.mean[which(skcm.expr.mean$upreg==0), 1]      # 21 genes < in TM
save(genelist.up, file="./Output/TCGA_SKCM_metastatic_genelist_upreg.RData")
save(genelist.down, file="./Output/TCGA_SKCM_metastatic_genelist_downreg.RData")


#--------- METASTATIC SCORE ----------------------------------------------------------------------------------------
uvm.genelist.up <- c("HTR2B", "ECM1", "RAB31", "CDH1","PRAME")
uvm.genelist.down <- c("FXR1","LTA4H","EIF1B","ID2","ROBO1","LMCD1","SATB1","MTUS1")


# Order matrix by upreg then downreg
tp.skcm.genelist.expr <- mat.tp.skcm.expr[,match(order.updown, colnames(mat.tp.skcm.expr))]
tm.skcm.genelist.expr <- mat.tm.skcm.expr[,match(order.updown, colnames(mat.tm.skcm.expr))]
uvm.genelist.expr <- mat.uvm.expr[,match(genes_uvm, colnames(mat.uvm.expr))]

# Metastatic potential score per sample = avg expression upreg genes / avg expression downreg genes 
tp.skcm.metscore <- data.frame(tcga_participant_barcode=noquote(rownames(tp.skcm.genelist.expr)),
                               cohort="SKCM TP", 
                               metscore=(rowMeans(tp.skcm.genelist.expr[,1:41]/rowMeans(tp.skcm.genelist.expr[,42:63]), 
                                              na.rm=TRUE)))
tm.skcm.metscore <- data.frame(tcga_participant_barcode=noquote(rownames(tm.skcm.genelist.expr)),
                               cohort = "SKCM TM",
                               metscore=(rowMeans(tm.skcm.genelist.expr[,1:41]/rowMeans(tm.skcm.genelist.expr[,42:63]), 
                                                      na.rm=TRUE)))
uvm.metscore <- data.frame(tcga_participant_barcode=noquote(rownames(mat.uvm.expr)),
                           cohort="UVM TP",
                           metscore=(rowMeans(uvm.genelist.expr[,1:5]/rowMeans(uvm.genelist.expr[,6:13]),  
                                                   na.rm=TRUE)))

save(tp.skcm.metscore, file="./Output/TCGA_SKCM_TP_metastatic_score_mel.RData")
save(tm.skcm.metscore, file="./Output/TCGA_SKCM_TM_metastatic_score_mel.RData")
save(uvm.metscore, file="./Output/TCGA_UVM_TP_metastatic_score_mel.RData")

# Overlapping histogram of metscore distributions
pdf("./Figures/TCGA_SKCM_UVM_metscore_hist_small_cutoff.pdf", w=8, h=6)
hist(tm.skcm.metscore$metscore, 
     col="#1F78B4B3",  
     main="Metastatic score distribution", 
     xlab="Metastatic score",
     xlim=c(0.4,1.8),
     ylim=c(0,100))
hist(tp.skcm.metscore$metscore, 
     col="#fb9a99B3", 
     main="Metastatic score distribution SKCM TP", 
     xlab="Metastatic score",
     xlim=c(0.4,1.8),
     ylim=c(0,100),
     add=TRUE)
hist(uvm.metscore$metscore,
     col="#B2DF8AB3", 
     main="Metastatic score distribution UVM TP", 
     xlab="Metastatic score",
     breaks=14,
     xlim=c(0.4,1.8),
     ylim=c(0,100),
     add=TRUE)
legend("topright", c("SKCM TP", "SKCM TM", "UVM TP"), col=c("#fb9a99B3","#1F78B4B3", "#B2DF8AB3"), lwd=10, bty="n")
abline(v=1, col="black", lty=2, lwd=2)
dev.off()


#--------- ~ Cutoff ---------

# Cutoff = histogram intersects ~= 1.0

tp.skcm.metscore$met_potential <- NA
tp.skcm.metscore$met_potential <- ifelse(tp.skcm.metscore$metscore > 1, "high", "low")

tm.skcm.metscore$met_potential <- NA
tm.skcm.metscore$met_potential <- ifelse(tm.skcm.metscore$metscore > 1, "high", "low")

uvm.metscore$met_potential <- NA
uvm.metscore$met_potential <- ifelse(uvm.metscore$metscore > 1, "high", "low")

# High vs. Low by quantiles
tp.skcm.metscore$met_potential_quant <- NA
for (i in 1:nrow(tp.skcm.metscore)) {
    if (tp.skcm.metscore$metscore[i] < quantile(tp.skcm.metscore$metscore)[2]) {
        tp.skcm.metscore$met_potential_quant[i] <- "low"
    }
    else if (tp.skcm.metscore$metscore[i] > quantile(tp.skcm.metscore$metscore)[4]) {
        tp.skcm.metscore$met_potential_quant[i] <- "high"
    }
    else {
        tp.skcm.metscore$met_potential_quant[i] <- "intermediate"
    }
    
}
tm.skcm.metscore$met_potential_quant <- NA
for (i in 1:nrow(tm.skcm.metscore)) {
    if (tm.skcm.metscore$metscore[i] < quantile(tm.skcm.metscore$metscore)[2]) {
        tm.skcm.metscore$met_potential_quant[i] <- "low"
    }
    else if (tm.skcm.metscore$metscore[i] > quantile(tm.skcm.metscore$metscore)[4]) {
        tm.skcm.metscore$met_potential_quant[i] <- "high"
    }
    else {
        tm.skcm.metscore$met_potential_quant[i] <- "intermediate"
    }
}

# Create table of low vs high potential counts with early or late disease progression SKCM
earlyLate_hilow <- data.frame(data=c("early","late"), low=NA, high=NA)
rownames(earlyLate_hilow) <- earlyLate_hilow[,1] 
earlyLate_hilow <- earlyLate_hilow[,-1]
earlyLate_hilow[1,1] <- nrow(tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$earlyLate == "early" & tp.skcm.metscore.stage$met_potential == "low"),]) + 
    nrow(tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate == "early" & tm.skcm.metscore.stage$met_potential == "low"),])
earlyLate_hilow[1,2] <- nrow(tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$earlyLate == "early" & tp.skcm.metscore.stage$met_potential == "high"),]) + 
    nrow(tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate == "early" & tm.skcm.metscore.stage$met_potential == "high"),])
earlyLate_hilow[2,1] <- nrow(tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$earlyLate == "late" & tp.skcm.metscore.stage$met_potential == "low"),]) + 
    nrow(tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate == "late" & tm.skcm.metscore.stage$met_potential == "low"),])
earlyLate_hilow[2,2] <- nrow(tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$earlyLate == "late" & tp.skcm.metscore.stage$met_potential == "high"),]) + 
    nrow(tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate == "late" & tm.skcm.metscore.stage$met_potential == "high"),])

# Fisher's exact test SKCM
fisher.test(earlyLate_hilow)

# Create table of low vs high potential counts with early or late disease progression UVM
uvm_earlyLate_hilow <- data.frame(data=c("early","late"), low=NA, high=NA)
rownames(uvm_earlyLate_hilow) <- uvm_earlyLate_hilow[,1] 
uvm_earlyLate_hilow <- uvm_earlyLate_hilow[,-1]
uvm_earlyLate_hilow[1,1] <- nrow(uvm.metscore.stage[which(uvm.metscore.stage$earlyLate == "early" & uvm.metscore.stage$met_potential == "low"),])
uvm_earlyLate_hilow[1,2] <- nrow(uvm.metscore.stage[which(uvm.metscore.stage$earlyLate == "early" & uvm.metscore.stage$met_potential == "high"),])
uvm_earlyLate_hilow[2,1] <- nrow(uvm.metscore.stage[which(uvm.metscore.stage$earlyLate == "late" & uvm.metscore.stage$met_potential == "low"),])
uvm_earlyLate_hilow[2,2] <- nrow(uvm.metscore.stage[which(uvm.metscore.stage$earlyLate == "late" & uvm.metscore.stage$met_potential == "high"),])

# Fisher's exact test UVM
fisher.test(uvm_earlyLate_hilow)


# low (0-0.8), med (0.8-1), high(>1)
tp.skcm.metscore$met_potential_med <- NA
tp.skcm.metscore$met_potential_med <- ifelse(tp.skcm.metscore$metscore < 1 & tp.skcm.metscore$metscore > 0.8, "med", tp.skcm.metscore$met_potential)

tm.skcm.metscore$met_potential_med <- NA
tm.skcm.metscore$met_potential_med <- ifelse(tm.skcm.metscore$metscore < 1 & tm.skcm.metscore$metscore > 0.8, "med", tm.skcm.metscore$met_potential)



#--------- //SKCM TP, TM, VS. UVM ----------------------------------------------------------------------------------------

metscore <- rbind(tp.skcm.metscore, tm.skcm.metscore)
#metscore <- rbind(tp.skcm.metscore, tm.skcm.metscore, uvm.metscore)
order2 <- factor(metscore$cohort, levels=c("SKCM TP", "SCKM TM", "UVM"))

pdf("./Figures/TCGA_SKCM_metscore_mel_boxplot.pdf", w=8, h=6)
ggplot(metscore, aes(x=cohort, y=metscore, fill=cohort)) + 
    geom_boxplot() +
    scale_fill_manual(values=colors4) +
    ggtitle("Metastatic score SKCM TP vs. TM") +
    xlab("Cohort") +
    ylab("Metastatic potential score") +
    scale_y_continuous(breaks=seq(0,10,0.1)) +
    coord_cartesian(ylim=c(0,2), expand = FALSE )
dev.off()

t.test(tp.skcm.metscore$metscore, tm.skcm.metscore$metscore)


metscore.anova <- aov(metscore~cohort,  data=metscore)
summary(metscore.anova)

compare <- data.frame(c1=c("SKCM TP","SKCM TP","SKCM TM"), c2=c("SKCM TM","UVM","UVM"), pval=NA)
compare[1,3] <- t.test(tp.skcm.metscore$met_score, tm.skcm.metscore$met_score)$p.value
compare[2,3] <- t.test(tp.skcm.metscore$met_score, uvm.metscore$met_score)$p.value
compare[3,3] <- t.test(tm.skcm.metscore$met_score, uvm.metscore$met_score)$p.value

compare$padj <- p.adjust(compare$pval, method="BH")
compare

#--------- CLINICAL DATA ----------------------------------------------------------------------------------------

#--------- ~ Metscore vs. stage ---------
tp.skcm.metscore.stage <- merge(tp.stage, tp.skcm.metscore, by="tcga_participant_barcode")
tm.skcm.metscore.stage <- merge(tm.stage, tm.skcm.metscore, by="tcga_participant_barcode")
uvm.metscore.stage <- merge(stage, uvm.metscore, by="tcga_participant_barcode")

save(tp.skcm.metscore.stage, file="./Output/TCGA_SKCM_TP_metastatic_score_stage.RData")
save(tm.skcm.metscore.stage, file="./Output/TCGA_SKCM_TM_metastatic_score_stage.RData")
save(uvm.metscore.stage, file="./Output/TCGA_UVM_TP_metastatic_score_stage.RData")

# remove stage NA
tp.skcm.metscore.stage <- tp.skcm.metscore.stage[!is.na(tp.skcm.metscore.stage$stage),]            # 4 NA removed
tm.skcm.metscore.stage <- tm.skcm.metscore.stage[!is.na(tm.skcm.metscore.stage$stage),]            # 34 NA removed
uvm.metscore.stage <- uvm.metscore.stage[!is.na(uvm.metscore.stage$stage),]                 # 0 NA removed


# Boxplot by stage
tp.st <- ggplot(tp.skcm.metscore.stage, aes(x=stage, y=metscore, fill=stage)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("SKCM TP") +
    xlab("Cancer stage") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian( ylim=c(0.4,1.8), expand = TRUE) +
    theme_bw()+
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=13)) +
    scale_x_discrete(labels=c("I\nn=2", "II\nn=66", "III\nn=28", "IV\nn=3"))

tm.st <- ggplot(tm.skcm.metscore.stage, aes(x=stage, y=metscore, fill=stage)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("SKCM TM") +
    xlab("Cancer stage") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian( ylim=c(0.4,1.8), expand = TRUE) + 
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=13)) +
    scale_x_discrete(labels=c("I\nn=75", "II\nn=73", "III\nn=141", "IV\nn=20"))


#pdf("./Figures/TCGA_UVM_TP_metscore_mel_stage_boxplot.pdf", w=8, h=6)
uvm.st <- ggplot(uvm.metscore.stage, aes(x=stage, y=metscore, fill=stage)) + 
    geom_boxplot() +
    scale_fill_manual(values=colors7) + 
    geom_boxplot() +
    ggtitle("UVM TP") +
    xlab("Cancer stage") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian( ylim=c(0.4,1.8), expand = TRUE) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=13)) +
    scale_x_discrete(labels=c("II\nn=36", "III\nn=40", "IV\nn=4")) +
    annotate("text", x=1, y=1.8, label=paste0("p=",uvm.anov.p))


# Anova
summary(tp.skcm.metscore.stage.anova <- aov(metscore~stage,  data=tp.skcm.metscore.stage))
summary(tm.skcm.metscore.stage.anova <- aov(metscore~stage,  data=tm.skcm.metscore.stage))
uvm.anova <- summary(uvm.metscore.stage.anova <- aov(metscore~stage,  data=uvm.metscore.stage))
uvm.anov.p <- round(uvm.anova[[1]][5][[1]][1], digits=4)

# Boxplot by disease progression
tp.pr <- ggplot(tp.skcm.metscore.stage, aes(x=earlyLate, y=metscore, fill=earlyLate)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("") +
    xlab("Disease progression") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian( ylim=c(0.4,1.8), expand = TRUE) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=13)) +
    scale_x_discrete(labels=c("Early\nn=68","Late\nn=31"))

tm.pr <- ggplot(tm.skcm.metscore.stage, aes(x=earlyLate, y=metscore, fill=earlyLate)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("") +
    xlab("Disease progression") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian( ylim=c(0.4,1.8), expand = TRUE) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=13)) +
    scale_x_discrete(labels=c("Early\nn=148","Late\nn=161"))

uvm.pr <- ggplot(uvm.metscore.stage, aes(x=earlyLate, y=metscore, fill=earlyLate)) + 
    geom_boxplot() +
    scale_fill_manual(values=colors4) +
    ggtitle("") +
    xlab("Disease progression") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian( ylim=c(0.4,1.8), expand = TRUE) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=13)) +
    scale_x_discrete(labels=c("Early\nn=36","Late\nn=44")) +
    annotate("text", x=0.82, y=1.8, label=paste0("p=",uvm.p))

# Wilcoxon test
tp.p <- round(wilcox.test(tp.skcm.metscore.stage$metscore[which(tp.skcm.metscore.stage$earlyLate=="early")], tp.skcm.metscore.stage$metscore[which(tp.skcm.metscore.stage$earlyLate=="late")])$p.value, digits=4)
tm.p <- round(wilcox.test(tm.skcm.metscore.stage$metscore[which(tm.skcm.metscore.stage$earlyLate=="early")], tm.skcm.metscore.stage$metscore[which(tm.skcm.metscore.stage$earlyLate=="late")])$p.value, digits=4)
uvm.p <- round(wilcox.test(uvm.metscore.stage$metscore[which(uvm.metscore.stage$earlyLate=="early")], uvm.metscore.stage$metscore[which(uvm.metscore.stage$earlyLate=="late")])$p.value, digits=4)


# Output all boxplots asone plot
pdf("./Figures/TCGA_metscore_mel_stage_prog_boxplot2.pdf", w=8, h=6)
ggarrange(tp.st, tm.st, uvm.st, tp.pr, tm.pr, uvm.pr, 
          labels = c("a)", "b)", "c)", "", "", ""),
          ncol = 3, nrow = 2)
dev.off()



#--------- ~ Met potential vs. stage ---------

tp.skcm.metpotential.low <- tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$met_potential=="low"),c(1,35,36,38,39)]
melt.tp.skcm.metpotential.low <- melt(tp.skcm.metpotential.low)
colnames(melt.tp.skcm.metpotential.low) <- c("Sample", "Stage", "earlyLate", "Potential", "variable", "value")

tp.skcm.metpotential.high <- tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$met_potential=="high"),c(1,35,36,38,39)]
melt.tp.skcm.metpotential.high <- melt(tp.skcm.metpotential.high)
colnames(melt.tp.skcm.metpotential.high) <- c("Sample", "Stage", "earlyLate", "Potential", "variable", "value")

pdf("./Figures/TCGA_SKCM_TP_metpotential_low_stage.pdf", w=6, h=6)
ggplot(melt.tp.skcm.metpotential.low, aes(x=Stage, y=value)) + 
    scale_fill_manual(values=colors3) +
    geom_boxplot() +
    ggtitle("Low metastatic potential SKCM TP by stage") +
    xlab("Stage") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,100,0.1)) +
    coord_cartesian(ylim=c(0,1.5), expand=FALSE)
dev.off()

pdf("./Figures/TCGA_SKCM_TP_metpotential_high_stage.pdf", w=6, h=6)
ggplot(melt.tp.skcm.metpotential.high, aes(x=Stage, y=value)) + 
    scale_fill_manual(values=colors3) +
    geom_boxplot() +
    ggtitle("High metastatic potential SKCM TP by stage") +
    xlab("Stage") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,100,0.1)) +
    coord_cartesian(ylim=c(0,1.7), expand=FALSE)
dev.off()

# Anova
tp.skcm.metpotential.low.anova <- aov(metscore~stage,  data=tp.skcm.metpotential.low)
summary(tp.skcm.metpotential.low.anova)

tp.skcm.metpotential.high.anova <- aov(metscore~stage,  data=tp.skcm.metpotential.high)
summary(tp.skcm.metpotential.high.anova)

# Pairwise t-test SKCM TP low metpotential
stage1 <- tp.skcm.metpotential.low[which(tp.skcm.metpotential.low$stage=="I"),4]
stage2 <- tp.skcm.metpotential.low[which(tp.skcm.metpotential.low$stage=="II"),4]
stage3 <- tp.skcm.metpotential.low[which(tp.skcm.metpotential.low$stage=="III"),4]

compare1 <- data.frame(c1=c("Stage I","Stage I","Stage II"), c2=c("Stage II","Stage III","Stage III"), pval=NA)
compare1[1,3] <- t.test(stage1, stage2)$p.value
compare1[2,3] <- t.test(stage1, stage3)$p.value
compare1[3,3] <- t.test(stage2, stage3)$p.value

compare1$padj <- p.adjust(compare1$pval, method="BH")
compare1


# SKCM TM

tm.skcm.metpotential.low <- tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$met_potential=="low"),c(1,35,36,38,39)]
melt.tm.skcm.metpotential.low <- melt(tm.skcm.metpotential.low)
colnames(melt.tm.skcm.metpotential.low) <- c("Sample", "Stage", "earlyLate", "Potential", "variable", "value")

tm.skcm.metpotential.high <- tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$met_potential=="high"),c(1,35,36,38,39)]
melt.tm.skcm.metpotential.high <- melt(tm.skcm.metpotential.high)
colnames(melt.tm.skcm.metpotential.high) <- c("Sample", "Stage", "earlyLate", "Potential", "variable", "value")

pdf("./Figures/TCGA_SKCM_TM_metpotential_low_stage.pdf", w=6, h=6)
ggplot(melt.tm.skcm.metpotential.low, aes(x=Stage, y=value)) + 
    scale_fill_manual(values=colors3) +
    geom_boxplot() +
    ggtitle("Low metastatic potential SKCM TM by stage") +
    xlab("Stage") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,100,0.1))
#coord_cartesian(ylim=c(0,75), expand=FALSE)
dev.off()

pdf("./Figures/TCGA_SKCM_TM_metpotential_high_stage.pdf", w=6, h=6)
ggplot(melt.tm.skcm.metpotential.high, aes(x=Stage, y=value)) + 
    scale_fill_manual(values=colors3) +
    geom_boxplot() +
    ggtitle("High metastatic potential SKCM TM by stage") +
    xlab("Stage") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,100,0.1))
#coord_cartesian(ylim=c(0,75), expand=FALSE)
dev.off()

# Anova
tm.skcm.metpotential.low.anova <- aov(metscore~stage,  data=tm.skcm.metpotential.low)
summary(tm.skcm.metpotential.low.anova)

tm.skcm.metpotential.high.anova <- aov(metscore~stage,  data=tm.skcm.metpotential.high)
summary(tm.skcm.metpotential.high.anova)


# UVM
uvm.metpotential.low <- uvm.metscore.stage[which(uvm.metscore.stage$met_potential=="low"),c(1,35,36,38,39)]
melt.uvm.metpotential.low <- melt(uvm.metpotential.low)
colnames(melt.uvm.metpotential.low) <- c("Sample", "Stage", "earlyLate", "Potential", "variable", "value")

uvm.metpotential.high <- uvm.metscore.stage[which(uvm.metscore.stage$met_potential=="high"),c(1,35,36,38,39)]
melt.uvm.metpotential.high <- melt(uvm.metpotential.high)
colnames(melt.uvm.metpotential.high) <- c("Sample", "Stage", "earlyLate", "Potential", "variable", "value")

pdf("./Figures/TCGA_UVM_TP_metpotential_low_stage.pdf", w=6, h=6)
ggplot(melt.uvm.metpotential.low, aes(x=Stage, y=value)) + 
    scale_fill_manual(values=colors3) +
    geom_boxplot() +
    ggtitle("Low metastatic potential UVM TP by stage") +
    xlab("Stage") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,100,0.1)) +
    coord_cartesian(ylim=c(0,1), expand=FALSE)
dev.off()

pdf("./Figures/TCGA_UVM_TP_metpotential_high_stage.pdf", w=6, h=6)
ggplot(melt.uvm.metpotential.high, aes(x=Stage, y=value)) + 
    scale_fill_manual(values=colors3) +
    geom_boxplot() +
    ggtitle("High metastatic potential SKCM TP by stage") +
    xlab("Stage") +
    ylab("Metastatic score") +
    scale_y_continuous(breaks=seq(0,100,0.1)) +
    coord_cartesian(ylim=c(1,1.4), expand=FALSE)
dev.off()

# Anova
uvm.metpotential.low.anova <- aov(metscore~stage,  data=uvm.metpotential.low)
summary(uvm.metpotential.low.anova)

uvm.metpotential.high.anova <- aov(metscore~stage,  data=uvm.metpotential.high)
summary(uvm.metpotential.high.anova)

#--------- METSCORE MODEL USING SIGS ---------------------------------------------------------------------------
p <- tp.skcm.metscore.stage

summary(tp.metscore.mod <- lm(p$metscore~p$S7+p$S1+p$S11+p$S6+p$S23))
summary(tp.metscore.mod <- glm(p$met_potential~p$S7+p$S1+p$S11+p$S6+p$S23))
summary(tp.metscore.mod <- lm(p$metscore~p$S7))





#--------- TOTAL MUTATIONS CONTRIBUTION ------------------------------------------------------------------------------------

# Total number of mutations per sample
skcm.total.muts <- as.data.frame(table(snvs$Tumor_Sample_Barcode))
colnames(skcm.total.muts) <- c("tcga_participant_barcode","Total_Muts")
skcm.total.muts[,1] <- sapply(skcm.total.muts$tcga_participant_barcode, function(x) substr(x, 1, 12))

tp.skcm.total.muts <- merge(tp.skcm.metscore.stage, skcm.total.muts, by="tcga_participant_barcode", all=TRUE)
tp.skcm.total.muts <- tp.skcm.total.muts[!is.na(tp.skcm.total.muts$patient_id),c(-2,-3,-37,-40,-41)]

# Replace signature % contribution columns with total mutation contributions (% sig contribution * total muts)
for (i in 2:32) {
    tp.skcm.total.muts[,i] <- tp.skcm.total.muts[,i]*tp.skcm.total.muts[,37]
}

save(tp.skcm.total.muts, file="./Output/TCGA_SKCM_TP_totalmutsigs.RData")

tm.skcm.total.muts <- merge(tm.skcm.metscore.stage, skcm.total.muts, by="tcga_participant_barcode", all=TRUE)
tm.skcm.total.muts <- tm.skcm.total.muts[!is.na(tm.skcm.total.muts$patient_id),c(-2,-3,-37,-40,-41)]

# Replace signature % contribution columns with total mutation contributions (% sig contribution * total muts)
for (i in 2:32) {
    tm.skcm.total.muts[,i] <- tm.skcm.total.muts[,i]*tm.skcm.total.muts[,37]
}

save(tm.skcm.total.muts, file="./Output/TCGA_SKCM_TM_totalmutsigs.RData")

# S1
s1 <- tp.skcm.total.muts[,c(1:2,35,36)]
quantile(s1$S1)
s1.out <- s1[which(s1$S1>quantile(s1$S1)[2][[1]] & s1$S1<quantile(s1$S1)[4][[1]]),]
cor.test(s1$S1, s1$metscore)
test <- cor.test(s1$S1, s1$metscore)
pval <- round(test$p.value, digits=4)
cor <- round(test$estimate[[1]], digits=7)

s1.scat <- ggplot(s1, aes(x=S1, y=metscore)) +
    geom_point(shape=16, size=5, color="#1F78B4B3") + 
    geom_smooth(color="black", method=lm) + 
    ggtitle("S1") +
    xlab("Total mutations contribution") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,200,10)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,137.7), ylim=c(0.6,1.6), expand=FALSE) +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13)) +
    annotate("text", x=110, y=c(1.55, 1.48), label = c(paste0("Correlation: ", cor), paste0("p-value: ", pval)))

s1$met_potential <- factor(s1$met_potential, levels = c("low","high"),ordered = TRUE)

s1.viol <- ggplot(s1, aes(x=met_potential, y=S1, fill=met_potential)) +
    geom_violin(color=NA) +
    scale_fill_manual(values=colors3) +
    ggtitle("") +
    xlab("Metastatic potential") +
    ylab("Total mutations contribution") +
    scale_y_continuous(breaks=seq(0,200,25)) +
    coord_cartesian(ylim=c(0,140), expand=FALSE) +
    theme_bw() +
    stat_summary(fun.y = "median", geom = "point", shape = 18, size = 5, color = "black") +
    scale_x_discrete(labels=c("Low\nn=45", "High\nn=53")) + 
    theme(legend.position="none",
          axis.text = element_text(size=12),
          axis.title = element_text(size=13))

wilcox.test(s1$S1[which(s1$met_potential=="low")], s1$S1[which(s1$met_potential=="high")])


# S7, 1 outlier removed
s7 <- tp.skcm.total.muts[,c(1,8,35,36)]
s7 <- s7[which(s7$S7<4000),]
quantile(s7$S7)
cor.test(s7$S7, s7$metscore)
test <- cor.test(s7$S7, s7$metscore)
pval <- round(test$p.value, digits=4)
cor <- round(test$estimate[[1]], digits=7)

s7.scat <- ggplot(s7, aes(x=S7, y=metscore)) +
    geom_point(shape=16, size=5, color="#1F78B4B3") + 
    geom_smooth(color="black", method=lm) +
    ggtitle("S7") +
    xlab("Total mutations contribution") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,10000,250)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,2575), ylim=c(0.6,1.6), expand=FALSE)+
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13)) +
    annotate("text", x=2100, y=c(1.55, 1.48),  label = c(paste0("Correlation: ", cor), paste0("p-value: ", pval)))

s7$met_potential <- factor(s7$met_potential, levels = c("low","high"),ordered = TRUE)

s7.viol <- ggplot(s7, aes(x=met_potential, y=S7, fill=met_potential)) +
    geom_violin(color=NA) +
    scale_fill_manual(values=colors3) +
    ggtitle("") +
    xlab("Metastatic potential") +
    ylab("Total mutations contribution") +
    scale_y_continuous(breaks=seq(0,100000,500)) +
    coord_cartesian(ylim=c(0,2700), expand=FALSE) +
    theme_bw() +
    stat_summary(fun.y = "median", geom = "point", shape = 18, size = 5, color = "black") +
    scale_x_discrete(labels=c("Low\nn=45", "High\nn=53")) + 
    theme(legend.position="none",
          axis.text = element_text(size=12),
          axis.title = element_text(size=13))

wilcox.test(s7$S7[which(s7$met_potential=="low")], s7$S7[which(s7$met_potential=="high")])

# Output S1 and S7 as one plot
pdf("./Figures/TCGA_SKCM_TP_metscore_pot_totalmuts.pdf", w=8, h=6)
ggarrange(s1.scat, s1.viol, s7.scat, s1.viol, 
          widths = c(2, 1),
          nrow = 2,
          ncol = 2,
          labels = c("a)"))
dev.off()


# S11, 1 outlier removed
s11 <- tp.skcm.total.muts[,c(1,12,35,36)]
s11 <- s11[which(s11$S11<250),]
quantile(s11$S11)
cor.test(s11$S11, s11$metscore)

pdf("./Figures/TCGA_SKCM_TP_metscore_S11_scatter.pdf", w=8, h=6)
ggplot(s11, aes(x=S11, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TP total S11 mutation contribution vs. metastatic score") +
    xlab("Total mutations contributed by S11") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,300,10)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,96), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=83, y=c(1.47, 1.44), label = c("Correlation: -0.06027398", "p-value: 0.5576"))
dev.off()

s11$met_potential <- factor(s11$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TP_metpotential_S11_boxplot.pdf", w=6, h=6)
ggplot(s11, aes(x=met_potential, y=S11)) +
    geom_boxplot() +
    ggtitle("SKCM TP total S11 mutation contribution vs. metastatic potential") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S11") +
    scale_y_continuous(breaks=seq(0,300,10)) +
    coord_cartesian(ylim=c(0,100), expand=FALSE)
dev.off()


# S23, 1 outlier removed
s23 <- tp.skcm.total.muts[,c(1,24,35,36)]
s23 <- s23[which(s23$S23<1400),]
quantile(s23$S23)
cor.test(s23$S23, s23$metscore)

pdf("./Figures/TCGA_SKCM_TP_metscore_S23_scatter.pdf", w=8, h=6)
ggplot(s23, aes(x=S23, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TP total S23 mutation contribution vs. metastatic score") +
    xlab("Total mutations contributed by S23") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,300,10)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,233.5), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=205, y=c(1.47, 1.44), label = c("Correlation: 0.08380639", "p-value: 0.4144"))
dev.off()

s23$met_potential <- factor(s23$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TP_metpotential_S23_boxplot.pdf", w=6, h=6)
ggplot(s23, aes(x=met_potential, y=S23)) +
    geom_boxplot() +
    ggtitle("SKCM TP total S23 mutation contribution vs. metastatic potential") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S23") +
    scale_y_continuous(breaks=seq(0,300,10)) +
    coord_cartesian(ylim=c(0,235), expand=FALSE)
dev.off()


# S6
s6 <- tp.skcm.total.muts[,c(1,7,35,36)]
quantile(s6$S6)
cor.test(s6$S6, s6$metscore)

pdf("./Figures/TCGA_SKCM_TP_metscore_S6_scatter.pdf", w=8, h=6)
ggplot(s6, aes(x=S6, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TP total S6 mutation contribution vs. metastatic score") +
    xlab("Total mutations contributed by S6") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,300,10)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,72), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=62, y=c(1.47, 1.44), label = c("Correlation: -0.07751905", "p-value: 0.448"))
dev.off()

s6$met_potential <- factor(s6$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TP_metpotential_S6_boxplot.pdf", w=6, h=6)
ggplot(s6, aes(x=met_potential, y=S6)) +
    geom_boxplot() +
    ggtitle("SKCM TP total S6 mutation contribution vs. metastatic potential") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S6") +
    scale_y_continuous(breaks=seq(0,300,10)) +
    coord_cartesian(ylim=c(0,75), expand=FALSE)
dev.off()


# S7/S1 1 outlier removed
s7.1 <- tp.skcm.total.muts[,c(1:2,8,35,36)]
s7.1$S7.1 <- s7.1$S7/s7.1$S1
s7.1 <- s7.1[which(s7.1$S7.1<500),]
s7.1 <- s7.1[complete.cases(s7.1),] # remove NaN (1 sample)
s7.1 <- s7.1[which(s7.1$S7.1 != Inf),] # remove Inf (11 samples)

quantile(s7.1$S7.1, na.rm=TRUE)
cor.test(s7.1$S7.1, s7.1$metscore)

pdf("./Figures/TCGA_SKCM_TP_metscore_S7.1_scatter.pdf", w=8, h=6)
ggplot(s7.1, aes(x=S7.1, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TP total mutation contribution S7/S1 vs. metastatic score") +
    xlab("Total mutations contributed by S7/S1") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,1000,10)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,205.3), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=180, y=c(1.47, 1.44), label = c("Correlation: 0.05790048 ", "p-value: 0.5986"))
dev.off()

s7.1$met_potential <- factor(s7.1$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TP_metpotential_S7.1_boxplot.pdf", w=6, h=6)
ggplot(s7.1, aes(x=met_potential, y=S7.1)) +
    geom_boxplot() +
    ggtitle("SKCM TP total mutation contribution S7/S1 vs. metastatic potential") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S7/S1") +
    scale_y_continuous(breaks=seq(0,200,10)) +
    coord_cartesian(ylim=c(0,206), expand=FALSE)
dev.off()



#--------- MATCHED NORMAL SAMPLES ------------------------------------------------------------------------------------------

# Blood derived normal (NB), Solid tissue normal (NT), Buccal cell normal (NBC),
# EBV immortalized normal (NEBV), Bone marrow normal (NBM)
norm.skcm <- Samples.mRNASeq(format = "csv",
                                 gene = order.updown,
                                 cohort = "SKCM",
                                 sample_type = c("NT"),
                                 protocol = "RSEM",
                                 page_size=6000
)
# 68 rows of NT samples, all with same participant barcode: TCGA-GN-A4U8 (matches a TM sample)



#--------- METASTATIC SITES ----------------------------------------------------------------------------------------------

# Generate df of mutsig, metscore, metsite, and stage, SKCM TM
tm.metsite <- merge(tm.clin[,c(2,23)], tm.skcm.metscore.stage[,c(1,4:36,38:39)], by="tcga_participant_barcode")

# Group metastatic sites
tm.metsite$metsite_group <- tm.metsite$distant_metastasis_anatomic_site

for (i in 1:nrow(tm.metsite)) {
    if (is.na(tm.metsite$metsite_group[i])) {
        tm.metsite$metsite_group[i] <- NA
    } else if (grepl("lymph", tm.metsite$metsite_group[i], ignore.case = TRUE)) {
        tm.metsite$metsite_group[i] <- "lymph node"
    } else if (grepl("skin", tm.metsite$metsite_group[i], ignore.case = TRUE)) {
        tm.metsite$metsite_group[i] <- "skin"
    } else if (grepl("brain", tm.metsite$metsite_group[i], ignore.case = TRUE)) {
        tm.metsite$metsite_group[i] <- "brain"
    } else if (grepl("lung", tm.metsite$metsite_group[i], ignore.case = TRUE)) {
        tm.metsite$metsite_group[i] <- "lung"
    } else if (grepl("soft tissue", tm.metsite$metsite_group[i], ignore.case = TRUE)) {
        tm.metsite$metsite_group[i] <- "soft tissue"
    } else if (grepl("bone", tm.metsite$metsite_group[i], ignore.case = TRUE) | grepl("skull", tm.metsite$metsite_group[i], ignore.case = TRUE)) {
        tm.metsite$metsite_group[i] <- "bone"
    } else if (grepl("abdomen", tm.metsite$metsite_group[i], ignore.case = TRUE) | grepl("colon", tm.metsite$metsite_group[i], ignore.case = TRUE) | 
        grepl("duodenum", tm.metsite$metsite_group[i], ignore.case = TRUE) | grepl("jejunum", tm.metsite$metsite_group[i], ignore.case = TRUE) |
        grepl("intestine", tm.metsite$metsite_group[i], ignore.case = TRUE) | grepl("bowel", tm.metsite$metsite_group[i], ignore.case = TRUE) ) {
        tm.metsite$metsite_group[i] <- "GI tract"
    }
    else {
        tm.metsite$metsite_group[i] <- "other"
    }
}

# ---- ~ Metscore ----

cols8 <- c(brewer.pal(8, "Paired"), "#000000")

ggplot(tm.metsite, aes(x=metsite_group, y=metscore, fill=metsite_group)) +
    geom_boxplot() +
    scale_fill_manual(values=cols8) +
    ggtitle("") +
    xlab("Metastatic site") +
    ylab("Metastatic score") +
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size=12),
          axis.title = element_text(size=13)) +
    annotate("text", x=1, y=1.8, label=paste0("p=",pval))

# Anova
an <- summary(metsite.anova <- aov(metscore~metsite_group,  data=tm.metsite))
pval <- round(an[[1]][["Pr(>F)"]][[1]], digits=4)
pairwise <- as.data.frame(TukeyHSD(metsite.anova)[[1]])
pairwise$sign <- NA
colnames(pairwise)[4] <- "padj"
pairwise$sign <- ifelse(pairwise$padj < 0.05, "*","")

write.csv(pairwise, file="./Output/TCGA_SKCM_TM_metsite_anova_pairwise.csv")

# compare lymph node to everything else
tm.metsite$lymph_met <- ifelse(tm.metsite$metsite_group=="lymph node", "lymph node", "other")
tm.metsite$brain_met <- ifelse(tm.metsite$metsite_group=="brain", "brain", "other")

pval_lymph <- round(wilcox.test(tm.metsite$metscore[which(tm.metsite$lymph_met=="lymph node")], tm.metsite$metscore[which(tm.metsite$lymph_met=="other")])$p.value, digits=4)

ggplot(tm.metsite, aes(x=lymph_met, y=metscore, fill=lymph_met)) +
    geom_boxplot() +
    scale_fill_manual(values=cols8[c(5,6)]) +
    ggtitle("") +
    xlab("Metastatic site") +
    ylab("Metastatic score") +
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size=12),
          axis.title = element_text(size=13)) +
    annotate("text", x=1, y=1.8, label=paste0("p=",pval_lymph))


pval_brain <- round(wilcox.test(tm.metsite$metscore[which(tm.metsite$brain_met=="brain")], tm.metsite$metscore[which(tm.metsite$brain_met=="other")])$p.value, digits=4)

ggplot(tm.metsite, aes(x=brain_met, y=metscore, fill=brain_met)) +
    geom_boxplot() +
    scale_fill_manual(values=cols8[c(6,5)]) +
    ggtitle("") +
    xlab("Metastatic site") +
    ylab("Metastatic score") +
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size=12),
          axis.title = element_text(size=13)) +
    annotate("text", x=1, y=1.8, label=paste0("p=",pval_brain))


# ---- ~ Metpotential ----

sites <- c("bone","brain", "GI", "lung", "lymph", "other", "skin", "soft")
metsite_pot <- data.frame(data=sites, low=NA, high=NA)

for (i in 1:length(sites)) {
    metsite_pot[i,2] <- nrow(tm.metsite[which(tm.metsite$metsite_group == sites[i] & tm.metsite$met_potential == "low"),])
    metsite_pot[i,3] <- nrow(tm.metsite[which(tm.metsite$metsite_group == sites[i] & tm.metsite$met_potential == "high"),])
}

# Generates contingency tables and runs Fisher's exact test --> outputs into data frame
fisher <- data.frame(metsite=metsite_pot$data, OR=NA, pval=NA, padj=NA)
for (i in 1:nrow(metsite_pot)) {
    p <- metsite_pot[i,2]
    m <- metsite_pot[i,3]
    sump <- colSums(metsite_pot[2])-p
    summ <- colSums(metsite_pot[3])-m
    mat <- array(c(p, sump, m, summ), c(2,2))
    fisher[i,2] <- fisher.test(mat)$estimate
    fisher[i,3] <- fisher.test(mat)$p.value
}
fisher$padj <- p.adjust(fisher$pval, method="BH")  
write.csv(fisher, file="./Output/TCGA_SKCM_fisher.csv")   


#--------- S1 VS. AGE ------------------------------------------------------------------------------------------

# Correlation test per cohort
cor.test(tp.skcm.totalmuts$S1, tp.clin$age_at_initial_pathologic_diagnosis)
cor.test(tm.skcm.totalmuts$S1, tm.clin$age_at_initial_pathologic_diagnosis)
cor.test(uvm.totalmuts$S1, clin$age_at_initial_pathologic_diagnosis)

# Plot SKCM TP
ggplot(tp.clin, aes(x=S1, y=age_at_initial_pathologic_diagnosis)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) +
    xlab("S1") +
    ylab("Age") +
    scale_x_continuous(breaks=seq(0,1000,10))+
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,205.3), ylim=c(0.6,1.5), expand=FALSE) +
    annotate("text", x=180, y=c(1.47, 1.44), label = c("Correlation: 0.05790048 ", "p-value: 0.5986"))

ggplot(tp.clin, aes(x=age_at_initial_pathologic_diagnosis, y=S1)) +
    geom_point(shape=16, size=5, color="#1F78B4B3") + 
    geom_smooth(color="black", method=lm) + 
    ggtitle("SKCM TP") +
    xlab("Age") +
    ylab("S1") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13))

# No significant correlations between S1 and age
