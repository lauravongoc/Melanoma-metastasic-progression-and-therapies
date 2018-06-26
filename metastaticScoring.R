# Metastatic Scoring
# DT Laura Vo Ngoc
# Start: 28/05/2018

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
load("./Output/TCGA_UVM_TP_metastatic_score_mel.RData")     # UVM TP


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
genes_uvm <- c("HTR2B","ECM1","RAB31","CDH1","FXR1","LTA4H","EIF1B","ID2","ROBO1","LMCD1","SATB1","MTUS1","PRAME")

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
uvm.genelist.up <- c("HTR2B", "ECM1", "RAB31", "CDH1")
uvm.genelist.down <- c("FXR1","LTA4H","EIF1B","ID2","ROBO1","LMCD1","SATB1","MTUS1","PRAME")


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
                           metscore=(rowMeans(uvm.genelist.expr[,1:4]/rowMeans(tm.skcm.genelist.expr[,5:13]),  
                                                   na.rm=TRUE)))

save(tp.skcm.metscore, file="./Output/TCGA_SKCM_TP_metastatic_score_mel.RData")
save(tm.skcm.metscore, file="./Output/TCGA_SKCM_TM_metastatic_score_mel.RData")
save(uvm.metscore, file="./Output/TCGA_UVM_TP_metastatic_score_mel.RData")

# Overlapping histogram of metscore distributions
pdf("./Figures/TCGA_SKCM_UVM_metscore_hist.pdf", w=8, h=6)
hist(tm.skcm.metscore$metscore, 
     col="#1F78B4B3", 
     main="Metastatic score distribution", 
     xlab="Metastatic score",
     xlim=c(0.5,1.8),
     ylim=c(0,100))
hist(tp.skcm.metscore$metscore, 
     col="#fb9a99B3", 
     main="Metastatic score distribution SKCM TP", 
     xlab="Metastatic score",
     xlim=c(0.5,1.8),
     ylim=c(0,100),
     add=TRUE)
hist(uvm.metscore$metscore,
     col="#B2DF8AB3", 
     main="Metastatic score distribution UVM TP", 
     xlab="Metastatic score",
     xlim=c(0.5,1.8),
     ylim=c(0,100),
     add=TRUE)
legend("topright", c("SKCM TP", "SKCM TM", "UVM TP"), col=c("#fb9a99B3","#1F78B4B3", "#B2DF8AB3"), lwd=10)
dev.off()


#--------- ~ Cutoff ---------

# Cutoff = histogram intersects ~= 1.0

tp.skcm.metscore$met_potential <- NA
tp.skcm.metscore$met_potential <- ifelse(tp.skcm.metscore$metscore > 1, "high", "low")

tm.skcm.metscore$met_potential <- NA
tm.skcm.metscore$met_potential <- ifelse(tm.skcm.metscore$metscore > 1, "high", "low")

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

# Create table of low vs high potential counts with early or late disease progression
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

# Fisher's exact test
fisher.test(earlyLate_hilow)

# Create table of low vs high potential counts with early or late disease progression
earlyLate_hiLow2 <- data.frame(data=c("early","late"), low=NA, high=NA)
rownames(earlyLate_hiLow2) <- earlyLate_hiLow2[,1] 
earlyLate_hiLow2 <- earlyLate_hiLow2[,-1]
earlyLate_hiLow2[1,1] <- nrow(tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$earlyLate == "early" & tp.skcm.metscore.stage$met_potential_quant == "low"),]) + 
    nrow(tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate == "early" & tm.skcm.metscore.stage$met_potential_quant == "low"),])
earlyLate_hiLow2[1,2] <- nrow(tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$earlyLate == "early" & tp.skcm.metscore.stage$met_potential_quant == "high"),]) + 
    nrow(tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate == "early" & tm.skcm.metscore.stage$met_potential_quant == "high"),])
earlyLate_hiLow2[2,1] <- nrow(tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$earlyLate == "late" & tp.skcm.metscore.stage$met_potential_quant == "low"),]) + 
    nrow(tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate == "late" & tm.skcm.metscore.stage$met_potential_quant == "low"),])
earlyLate_hiLow2[2,2] <- nrow(tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$earlyLate == "late" & tp.skcm.metscore.stage$met_potential_quant == "high"),]) + 
    nrow(tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate == "late" & tm.skcm.metscore.stage$met_potential_quant == "high"),])

# Fisher's exact test
fisher.test(earlyLate_hiLow2)

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

#--------- ~ Stage ---------
tp.skcm.metscore.stage <- merge(tp.stage, tp.skcm.metscore, by="tcga_participant_barcode")
tm.skcm.metscore.stage <- merge(tm.stage, tm.skcm.metscore, by="tcga_participant_barcode")
#uvm.metscore.stage <- merge(stage, uvm.metscore, by="tcga_participant_barcode")

save(tp.skcm.metscore.stage, file="./Output/TCGA_SKCM_TP_metastatic_score_stage.RData")
save(tm.skcm.metscore.stage, file="./Output/TCGA_SKCM_TM_metastatic_score_stage.RData")

# remove stage NA
tp.skcm.metscore.stage <- tp.skcm.metscore.stage[!is.na(tp.skcm.metscore.stage$stage),]            # 4 NA removed
tm.skcm.metscore.stage <- tm.skcm.metscore.stage[!is.na(tm.skcm.metscore.stage$stage),]            # 34 NA removed
#uvm.metscore.stage <- uvm.metscore.stage[!is.na(uvm.metscore.stage$stage),]                 # 0 NA removed


# Boxplot by stage
pdf("./Figures/TCGA_SKCM_TP_metscore_mel_stage_boxplot.pdf", w=10, h=6)
ggplot(tp.skcm.metscore.stage, aes(x=stage, y=metscore, fill=stage)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("Metastatic score SKCM TP by cancer stage") +
    xlab("Cancer stage") +
    ylab("Metastatic potential score") +
    scale_y_continuous(breaks=seq(0,10,0.1)) +
    coord_cartesian( ylim=c(0,2), expand = FALSE )
dev.off()

pdf("./Figures/TCGA_SKCM_TM_metscore_mel_stage_boxplot.pdf", w=10, h=6)
ggplot(tm.skcm.metscore.stage, aes(x=stage, y=metscore, fill=stage)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("Metastatic score SKCM TM by cancer stage") +
    xlab("Cancer stage") +
    ylab("Metastatic potential score") +
    scale_y_continuous(breaks=seq(0,10,0.1)) +
    coord_cartesian( ylim=c(0,2), expand = FALSE )
dev.off()

#pdf("./Figures/TCGA_UVM_TP_metscore_mel_stage_boxplot.pdf", w=8, h=6)
ggplot(uvm.metscore.stage, aes(x=stage, y=met_score, fill=stage)) + 
    geom_boxplot() +
    scale_fill_manual(values=colors7) +
    ggtitle("Metastatic score UVM TP by cancer stage") +
    xlab("Cancer stage") +
    ylab("Metastatic potential score (log2(expression))") +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_cartesian( ylim=c(5,10), expand = FALSE )
#dev.off()

# Anova
tp.skcm.metscore.stage.anova <- aov(metscore~stage,  data=tp.skcm.metscore.stage)
summary(tp.skcm.metscore.stage.anova)

tm.skcm.metscore.stage.anova <- aov(metscore~stage,  data=tm.skcm.metscore.stage)
summary(tm.skcm.metscore.stage.anova)

#uvm.metscore.stage.anova <- aov(met_score~stage,  data=uvm.metscore.stage)
#summary(uvm.metscore.stage.anova)


# Boxplot by disease progression
pdf("./Figures/TCGA_SKCM_TP_metscore_mel_prog_boxplot.pdf", w=5, h=6)
ggplot(tp.skcm.metscore.stage, aes(x=earlyLate, y=metscore, fill=earlyLate)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("Metastatic score SKCM TP by disease progression") +
    xlab("Disease progression") +
    ylab("Metastatic potential score") +
    scale_y_continuous(breaks=seq(0,10,0.1)) +
    coord_cartesian( ylim=c(0,2), expand = FALSE )
dev.off()

pdf("./Figures/TCGA_SKCM_TM_metscore_mel_prog_boxplot.pdf", w=5, h=6)
ggplot(tm.skcm.metscore.stage, aes(x=earlyLate, y=metscore, fill=earlyLate)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("Metastatic score SKCM TM by disease progression") +
    xlab("Disease progression") +
    ylab("Metastatic potential score") +
    scale_y_continuous(breaks=seq(0,10,0.1)) +
    coord_cartesian( ylim=c(0,2), expand = FALSE )
dev.off()

pdf("./Figures/TCGA_UVM_TP_metscore_mel_prog_boxplot.pdf", w=5, h=6)
ggplot(uvm.metscore.stage, aes(x=earlyLate, y=met_score, fill=earlyLate)) + 
    geom_boxplot() +
    scale_fill_manual(values=colors4) +
    ggtitle("Metastatic score UVM TP by disease progression") +
    xlab("Disease progression") +
    ylab("Metastatic potential score (log2(expression))") +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_cartesian( ylim=c(0,10), expand = FALSE )
dev.off()

# T-test early vs late
tp.skcm.metscore.stage.early <- tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$earlyLate=="early"),]
tp.skcm.metscore.stage.late <- tp.skcm.metscore.stage[which(tp.skcm.metscore.stage$earlyLate=="late"),]
t.test(tp.skcm.metscore.stage.early$metscore, tp.skcm.metscore.stage.late$metscore)

tm.skcm.metscore.stage.early <- tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate=="early"),]
tm.skcm.metscore.stage.late <- tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate=="late"),]
t.test(tm.skcm.metscore.stage.early$metscore, tm.skcm.metscore.stage.late$metscore)

uvm.metscore.stage.early <- uvm.metscore.stage[which(uvm.metscore.stage$earlyLate=="early"),]
uvm.metscore.stage.late <- uvm.metscore.stage[which(uvm.metscore.stage$earlyLate=="late"),]
t.test(uvm.metscore.stage.early$met_score, uvm.metscore.stage.late$met_score)



#--------- ~ Met potential stage ---------

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


