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
load("./Output/TCGA_SKCM_TP_expression_matrix.RData")  # SKCM TP
load("./Output/TCGA_SKCM_TM_expression_matrix.RData")  # SKCM TM
load("./Output/TCGA_UVM_TP_expression_matrix.RData")   # UVM TP

# Metastatic score
load("./Output/TCGA_SKCM_TP_metastatic_score.RData")  # SKCM TP
load("./Output/TCGA_SKCM_TM_metastatic_score.RData")  # SKCM TM
load("./Output/TCGA_UVM_TP_metastatic_score.RData")   # UVM TP


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
# cBioPortal genes
genes <- c("MMP1","MMP2","MMP3","MMP7","MMP9","MMP10","MMP11","MMP12","MMP13","MMP14","MMP15","MMP16","MMP17",
           "MMP19","MMP21","MMP23B","MMP24","MMP25","MMP26","MMP27","MMP28","ITGB3","ITGAV","PTK2","CDH1","SPARC",
           "WFDC2")

tp.skcm.expr <- Samples.mRNASeq(format = "csv",
                           gene = genes,
                           cohort = "SKCM",
                           sample_type = "TP",
                           protocol = "RSEM",
                           page_size=6000
)
save(tp.skcm.expr, file="./Data/TCGA_SKCM_TP_expression.RData")

tm.skcm.expr <- Samples.mRNASeq(format = "csv",
                                gene = genes,
                                cohort = "SKCM",
                                sample_type = "TM",
                                protocol = "RSEM",
                                page_size=6000
)
save(tm.skcm.expr, file="./Data/TCGA_SKCM_TM_expression.RData")

uvm.expr <- Samples.mRNASeq(format = "csv",
                                gene = genes,
                                cohort = "UVM",
                                sample_type = "TP",
                                protocol = "RSEM",
                                page_size=6000
)
save(uvm.expr, file="./Data/TCGA_UVM_TP_expression.RData")

# Generate matrix of participant and gene expr
tp.skcm.expr$expression_log2 <- as.numeric(tp.skcm.expr$expression_log2)
mat.tp.skcm.expr <- cast(tcga_participant_barcode~gene,data=tp.skcm.expr,
                value="expression_log2", fun.aggregate=mean)
rownames(mat.tp.skcm.expr) <- mat.tp.skcm.expr$tcga_participant_barcode
mat.tp.skcm.expr <- mat.tp.skcm.expr[,-1]

tm.skcm.expr$expression_log2 <- as.numeric(tm.skcm.expr$expression_log2)
mat.tm.skcm.expr <- cast(tcga_participant_barcode~gene,data=tm.skcm.expr,
                         value="expression_log2", fun.aggregate=mean)
rownames(mat.tm.skcm.expr) <- mat.tm.skcm.expr$tcga_participant_barcode
mat.tm.skcm.expr <- mat.tm.skcm.expr[,-1]

uvm.expr$expression_log2 <- as.numeric(uvm.expr$expression_log2)
mat.uvm.expr <- cast(tcga_participant_barcode~gene,data=uvm.expr,
                         value="expression_log2", fun.aggregate=mean)
rownames(mat.uvm.expr) <- mat.uvm.expr$tcga_participant_barcode
mat.uvm.expr <- mat.uvm.expr[,-1]


save(mat.tp.skcm.expr, file="./Output/TCGA_SKCM_TP_expression_matrix.RData")
save(mat.tm.skcm.expr, file="./Output/TCGA_SKCM_TM_expression_matrix.RData")
save(mat.uvm.expr, file="./Output/TCGA_UVM_TP_expression_matrix.RData")

hist(mat.uvm.expr$met_score)

# Heatmap overview of each cohort
tp.skcm <- as.data.frame(mat.tp.skcm.expr[,-28])
tp.skcm[is.na(tp.skcm)] <- 0                        # Replace NA with 0
pdf("./Figures/TCGA_SKCM_TP_expression_heatmap.pdf", w=10, h=6, onefile=FALSE)
tp.skcm.map <-pheatmap(t(tp.skcm), show_colnames = FALSE, main="Metastatis genes expression SKCM TP")
dev.off()

tm.skcm <- as.data.frame(mat.tm.skcm.expr[,-28])
tm.skcm[is.na(tm.skcm)] <- 0                        # Replace NA with 0
pdf("./Figures/TCGA_SKCM_TM_expression_heatmap.pdf", w=10, h=6, onefile=FALSE)
tm.skcm.map <-pheatmap(t(tm.skcm), show_colnames = FALSE, main="Metastatis genes expression SKCM TM")
dev.off()

uvm <- as.data.frame(mat.uvm.expr[,-28])
uvm[is.na(uvm)] <- 0                        # Replace NA with 0
pdf("./Figures/TCGA_UVM_TP_expression_heatmap.pdf", w=10, h=6, onefile=FALSE)
uvm.map <-pheatmap(t(uvm), show_colnames = FALSE, main="Metastatis genes expression UVM TP")
dev.off()


#--------- METASTATIC SCORE ----------------------------------------------------------------------------------------

# Metastatic potential score per sample = average expression of genes
tp.skcm.metscore <- data.frame(tcga_participant_barcode=noquote(rownames(mat.tp.skcm.expr)),
                       met_score=rowMeans(mat.tp.skcm.expr, na.rm=TRUE))
tm.skcm.metscore <- data.frame(tcga_participant_barcode=noquote(rownames(mat.tm.skcm.expr)),
                               met_score=rowMeans(mat.tm.skcm.expr, na.rm=TRUE))
uvm.metscore <- data.frame(tcga_participant_barcode=noquote(rownames(mat.uvm.expr)),
                               met_score=rowMeans(mat.uvm.expr, na.rm=TRUE))

save(tp.skcm.metscore, file="./Output/TCGA_SKCM_TP_metastatic_score.RData")
save(tm.skcm.metscore, file="./Output/TCGA_SKCM_TM_metastatic_score.RData")
save(uvm.metscore, file="./Output/TCGA_UVM_TP_metastatic_score.RData")


#--------- SKCM TP, TM, VS. UVM ----------------------------------------------------------------------------------------

tp.skcm.metscore$cohort <- "SKCM TP"
tm.skcm.metscore$cohort <- "SKCM TM"
uvm.metscore$cohort <- "UVM TP"

metscore <- rbind(tp.skcm.metscore, tm.skcm.metscore, uvm.metscore)
order2 <- factor(metscore$cohort, levels=c("SKCM TP", "SCKM TM", "UVM"))

pdf("./Figures/TCGA_SKCM_UVM_metscore_boxplot.pdf", w=8, h=6)
ggplot(metscore, aes(x=cohort, y=met_score)) + 
    geom_boxplot() +
    scale_fill_manual(values=colors4) +
    ggtitle("Metastatic score SKCM TP, SKCM TM, and UVM TP") +
    xlab("Cohort") +
    ylab("Metastatic potential score") +
    scale_y_continuous(breaks=seq(0,10,0.5)) +
    coord_cartesian( ylim=c(0,10), expand = FALSE )
dev.off()

metscore.anova <- aov(met_score~cohort,  data=metscore)
summary(metscore.anova)

compare <- data.frame(c1=c("SKCM TP","SKCM TP","SKCM TM"), c2=c("SKCM TM","UVM","UVM"), pval=NA)
compare[1,3] <- t.test(tp.skcm.metscore$met_score, tm.skcm.metscore$met_score)$p.value
compare[2,3] <- t.test(tp.skcm.metscore$met_score, uvm.metscore$met_score)$p.value
compare[3,3] <- t.test(tm.skcm.metscore$met_score, uvm.metscore$met_score)$p.value

compare$padj <- p.adjust(compare$pval, method="BH")


#--------- CLINICAL DATA ----------------------------------------------------------------------------------------

#--------- ~ Stage ---------
tp.skcm.metscore.stage <- merge(tp.stage, tp.skcm.metscore, by="tcga_participant_barcode")
tm.skcm.metscore.stage <- merge(tm.stage, tm.skcm.metscore, by="tcga_participant_barcode")
uvm.metscore.stage <- merge(stage, uvm.metscore, by="tcga_participant_barcode")

tp.skcm.metscore.stage <- tp.skcm.metscore.stage[!is.na(tp.skcm.metscore.stage$stage),]     # 4 NA removed
tm.skcm.metscore.stage <- tm.skcm.metscore.stage[!is.na(tm.skcm.metscore.stage$stage),]     # 34 NA removed
uvm.metscore.stage <- uvm.metscore.stage[!is.na(uvm.metscore.stage$stage),]                 # 0 NA removed

# Boxplot by stage
pdf("./Figures/TCGA_SKCM_TP_metscore_stage_boxplot.pdf", w=10, h=6)
ggplot(tp.skcm.metscore.stage, aes(x=stage, y=met_score, fill=stage)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("Metastatic score SKCM TP by cancer stage") +
    xlab("Cancer stage") +
    ylab("Metastatic potential score (log2(expression))") +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_cartesian( ylim=c(0,10), expand = FALSE )
dev.off()

pdf("./Figures/TCGA_SKCM_TM_metscore_stage_boxplot.pdf", w=10, h=6)
ggplot(tm.skcm.metscore.stage, aes(x=stage, y=met_score, fill=stage)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("Metastatic score SKCM TM by cancer stage") +
    xlab("Cancer stage") +
    ylab("Metastatic potential score (log2(expression))") +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_cartesian( ylim=c(0,10), expand = FALSE )
dev.off()

pdf("./Figures/TCGA_UVM_TP_metscore_stage_boxplot.pdf", w=8, h=6)
ggplot(uvm.metscore.stage, aes(x=stage, y=met_score, fill=stage)) + 
    geom_boxplot() +
    scale_fill_manual(values=colors7) +
    ggtitle("Metastatic score UVM TP by cancer stage") +
    xlab("Cancer stage") +
    ylab("Metastatic potential score (log2(expression))") +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_cartesian( ylim=c(0,10), expand = FALSE )
dev.off()

# Anova
tp.skcm.metscore.stage.anova <- aov(met_score~stage,  data=tp.skcm.metscore.stage)
summary(tp.skcm.metscore.stage.anova)

tm.skcm.metscore.stage.anova <- aov(met_score~stage,  data=tm.skcm.metscore.stage)
summary(tm.skcm.metscore.stage.anova)

uvm.metscore.stage.anova <- aov(met_score~stage,  data=uvm.metscore.stage)
summary(uvm.metscore.stage.anova)


# Boxplot by disease progression
pdf("./Figures/TCGA_SKCM_TP_metscore_prog_boxplot.pdf", w=5, h=6)
ggplot(tp.skcm.metscore.stage, aes(x=earlyLate, y=met_score, fill=earlyLate)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("Metastatic score SKCM TP by disease progression") +
    xlab("Disease progression") +
    ylab("Metastatic potential score (log2(expression))") +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_cartesian( ylim=c(0,10), expand = FALSE )
dev.off()

pdf("./Figures/TCGA_SKCM_TM_metscore_prog_boxplot.pdf", w=5, h=6)
ggplot(tm.skcm.metscore.stage, aes(x=earlyLate, y=met_score, fill=earlyLate)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("Metastatic score SKCM TM by disease progression") +
    xlab("Disease progression") +
    ylab("Metastatic potential score (log2(expression))") +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_cartesian( ylim=c(0,10), expand = FALSE )
dev.off()

pdf("./Figures/TCGA_UVM_TP_metscore_prog_boxplot.pdf", w=5, h=6)
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
t.test(tp.skcm.metscore.stage.early$met_score, tp.skcm.metscore.stage.late$met_score)

tm.skcm.metscore.stage.early <- tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate=="early"),]
tm.skcm.metscore.stage.late <- tm.skcm.metscore.stage[which(tm.skcm.metscore.stage$earlyLate=="late"),]
t.test(tm.skcm.metscore.stage.early$met_score, tm.skcm.metscore.stage.late$met_score)

uvm.metscore.stage.early <- uvm.metscore.stage[which(uvm.metscore.stage$earlyLate=="early"),]
uvm.metscore.stage.late <- uvm.metscore.stage[which(uvm.metscore.stage$earlyLate=="late"),]
t.test(uvm.metscore.stage.early$met_score, uvm.metscore.stage.late$met_score)