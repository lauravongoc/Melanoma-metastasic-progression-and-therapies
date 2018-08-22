# Scwann cell markers
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

# Mutsigs
load("./Output/TCGA_SKCM_TP_weights_cut0.00.RData")     #SKCM TP
load("./Output/TCGA_SKCM_TM_weights_cut0.00.RData")     #SKCM TM
load("./Output/TCGA_UVM_TP_weights_cut0.00.RData")      #UVM

# Metastatic score
load("./Output/TCGA_SKCM_TP_metastatic_score_mel.RData")  # SKCM TP
load("./Output/TCGA_SKCM_TM_metastatic_score_mel.RData")  # SKCM TM
load("./Output/TCGA_UVM_TP_metastatic_score_mel.RData")   # UVM TP

# Survival
load("./Output/TCGA_SKCM_TP_survival.RData")        # SKCM TP
load("./Output/TCGA_SKCM_TM_survival.RData")        # SKCM TM
load("./Output/TCGA_UVM_TP_survival.RData")         # UVM

# Genelist
load("./Output/Schwann_genelist.RData")

# Expression data Schwann markers
load("./Data/TCGA_SKCM_TP_Schwann_expression.RData")    # SKCM TP 
load("./Data/TCGA_SKCM_TM_Schwann_expression.RData")    # SKCM TM
load("./Data/TCGA_UVM_Schwann_expression.RData")        # UVM

# Total mutsig mutations
load("./Output/TCGA_SKCM_TP_totalmutsigs.RData")    # SKCM TP
load("./Output/TCGA_SKCM_TM_totalmutsigs.RData")    # SKCM TM



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

colors.trans <- c("#1F78B4B3", "#fb9a99B3", "#B2DF8AB3")


#--------- EXPRESSION DATA -------------------------------------------------------------------------------------------------
schwann.genes <- c("GFAP", "SOX10", "GAP43", "P75NTR")                  # No expression data P75NTR
save(schwann.genes, file="./Output/Schwann_genelist.RData")

# SKCM TP - 103 patients
tp.skcm.expr.schwann <- Samples.mRNASeq(format = "csv",
                                     gene = schwann.genes,
                                     cohort = "SKCM",
                                     sample_type = "TP",
                                     protocol = "RSEM",
                                     page_size=6000)
tp.skcm.expr.schwann$expression_log2 <- as.numeric(tp.skcm.expr.schwann$expression_log2)
save(tp.skcm.expr.schwann, file="./Data/TCGA_SKCM_TP_Schwann_expr.RData")


# SKCM TM -- 368 patients
tm.skcm.expr.schwann <- Samples.mRNASeq(format = "csv",
                                        gene = schwann.genes,
                                        cohort = "SKCM",
                                        sample_type = "TM",
                                        protocol = "RSEM",
                                        page_size=6000)
tm.skcm.expr.schwann$expression_log2 <- as.numeric(tm.skcm.expr.schwann$expression_log2)
save(tm.skcm.expr.schwann, file="./Data/TCGA_SKCM_TM_Schwann_expr.RData")

# UVM TP -- 80 patients
uvm.expr.schwann <- Samples.mRNASeq(format = "csv",
                                        gene = schwann.genes,
                                        cohort = "UVM",
                                        sample_type = "TP",
                                        protocol = "RSEM",
                                        page_size=6000)
uvm.expr.schwann$expression_log2 <- as.numeric(uvm.expr.schwann$expression_log2)
save(uvm.expr.schwann, file="./Data/TCGA_UVM_TP_Schwann_expr.RData")



# ---- ~ Distribution ----

# Overlapping histogram of metscore distributions
pdf("./Figures/TCGA_SKCM_UVM_schwann_expr_hist.pdf", w=8, h=6)
hist(tm.skcm.expr.schwann$expression_log2, 
     col="#1F78B4B3",  
     main="Schwann marker expression distribution", 
     xlab="Schwann marker expression",
     xlim=c(-2,16),
     ylim=c(0,270))
hist(tp.skcm.expr.schwann$expression_log2, 
     col="#fb9a99B3", 
     main="Schwann marker expression distribution SKCM TP", 
     xlab="Schwann marker expression",
     xlim=c(-2,16),
     ylim=c(0,270),
     add=TRUE)
hist(uvm.expr.schwann$expression_log2,
     col="#B2DF8AB3", 
     main="Schwann marker expression distribution UVM TP", 
     xlab="Schwann marker expression",
     #breaks=14,
     xlim=c(-2,16),
     ylim=c(0,270),
     add=TRUE)
legend("topright", c("SKCM TP", "SKCM TM", "UVM TP"), col=c("#fb9a99B3","#1F78B4B3", "#B2DF8AB3"), lwd=10, bty="n")
dev.off()


#--------- METSCORE -------------------------------------------------------------------------------------------------

# ---- . SKCM TP ---- 
#168 outliers removed
tp.skcm.schwann.metscore <- merge(tp.skcm.expr.schwann, tp.skcm.metscore[c(1,3)], by="tcga_participant_barcode")
tp.q <- quantile(tp.skcm.schwann.metscore$expression_log2, na.rm=T)
tp.skcm.schwann.metscore <- tp.skcm.schwann.metscore[which(tp.skcm.schwann.metscore$expression_log2>tp.q[2] & tp.skcm.schwann.metscore$expression_log2<tp.q[4]),]
cor.test(tp.skcm.schwann.metscore$expression_log2, tp.skcm.schwann.metscore$metscore)
tp.pval <- round(cor.test(tp.skcm.schwann.metscore$expression_log2, tp.skcm.schwann.metscore$metscore)[3][[1]], digits=4) # limit to 6 digits
tp.cor <- round(cor.test(tp.skcm.schwann.metscore$expression_log2, tp.skcm.schwann.metscore$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

tp.schwann <- ggplot(tp.skcm.schwann.metscore, aes(x=expression_log2, y=metscore, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("SKCM TP (n=80)") +
    xlab("log2(Schwann marker expression)") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", tp.cor), paste0("p-value: ", tp.pval)))

# ---- ~~ GAP43 ----
tp.skcm.gap43 <- tp.skcm.schwann.metscore[which(tp.skcm.schwann.metscore$gene=="GAP43"),]
cor.test(tp.skcm.gap43$expression_log2, tp.skcm.gap43$metscore)
tp.pval <- round(cor.test(tp.skcm.gap43$expression_log2, tp.skcm.gap43$metscore)[3][[1]], digits=4) # limit to 6 digits
tp.cor <- round(cor.test(tp.skcm.gap43$expression_log2, tp.skcm.gap43$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tp.skcm.gap43, aes(x=expression_log2, y=metscore)) +
    geom_point(shape=16, size=5, color=colors.trans[1]) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("SKCM TP (n=80)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", tp.cor), paste0("p-value: ", tp.pval)))

# ---- ~~ GFAP ----
tp.skcm.gfap <- tp.skcm.schwann.metscore[which(tp.skcm.schwann.metscore$gene=="GFAP"),]
cor.test(tp.skcm.gfap$expression_log2, tp.skcm.gfap$metscore)
tp.pval <- round(cor.test(tp.skcm.gfap$expression_log2, tp.skcm.gfap$metscore)[3][[1]], digits=4) # limit to 6 digits
tp.cor <- round(cor.test(tp.skcm.gfap$expression_log2, tp.skcm.gfap$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tp.skcm.gfap, aes(x=expression_log2, y=metscore)) +
    geom_point(shape=16, size=5, color=colors.trans[2]) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("SKCM TP (n=80)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", tp.cor), paste0("p-value: ", tp.pval)))

# ---- ~~ SOX10 ----
tp.skcm.sox10 <- tp.skcm.schwann.metscore[which(tp.skcm.schwann.metscore$gene=="SOX10"),]
cor.test(tp.skcm.sox10$expression_log2, tp.skcm.sox10$metscore)
tp.pval <- round(cor.test(tp.skcm.sox10$expression_log2, tp.skcm.sox10$metscore)[3][[1]], digits=4) # limit to 6 digits
tp.cor <- round(cor.test(tp.skcm.sox10$expression_log2, tp.skcm.sox10$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tp.skcm.sox10, aes(x=expression_log2, y=metscore)) +
    geom_point(shape=16, size=5, color=colors.trans[3]) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("SKCM TP (n=80)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", tp.cor), paste0("p-value: ", tp.pval)))

# ---- ~~ GAP43/GFAP ----
tp.skcm.gap43gfap <- tp.skcm.schwann.metscore[which(tp.skcm.schwann.metscore$gene=="GAP43" | tp.skcm.schwann.metscore$gene=="GFAP"),]
cor.test(tp.skcm.gap43gfap$expression_log2, tp.skcm.gap43gfap$metscore)
tp.pval <- round(cor.test(tp.skcm.gap43gfap$expression_log2, tp.skcm.gap43gfap$metscore)[3][[1]], digits=4) # limit to 6 digits
tp.cor <- round(cor.test(tp.skcm.gap43gfap$expression_log2, tp.skcm.gap43gfap$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tp.skcm.gap43gfap, aes(x=expression_log2, y=metscore, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("SKCM TP (n=80)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", tp.cor), paste0("p-value: ", tp.pval)))


# ---- . SKCM TM ---- 
# 604 outliers removed
tm.skcm.schwann.metscore <- merge(tm.skcm.expr.schwann, tm.skcm.metscore[c(1,3)], by="tcga_participant_barcode")
tm.q <- quantile(tm.skcm.schwann.metscore$expression_log2, na.rm=T)
tm.skcm.schwann.metscore <- tm.skcm.schwann.metscore[which(tm.skcm.schwann.metscore$expression_log2>tm.q[2] & tm.skcm.schwann.metscore$expression_log2<tm.q[4]),]
cor.test(tm.skcm.schwann.metscore$expression_log2, tm.skcm.schwann.metscore$metscore)
tm.pval <- round(cor.test(tm.skcm.schwann.metscore$expression_log2, tm.skcm.schwann.metscore$metscore)[3][[1]], digits=4) # limit to 6 digits
tm.cor <- round(cor.test(tm.skcm.schwann.metscore$expression_log2, tm.skcm.schwann.metscore$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

tm.schwann <- ggplot(tm.skcm.schwann.metscore, aes(x=expression_log2, y=metscore, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("SKCM TM (n=300)") +
    xlab("log2(Schwann marker expression)") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.765,13.3), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", tm.cor), paste0("p-value: ", tm.pval)))

# ---- ~~ GAP43 ----
tm.skcm.gap43 <- tm.skcm.schwann.metscore[which(tm.skcm.schwann.metscore$gene=="GAP43"),]
cor.test(tm.skcm.gap43$expression_log2, tm.skcm.gap43$metscore)
tm.pval <- round(cor.test(tm.skcm.gap43$expression_log2, tm.skcm.gap43$metscore)[3][[1]], digits=4) # limit to 6 digits
tm.cor <- round(cor.test(tm.skcm.gap43$expression_log2, tm.skcm.gap43$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tm.skcm.gap43, aes(x=expression_log2, y=metscore)) +
    geom_point(shape=16, size=5, color=colors.trans[1]) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("SKCM TM (n=300)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", tm.cor), paste0("p-value: ", tm.pval)))

# ---- ~~ GFAP ----
tm.skcm.gfap <- tm.skcm.schwann.metscore[which(tm.skcm.schwann.metscore$gene=="GFAP"),]
cor.test(tm.skcm.gfap$expression_log2, tm.skcm.gfap$metscore)
tm.pval <- round(cor.test(tm.skcm.gfap$expression_log2, tm.skcm.gfap$metscore)[3][[1]], digits=4) # limit to 6 digits
tm.cor <- round(cor.test(tm.skcm.gfap$expression_log2, tm.skcm.gfap$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tm.skcm.gfap, aes(x=expression_log2, y=metscore)) +
    geom_point(shape=16, size=5, color=colors.trans[2]) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("SKCM TM (n=300)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", tm.cor), paste0("p-value: ", tm.pval)))

# ---- ~~ SOX10 ----
tm.skcm.sox10 <- tm.skcm.schwann.metscore[which(tm.skcm.schwann.metscore$gene=="SOX10"),]
cor.test(tm.skcm.sox10$expression_log2, tm.skcm.sox10$metscore)
tm.pval <- round(cor.test(tm.skcm.sox10$expression_log2, tm.skcm.sox10$metscore)[3][[1]], digits=4) # limit to 6 digits
tm.cor <- round(cor.test(tm.skcm.sox10$expression_log2, tm.skcm.sox10$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tm.skcm.sox10, aes(x=expression_log2, y=metscore)) +
    geom_point(shape=16, size=5, color=colors.trans[3]) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("SKCM TM (n=300)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", tm.cor), paste0("p-value: ", tm.pval)))


# ---- ~~ GAP43/GFAP ----
tm.skcm.gap43gfap <- tm.skcm.schwann.metscore[which(tm.skcm.schwann.metscore$gene=="GAP43" | tm.skcm.schwann.metscore$gene=="GFAP"),]
cor.test(tm.skcm.gap43gfap$expression_log2, tm.skcm.gap43gfap$metscore)
tm.pval <- round(cor.test(tm.skcm.gap43gfap$expression_log2, tm.skcm.gap43gfap$metscore)[3][[1]], digits=5) # limit to 6 digits
tm.cor <- round(cor.test(tm.skcm.gap43gfap$expression_log2, tm.skcm.gap43gfap$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tm.skcm.gap43gfap, aes(x=expression_log2, y=metscore, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("SKCM TM (n=300)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,12.85), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", tm.cor), paste0("p-value: ", tm.pval)))


# ---- . UVM ---- 
# 122 outliers removed
uvm.schwann.metscore <- merge(uvm.expr.schwann, uvm.metscore[c(1,3)], by="tcga_participant_barcode")
q <- quantile(uvm.schwann.metscore$expression_log2, na.rm=T)
uvm.schwann.metscore <- uvm.schwann.metscore[which(uvm.schwann.metscore$expression_log2>q[2] & uvm.schwann.metscore$expression_log2<q[4]),]
cor.test(uvm.schwann.metscore$expression_log2, uvm.schwann.metscore$metscore)
pval <- round(cor.test(uvm.schwann.metscore$expression_log2, uvm.schwann.metscore$metscore)[3][[1]], digits=4) # limit to 6 digits
cor <- round(cor.test(uvm.schwann.metscore$expression_log2, uvm.schwann.metscore$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

uvm.schwann <- ggplot(uvm.schwann.metscore, aes(x=expression_log2, y=metscore, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("UVM TP (n=71)") +
    xlab("log2(Schwann marker expression)") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(2.05,14.2), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.title = element_text(size=13),
          legend.text = element_text(size=12)) +
    annotate("text", x=11, y=c(1.65, 1.58), label = c(paste0("Correlation: ", cor), paste0("p-value: ", pval)))

# ---- ~~ GAP43 ----
uvm.gap43 <- uvm.schwann.metscore[which(uvm.schwann.metscore$gene=="GAP43"),]
cor.test(uvm.gap43$expression_log2, uvm.gap43$metscore)
pval <- round(cor.test(uvm.gap43$expression_log2, uvm.gap43$metscore)[3][[1]], digits=4) # limit to 6 digits
cor <- round(cor.test(uvm.gap43$expression_log2, uvm.gap43$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(uvm.gap43, aes(x=expression_log2, y=metscore)) +
    geom_point(shape=16, size=5, color=colors.trans[1]) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("UVM TP (n=81)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", cor), paste0("p-value: ", pval)))

# ---- ~~ GFAP ----
uvm.gfap <- uvm.schwann.metscore[which(uvm.schwann.metscore$gene=="GFAP"),]
cor.test(uvm.gfap$expression_log2, uvm.gfap$metscore)
pval <- round(cor.test(uvm.gfap$expression_log2, uvm.gfap$metscore)[3][[1]], digits=4) # limit to 6 digits
cor <- round(cor.test(uvm.gfap$expression_log2, uvm.gfap$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(uvm.gfap, aes(x=expression_log2, y=metscore)) +
    geom_point(shape=16, size=5, color=colors.trans[2]) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("UVM TP (n=81)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", cor), paste0("p-value: ", pval)))

# ---- ~~ SOX10 ----
uvm.sox10 <- uvm.schwann.metscore[which(uvm.schwann.metscore$gene=="SOX10"),]
cor.test(uvm.sox10$expression_log2, uvm.sox10$metscore)
pval <- round(cor.test(uvm.sox10$expression_log2, uvm.sox10$metscore)[3][[1]], digits=4) # limit to 6 digits
cor <- round(cor.test(uvm.sox10$expression_log2, uvm.sox10$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(uvm.sox10, aes(x=expression_log2, y=metscore)) +
    geom_point(shape=16, size=5, color=colors.trans[3]) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("UVM TP (n=81)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", cor), paste0("p-value: ", pval)))


# ---- ~~ GAP43/GFAP ----
uvm.gap43gfap <- uvm.schwann.metscore[which(uvm.schwann.metscore$gene=="GAP43" | uvm.schwann.metscore$gene=="GFAP"),]
cor.test(uvm.gap43gfap$expression_log2, uvm.gap43gfap$metscore)
pval <- round(cor.test(uvm.gap43gfap$expression_log2, uvm.gap43gfap$metscore)[3][[1]], digits=4) # limit to 6 digits
cor <- round(cor.test(uvm.gap43gfap$expression_log2, uvm.gap43gfap$metscore)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(uvm.gap43gfap, aes(x=expression_log2, y=metscore, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("UVM TP (n=81)") +
    xlab("Schwann marker expression") +
    ylab("Metastatic score") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0.6,1.7), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(1.65, 1.58), label = c(paste0("Correlation: ", cor), paste0("p-value: ", pval)))



# ---- ~ Plot ----
pdf("./Figures/TCGA_SKCM_Schwann_scatter.pdf", w=8, h=4)
ggarrange(tp.schwann, tm.schwann, ncol=2)
dev.off()

pdf("./Figures/TCGA_UVM_Schwann_scatter.pdf", w=5.2, h=4)
uvm.schwann
dev.off()

#--------- MUTSIGS -------------------------------------------------------------------------------------------------

# ---- ~ Total muts contribution ----

# SKCM TP -- 168 outliers removed

# S1
tp.skcm.schwann.s1 <- merge(tp.skcm.expr.schwann, tp.skcm.total.muts[c(1:2)], by="tcga_participant_barcode")
tp.q <- quantile(tp.skcm.schwann.s1$expression_log2, na.rm=T)
tp.skcm.schwann.s1 <- tp.skcm.schwann.s1[which(tp.skcm.schwann.s1$expression_log2>tp.q[2] & tp.skcm.schwann.s1$expression_log2<tp.q[4]),]
cor.test(tp.skcm.schwann.s1$expression_log2, tp.skcm.schwann.s1$S1)
s1.pval <- round(cor.test(tp.skcm.schwann.s1$expression_log2, tp.skcm.schwann.s1$S1)[3][[1]], digits=4) # limit to 6 digits
s1.cor <- round(cor.test(tp.skcm.schwann.s1$expression_log2, tp.skcm.schwann.s1$S1)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tp.skcm.schwann.s1, aes(x=expression_log2, y=S1, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("S1") +
    xlab("Schwann marker expression") +
    ylab("Total mutation contribution") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(-19,2000,10)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(-1,138), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(130, 120), label = c(paste0("Correlation: ", s1.cor), paste0("p-value: ", s1.pval)))

# S7
tp.skcm.schwann.s7 <- merge(tp.skcm.expr.schwann, tp.skcm.total.muts[c(1,8)], by="tcga_participant_barcode")
tp.q <- quantile(tp.skcm.schwann.s7$expression_log2, na.rm=T)
tp.skcm.schwann.s7 <- tp.skcm.schwann.s7[which(tp.skcm.schwann.s7$expression_log2>tp.q[2] & tp.skcm.schwann.s7$expression_log2<tp.q[4]),]
cor.test(tp.skcm.schwann.s7$expression_log2, tp.skcm.schwann.s7$S7)
s7.pval <- round(cor.test(tp.skcm.schwann.s7$expression_log2, tp.skcm.schwann.s7$S7)[3][[1]], digits=4) # limit to 6 digits
s7.cor <- round(cor.test(tp.skcm.schwann.s7$expression_log2, tp.skcm.schwann.s7$S7)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tp.skcm.schwann.s7, aes(x=expression_log2, y=S7, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("S7") +
    xlab("Schwann marker expression") +
    ylab("Total mutation contribution") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,9000,500)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0,3000), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(2900, 2700), label = c(paste0("Correlation: ", s7.cor), paste0("p-value: ", s7.pval)))


# S6
tp.skcm.schwann.s6 <- merge(tp.skcm.expr.schwann, tp.skcm.total.muts[c(1,7)], by="tcga_participant_barcode")
tp.q <- quantile(tp.skcm.schwann.s6$expression_log2, na.rm=T)
tp.skcm.schwann.s6 <- tp.skcm.schwann.s6[which(tp.skcm.schwann.s6$expression_log2>tp.q[2] & tp.skcm.schwann.s6$expression_log2<tp.q[4]),]
cor.test(tp.skcm.schwann.s6$expression_log2, tp.skcm.schwann.s6$S6)
s6.pval <- round(cor.test(tp.skcm.schwann.s6$expression_log2, tp.skcm.schwann.s6$S6)[3][[1]], digits=4) # limit to 6 digits
s6.cor <- round(cor.test(tp.skcm.schwann.s6$expression_log2, tp.skcm.schwann.s6$S6)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tp.skcm.schwann.s6, aes(x=expression_log2, y=S6, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("S6") +
    xlab("Schwann marker expression") +
    ylab("Total mutation contribution") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,9000,500)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0,3000), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(2900, 2700), label = c(paste0("Correlation: ", s6.cor), paste0("p-value: ", s6.pval)))


# SKCM TP -- 604 outliers removed

# S1
tm.skcm.schwann.s1 <- merge(tm.skcm.expr.schwann, tm.skcm.total.muts[c(1:2)], by="tcga_participant_barcode")
tm.q <- quantile(tm.skcm.schwann.s1$expression_log2, na.rm=T)
tm.skcm.schwann.s1 <- tm.skcm.schwann.s1[which(tm.skcm.schwann.s1$expression_log2>tm.q[2] & tm.skcm.schwann.s1$expression_log2<tm.q[4]),]
cor.test(tm.skcm.schwann.s1$expression_log2, tm.skcm.schwann.s1$S1)
s1.pval <- round(cor.test(tm.skcm.schwann.s1$expression_log2, tm.skcm.schwann.s1$S1)[3][[1]], digits=4) # limit to 6 digits
s1.cor <- round(cor.test(tm.skcm.schwann.s1$expression_log2, tm.skcm.schwann.s1$S1)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tm.skcm.schwann.s1, aes(x=expression_log2, y=S1, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("S1") +
    xlab("Schwann marker expression") +
    ylab("Total mutation contribution") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(-19,2000,10)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(-1,138), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(130, 120), label = c(paste0("Correlation: ", s1.cor), paste0("p-value: ", s1.pval)))

# S7
tm.skcm.schwann.s7 <- merge(tm.skcm.expr.schwann, tm.skcm.total.muts[c(1,8)], by="tcga_participant_barcode")
tm.q <- quantile(tm.skcm.schwann.s7$expression_log2, na.rm=T)
tm.skcm.schwann.s7 <- tm.skcm.schwann.s7[which(tm.skcm.schwann.s7$expression_log2>tm.q[2] & tm.skcm.schwann.s7$expression_log2<tm.q[4]),]
cor.test(tm.skcm.schwann.s7$expression_log2, tm.skcm.schwann.s7$S7)
s7.pval <- round(cor.test(tm.skcm.schwann.s7$expression_log2, tm.skcm.schwann.s7$S7)[3][[1]], digits=4) # limit to 6 digits
s7.cor <- round(cor.test(tm.skcm.schwann.s7$expression_log2, tm.skcm.schwann.s7$S7)[4][[1]][[1]], digits=7) # limit to 8 digits

ggplot(tm.skcm.schwann.s7, aes(x=expression_log2, y=S7, color=gene)) +
    geom_point(shape=16, size=5) + 
    geom_smooth(color="black", method=lm) + 
    scale_color_manual(values=colors.trans, name="Gene") +
    ggtitle("S7") +
    xlab("Schwann marker expression") +
    ylab("Total contribution") +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,9000,500)) +
    coord_cartesian(xlim=c(0.711,13.5), ylim=c(0,3000), expand=FALSE) +
    theme_bw() + 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.position = "none") +
    annotate("text", x=10, y=c(2900, 2700), label = c(paste0("Correlation: ", s7.cor), paste0("p-value: ", s7.pval)))


