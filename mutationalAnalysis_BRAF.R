# SKCM Cancer BRAF
# DT Laura Vo Ngoc
# Start: 15/06/2018

library("BSgenome.Hsapiens.UCSC.hg38")
library(deconstructSigs)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(TCGA2STAT)
library(TCGAbiolinks)
require(FirebrowseR)

#--------- WD & LOAD FILES --------------------------------------------------------------------------------------------
setwd("C:/Users/rockp/Desktop/UCL/Project/Project-Laura")

load("./Data/TCGA_SKCM_mutations.RData")            # Raw TCGA mutation data
load("./Output/TCGA_SKCM_snvs.RData")               # Selected SNVs

# Primary tumor
load("./Output/TCGA_SKCM_TP_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_SKCM_TP_weights_cut0.00.RData")     # Mutational signatures output

# Metastatic tumor
load("./Output/TCGA_SKCM_TM_mutsig_input.RData")            # Mutation analysis input
load("./Output/TCGA_SKCM_TM_weights_cut0.00.RData")     # Mutational signatures output

# Mutsigs, metscores, and stage
load("./Output/TCGA_SKCM_TP_metastatic_score_stage.RData")  # SKCM TP
load("./Output/TCGA_SKCM_TM_metastatic_score_stage.RData")  # SKCM TM
load("./Output/TCGA_UVM_TP_metastatic_score_mel.RData")     # UVM TP

# Mutation data
load("./Data/TCGA_SKCM_BRAFmut.RData")          # BRAF mutation data SKCM ---> mut.skcm
load("./Data/TCGA_UVM_ALLmut.RData")            # All mutation data UVM ---> mut.uvm

# Total mutsig mutations
load("./Output/TCGA_SKCM_TP_totalmutsigs.RData")    # SKCM TP
load("./Output/TCGA_SKCM_TM_totalmutsigs.RData")    # SKCM TM

# BRAF
load("./Output/TCGA_SKCM_BRAFmut.RData")        # BRAF mutation data SKCM df --> braf.skcm
load("./Output/TCGA_SKCM_TP_BRAFmut.RData")     # BRAF mutation data SKCM TP --> tp.braf
load("./Output/TCGA_SKCM_TM_BRAFmut.RData")     # BRAF mutation data SKCM TM --> tm.braf

#---------  DATA OF INTEREST ------------------------------------------------------------------------------------------

# Sys.setenv(TAR="C:/cygwin64/bin/tar", R_GZIPCMD="C:/cygwin64/bin/gzip")
# mut.skcm <- getTCGA(disease="SKCM", data.type="Mutation", type="somatic")       # can't get it to work --> Maria downloaded

braf.skcm <- data.frame(tcga_participant_barcode=names(mut.skcm), BRAF=mut.skcm)
save(braf.skcm, file="./Output/TCGA_SKCM_BRAFmut.RData")

tp.braf <- merge(tp.skcm.total.muts, braf.skcm, by="tcga_participant_barcode", all.x=TRUE)
tm.braf <- merge(tm.skcm.total.muts, braf.skcm, by="tcga_participant_barcode", all.x=TRUE)
save(tp.braf, file="./Output/TCGA_SKCM_TP_BRAFmut.RData")
save(tm.braf, file="./Output/TCGA_SKCM_TM_BRAFmut.RData")

# Split into BRAF WT and mutated cohorts
tp.braf.wt <- tp.braf[which(tp.braf$BRAF==0),]
tp.braf.mut <- tp.braf[which(tp.braf$BRAF==1),]

tm.braf.wt <- tm.braf[which(tm.braf$BRAF==0),]
tm.braf.mut <- tm.braf[which(tm.braf$BRAF==1),]



#---------  BRAF WT VS. MUT ------------------------------------------------------------------------------------------

#--------- ~ SKCM TP --------

#--------- ~~~ S7 --------
# WT plot -- S7 total muts vs. metscore -- 1 outlier removed
wt.s7 <- tp.braf.wt[which(tp.braf.wt$S7<2000),]
quantile(wt.s7$S7)
cor.test(wt.s7$S7, wt.s7$metscore)

pdf("./Figures/TCGA_SKCM_TP_BRAF_WT_S7_metscore_scatter.pdf", w=8, h=6)
ggplot(wt.s7, aes(x=S7, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TP BRAF WT") +
    xlab("Total mutations contributed by S7") +
    ylab("Metastatic score") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,2000,50)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,1976), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=1750, y=c(1.47, 1.44), label = c("Correlation: 0.4614446", "p-value: 0.02322"))
dev.off()

wt.s7$met_potential <- factor(wt.s7$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TP_BRAF_WT_S7_metpotential_boxplot.pdf", w=6, h=6)
ggplot(wt.s7, aes(x=met_potential, y=S7)) +
    geom_boxplot() +
    ggtitle("SKCM TP BRAF WT") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S7") +
    scale_y_continuous(breaks=seq(0,2000,100)) +
    coord_cartesian(ylim=c(0,2000), expand=FALSE)
#annotate("text", x=120, y=c(1.47, 1.44), label = c("Correlation: 0.1649576", "p-value: 0.1045"))
dev.off()

wilcox.test(wt.s7[which(mut.s7$met_potential=="low"),8],wt.s7[which(mut.s7$met_potential=="high"),8])

# Mut plot -- S7 total muts vs. metscore -- 2 outliers removed
mut.s7 <- tp.braf.mut[which(tp.braf.mut$S7<1800),]
quantile(mut.s7$S7)
cor.test(mut.s7$S7, mut.s7$metscore)

pdf("./Figures/TCGA_SKCM_TP_BRAF_mut_S7_metscore_scatter.pdf", w=8, h=6)
ggplot(mut.s7, aes(x=S7, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TP BRAF mut") +
    xlab("Total mutations contributed by S7") +
    ylab("Metastatic score") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,1000,50)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,711), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=620, y=c(1.47, 1.44), label = c("Correlation: 0.1054162", "p-value: 0.5658"))
dev.off()

mut.s7$met_potential <- factor(mut.s7$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TP_BRAF_mut_S7_metpotential_boxplot.pdf", w=6, h=6)
ggplot(mut.s7, aes(x=met_potential, y=S7)) +
    geom_boxplot() +
    ggtitle("SKCM TP BRAF mut") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S7") +
    scale_y_continuous(breaks=seq(0,2000,50)) +
    coord_cartesian(ylim=c(0,1000), expand=FALSE)
#annotate("text", x=120, y=c(1.47, 1.44), label = c("Correlation: 0.1649576", "p-value: 0.1045"))
dev.off()

wilcox.test(mut.s7[which(mut.s7$met_potential=="low"),8],mut.s7[which(mut.s7$met_potential=="high"),8])


#--------- ~~~ S1 --------
# WT plot -- S1 total muts vs. metscore
wt.s1 <- tp.braf.wt
max <- quantile(wt.s1$S1)[5][[1]]
# Correlation
cor.test(wt.s1$S1, wt.s1$metscore)
pval <- round(cor.test(wt.s1$S1, wt.s1$metscore)[3][[1]], digits=6) # limit to 6 digits
cor <- round(cor.test(wt.s1$S1, wt.s1$metscore)[4][[1]][[1]], digits=8) # limit to 8 digits
lab <- c(paste0("Correlation: ", cor), paste0("p-value: ", pval))


pdf("./Figures/TCGA_SKCM_TP_BRAF_WT_S1_metscore_scatter.pdf", w=8, h=6)
ggplot(wt.s1, aes(x=S1, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TP BRAF WT") +
    xlab("Total mutations contributed by S1") +
    ylab("Metastatic score") +
    #theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,max+0.1), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=max-10, y=c(1.47, 1.44), label = lab)
dev.off()

wt.s1$met_potential <- factor(wt.s1$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TP_BRAF_WT_S1_metpotential_boxplot.pdf", w=6, h=6)
ggplot(wt.s1, aes(x=met_potential, y=S1)) +
    geom_boxplot() +
    ggtitle("SKCM TP BRAF WT") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S1") +
    scale_y_continuous(breaks=seq(0,max+10,2)) +
    coord_cartesian(ylim=c(0,max+10), expand=FALSE)
dev.off()

wilcox.test(wt.s1[which(wt.s1$met_potential=="low"),2],wt.s1[which(wt.s1$met_potential=="high"),2])

# Mut plot -- S1 total muts vs. metscore
mut.s1 <- tp.braf.mut
max <- quantile(mut.s1$S1)[5][[1]]
# Correlation
cor.test(mut.s1$S1, mut.s1$metscore)
pval <- round(cor.test(mut.s1$S1, mut.s1$metscore)[3][[1]], digits=6) # limit to 6 digits
cor <- round(cor.test(mut.s1$S1, mut.s1$metscore)[4][[1]][[1]], digits=8) # limit to 8 digits
lab <- c(paste0("Correlation: ", cor), paste0("p-value: ", pval))


pdf("./Figures/TCGA_SKCM_TP_BRAF_mut_S1_metscore_scatter.pdf", w=8, h=6)
ggplot(mut.s1, aes(x=S1, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TP BRAF mut") +
    xlab("Total mutations contributed by S1") +
    ylab("Metastatic score") +
    #theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,max+0.1), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=max-10, y=c(1.47, 1.44), label = lab)
dev.off()

mut.s1$met_potential <- factor(mut.s1$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TP_BRAF_mut_S1_metpotential_boxplot.pdf", w=6, h=6)
ggplot(mut.s1, aes(x=met_potential, y=S1)) +
    geom_boxplot() +
    ggtitle("SKCM TP BRAF mut") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S1") +
    scale_y_continuous(breaks=seq(0,max+10,2)) +
    coord_cartesian(ylim=c(0,max+10), expand=FALSE)
dev.off()

wilcox.test(mut.s1[which(mut.s1$met_potential=="low"),2],mut.s1[which(mut.s1$met_potential=="high"),2])


#--------- ~~~ S11 --------

# WT plot -- S11 total muts vs. metscore -- 1 outlier removed
wt.s11 <- tp.braf.wt[which(tp.braf.wt$S11<200),]
max <- quantile(wt.s11$S11)[5][[1]]
# Correlation
cor.test(wt.s11$S11, wt.s11$metscore)
pval <- round(cor.test(wt.s11$S11, wt.s11$metscore)[3][[1]], digits=6) # limit to 6 digits
cor <- round(cor.test(wt.s11$S11, wt.s11$metscore)[4][[1]][[1]], digits=8) # limit to 8 digits
lab <- c(paste0("Correlation: ", cor), paste0("p-value: ", pval))


pdf("./Figures/TCGA_SKCM_TP_BRAF_WT_S11_metscore_scatter.pdf", w=8, h=6)
ggplot(wt.s11, aes(x=S11, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TP BRAF WT") +
    xlab("Total mutations contributed by S11") +
    ylab("Metastatic score") +
    #theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,2000,1)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,max+0.1), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=max-6, y=c(1.47, 1.44), label = lab)
dev.off()

wt.s11$met_potential <- factor(wt.s11$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TP_BRAF_WT_S11_metpotential_boxplot.pdf", w=6, h=6)
ggplot(wt.s11, aes(x=met_potential, y=S11)) +
    geom_boxplot() +
    ggtitle("SKCM TP BRAF WT") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S11") +
    scale_y_continuous(breaks=seq(0,max+10,1)) +
    coord_cartesian(ylim=c(0,max+4), expand=FALSE)
dev.off()

wilcox.test(wt.s11[which(wt.s11$met_potential=="low"),12],wt.s11[which(wt.s11$met_potential=="high"),12])

# Mut plot -- S11 total muts vs. metscore
mut.s11 <- tp.braf.mut
max <- quantile(mut.s11$S11)[5][[1]]
# Correlation
cor.test(mut.s11$S11, mut.s11$metscore)
pval <- round(cor.test(mut.s11$S11, mut.s11$metscore)[3][[1]], digits=6) # limit to 6 digits
cor <- round(cor.test(mut.s11$S11, mut.s11$metscore)[4][[1]][[1]], digits=8) # limit to 8 digits
lab <- c(paste0("Correlation: ", cor), paste0("p-value: ", pval))


pdf("./Figures/TCGA_SKCM_TP_BRAF_mut_S11_metscore_scatter.pdf", w=8, h=6)
ggplot(mut.s11, aes(x=S11, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TP BRAF mut") +
    xlab("Total mutations contributed by S11") +
    ylab("Metastatic score") +
    #theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,2000,2)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,max+0.1), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=max-11, y=c(1.47, 1.44), label = lab)
dev.off()

mut.s11$met_potential <- factor(mut.s11$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TP_BRAF_mut_S11_metpotential_boxplot.pdf", w=6, h=6)
ggplot(mut.s11, aes(x=met_potential, y=S11)) +
    geom_boxplot() +
    ggtitle("SKCM TP BRAF mut") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S11") +
    scale_y_continuous(breaks=seq(0,max+10,2)) +
    coord_cartesian(ylim=c(0,max+10), expand=FALSE)
dev.off()

wilcox.test(mut.s11[which(mut.s11$met_potential=="low"),12],mut.s11[which(mut.s11$met_potential=="high"),12])


#--------- ~ SKCM TM --------

#--------- ~~~ S7 --------
# WT plot -- S7 total muts vs. metscore
wt.s7 <- tm.braf.wt
max <- quantile(wt.s7$S7)[5][[1]]
# Correlation
cor.test(wt.s7$S7, wt.s7$metscore)
pval <- round(cor.test(wt.s7$S7, wt.s7$metscore)[3][[1]], digits=6) # limit to 6 digits
cor <- round(cor.test(wt.s7$S7, wt.s7$metscore)[4][[1]][[1]], digits=8) # limit to 8 digits
lab <- c(paste0("Correlation: ", cor), paste0("p-value: ", pval))


pdf("./Figures/TCGA_SKCM_TM_BRAF_WT_S7_metscore_scatter.pdf", w=8, h=6)
ggplot(wt.s7, aes(x=S7, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TM BRAF WT") +
    xlab("Total mutations contributed by S7") +
    ylab("Metastatic score") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,5000,100)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,max+0.1), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=max-500, y=c(1.47, 1.44), label = lab)
dev.off()

wt.s7$met_potential <- factor(wt.s7$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TM_BRAF_WT_S7_metpotential_boxplot.pdf", w=6, h=6)
ggplot(wt.s7, aes(x=met_potential, y=S7)) +
    geom_boxplot() +
    ggtitle("SKCM TM BRAF WT") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S7") +
    scale_y_continuous(breaks=seq(0,max+10,100)) +
    coord_cartesian(ylim=c(0,max+10), expand=FALSE)
dev.off()

wilcox.test(wt.s7[which(wt.s7$met_potential=="low"),8],wt.s7[which(wt.s7$met_potential=="high"),8])


# mut plot -- S7 total muts vs. metscore -- 1 outlier removed
mut.s7 <- tm.braf.mut[which(tm.braf.mut$S7<9000),]
max <- quantile(mut.s7$S7)[5][[1]]
# Correlation
cor.test(mut.s7$S7, mut.s7$metscore)
pval <- round(cor.test(mut.s7$S7, mut.s7$metscore)[3][[1]], digits=6) # limit to 6 digits
cor <- round(cor.test(mut.s7$S7, mut.s7$metscore)[4][[1]][[1]], digits=8) # limit to 8 digits
lab <- c(paste0("Correlation: ", cor), paste0("p-value: ", pval))


pdf("./Figures/TCGA_SKCM_TM_BRAF_mut_S7_metscore_scatter.pdf", w=8, h=6)
ggplot(mut.s7, aes(x=S7, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TM BRAF mut") +
    xlab("Total mutations contributed by S7") +
    ylab("Metastatic score") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,100000,100)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,max+0.1), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=max-500, y=c(1.47, 1.44), label = lab)
dev.off()

mut.s7$met_potential <- factor(mut.s7$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TM_BRAF_mut_S7_metpotential_boxplot.pdf", w=6, h=6)
ggplot(mut.s7, aes(x=met_potential, y=S7)) +
    geom_boxplot() +
    ggtitle("SKCM TM BRAF mut") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S7") +
    scale_y_continuous(breaks=seq(0,max+10,100)) +
    coord_cartesian(ylim=c(0,max+10), expand=FALSE)
dev.off()

wilcox.test(mut.s7[which(mut.s7$met_potential=="low"),8],mut.s7[which(mut.s7$met_potential=="high"),8])


#--------- ~~~ S1 --------
# WT plot -- S1 total muts vs. metscore
wt.s1 <- tm.braf.wt
max <- quantile(wt.s1$S1)[5][[1]]
# Correlation
cor.test(wt.s1$S1, wt.s1$metscore)
pval <- round(cor.test(wt.s1$S1, wt.s1$metscore)[3][[1]], digits=6) # limit to 6 digits
cor <- round(cor.test(wt.s1$S1, wt.s1$metscore)[4][[1]][[1]], digits=8) # limit to 8 digits
lab <- c(paste0("Correlation: ", cor), paste0("p-value: ", pval))


pdf("./Figures/TCGA_SKCM_TM_BRAF_WT_S1_metscore_scatter.pdf", w=8, h=6)
ggplot(wt.s1, aes(x=S1, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TM BRAF WT") +
    xlab("Total mutations contributed by S1") +
    ylab("Metastatic score") +
    #theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,5000,10)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,max+0.1), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=max-50, y=c(1.47, 1.44), label = lab)
dev.off()

wt.s1$met_potential <- factor(wt.s1$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TM_BRAF_WT_S1_metpotential_boxplot.pdf", w=6, h=6)
ggplot(wt.s1, aes(x=met_potential, y=S1)) +
    geom_boxplot() +
    ggtitle("SKCM TM BRAF WT") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S1") +
    scale_y_continuous(breaks=seq(0,max+10,10)) +
    coord_cartesian(ylim=c(0,max+10), expand=FALSE)
dev.off()

wilcox.test(wt.s1[which(wt.s1$met_potential=="low"),2],wt.s1[which(wt.s1$met_potential=="high"),2])


# mut plot -- S1 total muts vs. metscore
mut.s1 <- tm.braf.mut
max <- quantile(mut.s1$S1)[5][[1]]
# Correlation
cor.test(mut.s1$S1, mut.s1$metscore)
pval <- round(cor.test(mut.s1$S1, mut.s1$metscore)[3][[1]], digits=6) # limit to 6 digits
cor <- round(cor.test(mut.s1$S1, mut.s1$metscore)[4][[1]][[1]], digits=8) # limit to 8 digits
lab <- c(paste0("Correlation: ", cor), paste0("p-value: ", pval))


pdf("./Figures/TCGA_SKCM_TM_BRAF_mut_S1_metscore_scatter.pdf", w=8, h=6)
ggplot(mut.s1, aes(x=S1, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TM BRAF mut") +
    xlab("Total mutations contributed by S1") +
    ylab("Metastatic score") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,100000,50)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,max+0.1), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=max-100, y=c(1.47, 1.44), label = lab)
dev.off()

mut.s1$met_potential <- factor(mut.s1$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TM_BRAF_mut_S1_metpotential_boxplot.pdf", w=6, h=6)
ggplot(mut.s1, aes(x=met_potential, y=S1)) +
    geom_boxplot() +
    ggtitle("SKCM TM BRAF mut") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S1") +
    scale_y_continuous(breaks=seq(0,max+10,10)) +
    coord_cartesian(ylim=c(0,max+10), expand=FALSE)
dev.off()

wilcox.test(mut.s1[which(mut.s1$met_potential=="low"),2],mut.s1[which(mut.s1$met_potential=="high"),2])


#--------- ~~~ S11 --------
# WT plot -- S11 total muts vs. metscore
wt.s11 <- tm.braf.wt
max <- quantile(wt.s11$S11)[5][[1]]
# Correlation
cor.test(wt.s11$S11, wt.s11$metscore)
pval <- round(cor.test(wt.s11$S11, wt.s11$metscore)[3][[1]], digits=6) # limit to 6 digits
cor <- round(cor.test(wt.s11$S11, wt.s11$metscore)[4][[1]][[1]], digits=8) # limit to 8 digits
lab <- c(paste0("Correlation: ", cor), paste0("p-value: ", pval))


pdf("./Figures/TCGA_SKCM_TM_BRAF_WT_S11_metscore_scatter.pdf", w=8, h=6)
ggplot(wt.s11, aes(x=S11, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TM BRAF WT") +
    xlab("Total mutations contributed by S11") +
    ylab("Metastatic score") +
    #theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,5000,10)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,max+0.1), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=max-30, y=c(1.47, 1.44), label = lab)
dev.off()

wt.s11$met_potential <- factor(wt.s11$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TM_BRAF_WT_S11_metpotential_boxplot.pdf", w=6, h=6)
ggplot(wt.s11, aes(x=met_potential, y=S11)) +
    geom_boxplot() +
    ggtitle("SKCM TM BRAF WT") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S11") +
    scale_y_continuous(breaks=seq(0,max+10,10)) +
    coord_cartesian(ylim=c(0,max+10), expand=FALSE)
dev.off()

wilcox.test(wt.s11[which(wt.s11$met_potential=="low"),12],wt.s11[which(wt.s11$met_potential=="high"),12])


# mut plot -- S11 total muts vs. metscore
mut.s11 <- tm.braf.mut
max <- quantile(mut.s11$S11)[5][[1]]
# Correlation
cor.test(mut.s11$S11, mut.s11$metscore)
pval <- round(cor.test(mut.s11$S11, mut.s11$metscore)[3][[1]], digits=6) # limit to 6 digits
cor <- round(cor.test(mut.s11$S11, mut.s11$metscore)[4][[1]][[1]], digits=8) # limit to 8 digits
lab <- c(paste0("Correlation: ", cor), paste0("p-value: ", pval))


pdf("./Figures/TCGA_SKCM_TM_BRAF_mut_S11_metscore_scatter.pdf", w=8, h=6)
ggplot(mut.s11, aes(x=S11, y=metscore)) +
    geom_point(shape=1) + 
    geom_smooth(method=lm) + 
    ggtitle("SKCM TM BRAF mut") +
    xlab("Total mutations contributed by S11") +
    ylab("Metastatic score") +
    #theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    scale_x_continuous(breaks=seq(0,100000,10)) +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    coord_cartesian(xlim=c(0,max+0.1), ylim=c(0.6,1.5), expand=FALSE)+
    annotate("text", x=max-30, y=c(1.47, 1.44), label = lab)
dev.off()

mut.s11$met_potential <- factor(mut.s11$met_potential, levels = c("low","high"),ordered = TRUE)

pdf("./Figures/TCGA_SKCM_TM_BRAF_mut_S11_metpotential_boxplot.pdf", w=6, h=6)
ggplot(mut.s11, aes(x=met_potential, y=S11)) +
    geom_boxplot() +
    ggtitle("SKCM TM BRAF mut") +
    xlab("Metastatic potential") +
    ylab("Total mutations contributed by S11") +
    scale_y_continuous(breaks=seq(0,max+10,10)) +
    coord_cartesian(ylim=c(0,max+10), expand=FALSE)
dev.off()

wilcox.test(mut.s11[which(mut.s11$met_potential=="low"),12],mut.s11[which(mut.s11$met_potential=="high"),12])


# try to see what get with these numbers
# then go back to dataset with all the patients and score those that have expr data available
# to see if braf mut has effect
# then can include mutsigs