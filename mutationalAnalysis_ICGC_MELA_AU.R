# ICGC MELA-AU data analysis
# DT Laura Vo Ngoc
# Start: 19/07/2018

library("BSgenome.Hsapiens.UCSC.hg19")
library(deconstructSigs)
library(ggplot2)
library(MASS)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(TCGAbiolinks)
require(FirebrowseR)

#--------- WD & LOAD FILES --------------------------------------------------------------------------------------------
setwd("C:/Users/rockp/Desktop/UCL/Project/Project-Laura")

# Clinical data
load("./Data/ICGC_MELA_AU_clinical.RData")

# Sample data
load("./Data/ICGC_MELA_AU_sample.RData")

# Risk factor exposure data
load("./Data/ICGC_MELA_AU_exposure.RData")

# Therapy data
load("./Data/ICGC_MELA_AU_therapy.RData")

# SNV data (from mutations file)
load("./Output/ICGC_MELA_AU_50snvs.RData")
load("./Output/ICGC_MELA_AU_100snvs.RData")
load("./Output/ICGC_MELA_AU_100_snvs.RData")

# mutsig input
load("./Output/ICGC_MELA_AU_mutsig_input.RData")
load("./Output/ICGC_MELA_AU_100_mutsig_input.RData")
load("./Output/ICGC_MELA_AU_100_2_mutsig_input.RData")

# mutsig data
load("./Output/ICGC_MELA_AU_weights_cut0.00.RData")
load("./Output/ICGC_MELA_AU__50_weights_cut0.00.RData")
load("./Output/ICGC_MELA_AU__100_weights_cut0.00.RData")
load("./Output/ICGC_MELA_AU_100_2_weights_cut0.00.RData")


# Mutsigs + clinical + therapy data
load("./Output/ICGC_MELA_AU_clin_ther_mutsig.RData")


#---------  IMPORT DATA -----------------------------------------------------------------------------------------------

mela.clin <- read.table(file = './ICGC_MELA_AU/donor.MELA-AU.tsv', sep = '\t', header = TRUE, stringAsFactors=FALSE)
save(mela.clin, file="./Data/ICGC_MELA_AU_clinical.RData")

mela.expo <- read.table(file = './ICGC_MELA_AU/donor_exposure.MELA-AU.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
save(mela.expo, file="./Data/ICGC_MELA_AU_exposure.RData")

mela.ther <- read.table(file = './ICGC_MELA_AU/donor_therapy.MELA-AU.tsv', sep = '\t', header = TRUE)
save(mela.ther, file="./Data/ICGC_MELA_AU_therapy.RData")

mela.sample <- read.table(file = './ICGC_MELA_AU/sample.MELA-AU.tsv', sep = '\t', header = TRUE)
save(mela.sample, file="./Data/ICGC_MELA_AU_sample.RData")


mela50.snvs <- read.table(file = './ICGC_MELA_AU/mutfirst50.MELA-AU.short.tsv', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(mela50.snvs) <- c("icgc_mutation_id", "icgc_donor_id", "icgc_sample_id", "chromosome", "chromosome_start", "mutated_from_allele", "mutated_to_allele")
save(mela50.snvs, file="./Output/ICGC_MELA_AU_50snvs.RData")

mela100.snvs <- read.table(file = './ICGC_MELA_AU/mutnext100.MELA-AU.tsv', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(mela100.snvs) <- c("icgc_mutation_id", "icgc_donor_id", "icgc_sample_id", "chromosome", "chromosome_start", "mutated_from_allele", "mutated_to_allele")
save(mela100.snvs, file="./Output/ICGC_MELA_AU_100snvs.RData")

mela100.snvs <- read.table(file = './ICGC_MELA_AU/mutnext100.MELA-AU.tsv', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(mela100.snvs) <- c("icgc_mutation_id", "icgc_donor_id", "icgc_sample_id", "chromosome", "chromosome_start", "mutated_from_allele", "mutated_to_allele")
save(mela100.snvs, file="./Output/ICGC_MELA_AU_100snvs.RData")

mela.100.snvs <- read.table(file = './ICGC_MELA_AU/mutnext100.MELA-AU.2ndtry.tsv', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(mela.100.snvs) <- c("icgc_mutation_id", "icgc_donor_id", "icgc_sample_id", "chromosome", "chromosome_start", "mutated_from_allele", "mutated_to_allele")
save(mela.100.snvs, file="./Output/ICGC_MELA_AU_100_snvs.RData")


#--------- MUTATIONAL SIGNATURES ANALYSIS -------------------------------------------------------------------------

# Add "chr" to chromosome number in mela50.snvs
mela50.snvs$chromosome <- sapply(mela50.snvs$chromosome, function(x) paste0("chr",x))
mela100.snvs$chromosome <- sapply(mela100.snvs$chromosome, function(x) paste0("chr",x))
mela.100.snvs$chromosome <- sapply(mela.100.snvs$chromosome, function(x) paste0("chr",x))

# Copy original snvs var into another one
mela50.snvs.ori <- mela50.snvs
mela100.snvs.ori <- mela100.snvs
mela.100.snvs.ori <- mela.100.snvs

# Remove mitochondrial mutations (48 removed)
mela50.snvs <- mela50.snvs[which(mela50.snvs$chromosome!="chrMT"),]     # 48 removed
mela100.snvs <- mela100.snvs[which(mela100.snvs$chromosome!="chrMT"),]     # 36 removed
mela.100.snvs <- mela.100.snvs[which(mela.100.snvs$chromosome!="chrMT"),]     # 46 removed

# Convert to deconstructSigs input:
mela.sigs.input <- mut.to.sigs.input(mut.ref = mela50.snvs, 
                                   sample.id = "icgc_sample_id", 
                                   chr = "chromosome", 
                                   pos = "chromosome_start", 
                                   ref = "mutated_from_allele", 
                                   alt = "mutated_to_allele",
                                   bsg = BSgenome.Hsapiens.UCSC.hg19)
mela100.sigs.input <- mut.to.sigs.input(mut.ref = mela100.snvs, 
                                     sample.id = "icgc_sample_id", 
                                     chr = "chromosome", 
                                     pos = "chromosome_start", 
                                     ref = "mutated_from_allele", 
                                     alt = "mutated_to_allele",
                                     bsg = BSgenome.Hsapiens.UCSC.hg19)
mela.100.sigs.input <- mut.to.sigs.input(mut.ref = mela.100.snvs, 
                                        sample.id = "icgc_sample_id", 
                                        chr = "chromosome", 
                                        pos = "chromosome_start", 
                                        ref = "mutated_from_allele", 
                                        alt = "mutated_to_allele",
                                        bsg = BSgenome.Hsapiens.UCSC.hg19)

save(mela.sigs.input, file="./Output/ICGC_MELA_AU_mutsig_input.RData")
save(mela100.sigs.input, file="./Output/ICGC_MELA_AU_100_mutsig_input.RData")
save(mela.100.sigs.input, file="./Output/ICGC_MELA_AU_100_2_mutsig_input.RData")


# Signatures
mela.weights <- as.data.frame(t(sapply(rownames(mela.sigs.input), 
                                     function(x) whichSignatures(tumor.ref = mela.sigs.input,
                                                                 signatures.ref = signatures.cosmic,
                                                                 sample.id = x, 
                                                                 contexts.needed = TRUE,
                                                                 signature.cutoff = 0.00,              # default = 0.06
                                                                 tri.counts.method = 'default')$weights)), 
                            row.names=rownames(mela.sigs.input))
mela100.weights <- as.data.frame(t(sapply(rownames(mela100.sigs.input), 
                                       function(x) whichSignatures(tumor.ref = mela100.sigs.input,
                                                                   signatures.ref = signatures.cosmic,
                                                                   sample.id = x, 
                                                                   contexts.needed = TRUE,
                                                                   signature.cutoff = 0.00,              # default = 0.06
                                                                   tri.counts.method = 'default')$weights)), 
                              row.names=rownames(mela100.sigs.input))
mela.100.weights <- as.data.frame(t(sapply(rownames(mela.100.sigs.input), 
                                          function(x) whichSignatures(tumor.ref = mela.100.sigs.input,
                                                                      signatures.ref = signatures.cosmic,
                                                                      sample.id = x, 
                                                                      contexts.needed = TRUE,
                                                                      signature.cutoff = 0.00,              # default = 0.06
                                                                      tri.counts.method = 'default')$weights)), 
                                 row.names=rownames(mela.100.sigs.input))

mela50.weight <- as.data.frame(sapply(mela50.weights, function(x) unlist(x)), row.names = rownames(mela50.weights))
mela100.weight <- as.data.frame(sapply(mela100.weights, function(x) unlist(x)), row.names = rownames(mela100.weights))
mela.100.weight <- as.data.frame(sapply(mela.100.weights, function(x) unlist(x)), row.names = rownames(mela.100.weights))

# Change Signature.# to S#
colnames(mela50.weight) <- sapply(colnames(mela50.weight), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))
colnames(mela100.weight) <- sapply(colnames(mela100.weight), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))
colnames(mela.100.weight) <- sapply(colnames(mela.100.weight), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))

# Add unknown contribution column
mela50.weight$unk <- (1-rowSums(mela50.weight))
mela100.weight$unk <- (1-rowSums(mela100.weight))
mela.100.weight$unk <- (1-rowSums(mela.100.weight))



# Combine the two groups
mela.weight <- rbind(mela50.weight, mela.100.weight)   # 107 samples

save(mela.weight, file="./Output/ICGC_MELA_AU_weights_cut0.00.RData")        # 107 samples
save(mela50.weight, file="./Output/ICGC_MELA_AU_weights_50_cut0.00.RData")   # 50 samples
save(mela100.weight, file="./Output/ICGC_MELA_AU_100_weights_cut0.00.RData") # 43 samples
save(mela.100.weight, file="./Output/ICGC_MELA_AU_100_2_weights_cut0.00.RData") # 57 samples

# Sort by S1 (descending)
mela.sorted <- mela.weight[order(-mela.weight$S1),]
mela.sorted$icgc_sample_id <- rownames(mela.sorted)
mela.sorted <- merge(mela.sorted, mela.sample[,c(1,6)], by="icgc_sample_id", all.x=TRUE)
mela.sorted$icgc_sample_id <- mela.sorted$icgc_donor_id
colnames(mela.sorted)[1] <- "icgc_donor_id"
mela.sorted <- mela.sorted[,-33]

# Remove non-contributing signatures
mela.sorted <- as.matrix(mela.sorted)
mela.trim <- mela.sorted[,which(colSums(mela.sorted)>0)]
mela.melt <- melt(mela.sorted)
colnames(mela.melt) <- c("icgc_sample_id", "Signature", "weight")
mela.melt$weight <- sapply(mela.melt$weight, function(x) 100*x)

# Plot: boxplot signatures
pdf("./Figures/ICGC_MELA_AU_mutsigs_boxplot2.pdf", w=10, h=6)
ggplot(mela.melt, aes(x=Signature, y=weight, fill=Signature)) +
    geom_boxplot() +
    scale_fill_manual(values=colors2)+
    ggtitle("Mutational signature analysis ICGC MELA-AU (n=107)") +
    xlab("Signatures") +
    ylab("Contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,10)) +
    coord_cartesian(ylim=c(0,100), expand=FALSE) +
    theme_bw() +
    theme(legend.position="none",
          axis.title.x=element_blank())
dev.off()



# Anova
summary(mutsigs.anova <- aov(weight~Signature,  data=mela.melt))
pairwise <- as.data.frame(TukeyHSD(mutsigs.anova)[[1]])
pairwise$sign <- NA
colnames(pairwise)[4] <- "padj"
pairwise$sign <- ifelse(pairwise$padj < 0.05, "*","")

write.csv(pairwise, file="./Output/ICGC_MELA_mutsig_anova_pairwise.csv")

#--------- ~ Mean contributions per signature ---------

mela.means <- data.frame(Signatures=c(colnames(mela.trim)), value=c(colMeans(mela.trim)*100))
mela.fiveplus <- mela.means[which(mela.means$value>5),] 
others <- setdiff(rownames(mela.means), rownames(mela.fiveplus))
mela.others <- sum(mela.means$value[which(mela.means$Signatures %in% others)])
mela.five.others <- rbind(mela.fiveplus, c("oth", mela.others))

pdf("./Figures/ICGC_MELA_AU_mutsigs_avg.pdf", w=13, h=6)
ggplot(mela.means, aes(x=reorder(Signatures, -value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    ggtitle("Mean mutational signature contributions ICGC MELA-AU (n=93)") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,60), expand=FALSE) +
    theme_bw()
dev.off()

fiveplus <- ggplot(mela.five.others, aes(x=reorder(Signatures,-value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    #ggtitle("Mean mutational signature contributions ICGC MELA-AU (n=93)") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,60), expand=FALSE) +
    theme_bw()+
    scale_x_discrete(labels=c("S7","S2", "S8", "S5", "Others"))

pdf("./Figures/ICGC_MELA_AU_mutsigs_topavg.pdf", w=4, h=8)
fiveplus
dev.off()

save(mela.means, file="./Output/ICGC_MELA_AU_mutsigs_means.RData")


# Output boxplot and barplot as one plot
pdf("./Figures/TCGA_metscore_mel_stage_prog_boxplot.pdf", w=8, h=6)
ggarrange(sigbox, fiveplus, 
          #labels = c("a)", "b)", "c)", "", "", ""),
          ncol = 2, nrow = 1)
dev.off()


mela.means <- data.frame(Signatures=c(colnames(mela.sorted)), value=c(colMeans(mela.sorted)*100))
mela.means$Signatures <- factor(mela.means$Signatures, levels=order)

mela.means <- data.frame(Signatures=c(colnames(trim.tm)), value=c(colMeans(trim.tm)*100))
mela.means.2.5 <- mela.means[which(mela.means$value>2.5),]
mela.means.1 <- mela.means[which(mela.means$value>1),]
mela.means.oth <- setdiff(mela.means$Signatures, mela.means.1$Signatures)
sum.oth <- sum(mela.means$value[which(mela.means$Signatures %in% mela.means.oth)])
oth <- data.frame(Signatures="oth", value=sum.oth)
means.oth.mela <- rbind(mela.means.1, oth)

means.oth.tm$Signatures <- factor(means.oth.tm$Signatures, levels=order)
means.oth.tm$sample <- "MELA-AU"

colors10 <- brewer.pal(7, "Paired")
pdf("./Figures/TCGA_SKCM_tm_mutsigs_avg_2.pdf", w=10, h=6)
mela.sigs <- ggplot(means.oth.mela, aes(x=Signatures, y=value)) + 
    geom_bar(stat="identity", fill=colors10[7]) +
    #ggtitle("Mean mutational signature contributions SKCM primary tumor") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,60), expand=FALSE) +
    theme_bw() +
    theme(legend.position="none")
dev.off()

save(mela.means, file="./Output/TCGA_SKCM_TM_mutsigs_means.RData")

mela.mutsigs <- ggplot(means.oth.mela, aes(x=Signatures, y=value)) + 
    geom_bar(stat="identity", fill=colors10[7]) +
        #ggtitle("Mean mutational signature contributions SKCM primary tumor") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,2)) +
    coord_cartesian(ylim=c(0,20), expand=FALSE) +
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size=12),
          axis.title = element_text(size=13),
          legend.text = element_text(size=13))

mela.mutsigs1 <- ggplot(means.oth.mela, aes(x=Signatures, y=value)) + 
    geom_bar(stat="identity", fill=colors10[7]) +
    ggtitle("MELA-AU") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,2)) +
    coord_cartesian(ylim=c(55,60), expand=FALSE) +
    theme_bw() +
    theme(legend.position=c(0.87,0.7), 
          legend.title = element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=13))


# Get scka.mutsigs and skca.mutsigs1 from mutationAnalysis_ICGC_SKCA_BR.R
pdf("./Figures/ICGC_mutsigs_avg.pdf", w=10, h=5)
ggarrange(ggarrange(skca.mutsigs1, skca.mutsigs,
                    nrow=2,
                    heights = c(0.5,2)),
    ggarrange(mela.mutsigs1, mela.mutsigs,
                    nrow=2,
                    heights = c(0.5,2)),
    ncol=2)
dev.off()


#--------- TOTAL MUTATIONS CONTRIBUTION ------------------------------------------------------------------------------------

# Total number of mutations per sample
mela50.total.muts <- as.data.frame(table(mela50.snvs$icgc_donor_id))
mela.100.total.muts <- as.data.frame(table(mela.100.snvs$icgc_donor_id))
mela.total.muts <- rbind(mela50.total.muts,mela.100.total.muts)
colnames(mela.total.muts) <- c("icgc_donor_id","Total_Muts")

mela.total.muts <- merge(mela.metscore.stage, mela.total.muts, by="tcga_participant_barcode", all=TRUE)
mela.total.muts <- mela.total.muts[!is.na(mela.total.muts$patient_id),c(-2,-3,-37,-40,-41)]

# Replace signature % contribution columns with total mutation contributions (% sig contribution * total muts)
for (i in 2:32) {
    mela.total.muts[,i] <- mela.total.muts[,i]*mela.total.muts[,37]
}

save(mela.total.muts, file="./Output/TCGA_SKCM_TP_totalmutsigs.RData")

tm.mela.total.muts <- merge(tm.mela.metscore.stage, mela.total.muts, by="tcga_participant_barcode", all=TRUE)
tm.mela.total.muts <- tm.mela.total.muts[!is.na(tm.mela.total.muts$patient_id),c(-2,-3,-37,-40,-41)]

# Replace signature % contribution columns with total mutation contributions (% sig contribution * total muts)
for (i in 2:32) {
    tm.mela.total.muts[,i] <- tm.mela.total.muts[,i]*tm.mela.total.muts[,37]
}

save(tm.mela.total.muts, file="./Output/TCGA_SKCM_TM_totalmutsigs.RData")







#--------- SURVIVAL MODELS -----------------------------------------------------------------------------------------------

# Rename mutsigs table with donor id
temp <- mela.sample[,c(1,6)]
mela.mutsigs <- mela.weight
mela.mutsigs$icgc_sample_id <- rownames(mela.mutsigs)
mela.mutsigs <- merge(mela.mutsigs, temp, by="icgc_sample_id")
mela.mutsigs$icgc_sample_id <- mela.mutsigs$icgc_donor_id
colnames(mela.mutsigs)[1] <- "icgc_donor_id"
mela.mutsigs <- mela.mutsigs[, -33]

# Set up variable with all data
mela.clinall <- merge(mela.clin[,c(1,4:21)], mela.mutsigs, by="icgc_donor_id")
mela.clinall <- merge(mela.clin[,c(1,4:21)], mela.ther[,c(1,4:15)], by="icgc_donor_id")
mela.clinall <- merge(mela.clinall, mela.ther, by="icgc_donor_id", all.x=TRUE)

# Recode stage data
mela.clinall$stage <- as.character(mela.clinall$donor_tumour_stage_at_diagnosis)

# Combine IA/IB to I, IIA/IIB to II, IIIA/IIIB to III
# Set "" to NA, I/II to NA
for (i in 1:nrow(mela.clinall)) {
    if (mela.clinall$stage[i]=="IA/IB") {
        mela.clinall$stage[i] <- "I"
    }
    else if (mela.clinall$stage[i] == "IIA/IIB"){
        mela.clinall$stage[i] <- "II"
    }
    else if (mela.clinall$stage[i] == "IIIA/IIIB") {
        mela.clinall$stage[i] <- "III"
    }
    else if (mela.clinall$stage[i] == "IB/IIA" | mela.clinall$stage[i] == "") {
        mela.clinall$stage[i] <- NA
    }
    else {
        mela.clinall$stage[i] <- mela.clinall$stage[i]
    }
}
# Remove "A" "B" "C" stage specifications
for (i in 1:nrow(mela.clinall)) {
    if (grepl("a", mela.clinall$stage[i], ignore.case = TRUE) | grepl("b", mela.clinall$stage[i], ignore.case = TRUE) | 
        grepl("c", mela.clinall$stage[i], ignore.case = TRUE)) {
        mela.clinall$stage[i] <- substr(mela.clinall$stage[i], 1, nchar(mela.clinall$stage[i])-1)
    }
    else {
        mela.clinall$stage[i] <- mela.clinall$stage[i]
    }
}

# Turn vital status into binary
mela.clinall$donor_vital_status <- ifelse(mela.clinall$donor_vital_status == "deceased",1,0)
save(mela.clinall, file="./Output/ICGC_MELA_AU_clin_ther_mutsig.RData")

# Survival model
summary(mela.clinall.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~donor_sex+donor_age_at_diagnosis+stage+
                                      first_therapy_type+second_therapy_type+S7+S2+S8+S5+S3+S9+S1,
                                  data=mela.clinall))
summary(mela.trim.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~stage,
                                  data=mela.clinall))


summary(mela.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~stage+S5+donor_age_at_diagnosis,
                               data=mela.clinall))


summary(mela.mutsig.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~S1+S2+S3+S4+S5+S6+S7+S8+S9+S10+
                                     S11+S12+S13+S14+S15+S16+S17+S18+S19+S20+
                                     S21+S22+S23+S24+S25+S26+S27+S28+S29+S30+unk,
                               data=mela.clinall))
summary(mela.mutsig.trim.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~S6+S15,
                                 data=mela.clinall))


summary(mela.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~donor_age_at_diagnosis,
                                      data=mela.clinall))


mela.cox <- summary(mela.therapy.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~stage+S5+donor_age_at_diagnosis,
                                  data=mela.clinall))

g <- as.data.frame(mela.cox[[7]])
h <- as.data.frame(mela.cox[[8]])
mela.cox.mod <- cbind(rownames(g), g, h[,c(3:4)])
colnames(mela.cox.mod) <- c("var","coef", "HR", "SE", "z", "p", "CI_low", "CI_high")

melaau <- ggplot(mela.cox.mod, aes(x=factor(var, mela.cox.mod$var), y=coef)) + 
    scale_color_manual(values = c("(-Inf,0]" = "blue",
                                  "(0, Inf]" = "red")) +
    geom_errorbar(aes(ymin=coef-SE, ymax=coef+SE), width=0.05) +
    geom_point(aes(colour = cut(coef, c(-Inf, 0, Inf))),
               size = 3) +
    ggtitle("MELA-AU (n=94)") +
    ylab("log(HR)") + 
    scale_y_continuous(breaks=seq(-10,10,1)) +
    coord_cartesian(ylim=c(-1,5), expand=TRUE) + 
    theme_bw() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title = element_text(size=13),
          axis.text = element_text(size=12)) +
    scale_x_discrete(labels=c("Stage\nII", "Stage\nIII","Stage\nIV","S5", "Age")) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    annotate("text", x=1, y=-1, label=paste0("p=",round(mela.cox.mod$p[1], digits=3))) +
    annotate("text", x=2, y=-1, label=paste0("p=",round(mela.cox.mod$p[2], digits=3))) +
    annotate("text", x=3, y=-1, label=paste0("p=",round(mela.cox.mod$p[3], digits=3))) +
    annotate("text", x=4, y=-1, label=paste0("p=",round(mela.cox.mod$p[4], digits=3))) +
    annotate("text", x=5, y=-1, label=paste0("p=",round(mela.cox.mod$p[5], digits=3)))




# TCGA SKCM TP model
summary(mela.tcgatp.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~stage+S23+S11,
                                  data=mela.clinall))

# TCGA SKCM TM model
summary(mela.tcgatm.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~first_therapy_type+second_therapy_type+S7+S6,
                                  data=mela.clinall))

#--------- EXPOSURE VS SURVIVAL -------------------------------------------------------------------------
mela.clinex <- merge(mela.clinall, mela.expo[,c(1,5:9)], by="icgc_donor_id")

# Survival model (clinical, therapy, mutsigs, expo)
summary(mela.clinex.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~exposure_intensity+
                                     second_therapy_type,
                                 data=mela.clinex))







fitCoxbystage <- coxph(Surv(donor_survival_time, donor_vital_status)~stage,
                               data=mela.clinall)
predictions.stage <- predict(fitCoxbystage,
                        newdata = mela.clinall,
                        type="expected")

mela.clinall$probEventstage <- predictions.stage
boxplot(probEventstage~stage,data=mela.clinall)

