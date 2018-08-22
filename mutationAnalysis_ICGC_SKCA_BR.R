# DT Laura Vo Ngoc
# Start: 21/07/2018

library("BSgenome.Hsapiens.UCSC.hg19")
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

# Clinical data
load("./Data/ICGC_SKCA_BR_clinical.RData")

# Sample data
load("./Data/ICGC_SKCA_BR_sample.RData")

# Risk factor exposure data
load("./Data/ICGC_SKCA_BR_exposure.RData")

# Therapy data
load("./Data/ICGC_SKCA_BR_therapy.RData")

# SNV data (from mutations file)
load("./Output/ICGC_SKCA_BR_snvs.RData")

# Mutsig data
load("./Output/ICGC_SKCA_BR_weights_cut0.00.RData")     # Mutsig data
load("./Output/ICGC_SKCA_BR_mutsigs_means.RData")       # Mean mutsig contribution

# All combined data (clinical, therapy, mutsigs)
load("./Output/ICGC_SKCA_BR_clin_ther_mutsig.RData")



#---------  IMPORT DATA -----------------------------------------------------------------------------------------------

skca.clin <- read.table(file = './ICGC_SKCA_BR/donor.SKCA-BR.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
save(skca.clin, file="./Data/ICGC_SKCA_BR_clinical.RData")

skca.expo <- read.table(file = './ICGC_SKCA_BR/donor_exposure.SKCA-BR.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
save(skca.expo, file="./Data/ICGC_SKCA_BR_exposure.RData")

skca.ther <- read.table(file = './ICGC_SKCA_BR/donor_therapy.SKCA-BR.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
save(skca.ther, file="./Data/ICGC_SKCA_BR_therapy.RData")

skca.sample <- read.table(file = './ICGC_SKCA_BR/sample.SKCA-BR.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
save(skca.ther, file="./Data/ICGC_SKCA_BR_sample.RData")

skca.snvs <- read.table(file = './ICGC_SKCA_BR/simple_somatic_mutation.SKCA-BR.SB.shortcols.tsv', sep = '\t', header = FALSE)
colnames(skca.snvs) <- c("icgc_mutation_id", "icgc_donor_id", "project_code", "icgc_sample_id", "chromosome", "chromosome_start", "mutated_from_allele", "mutated_to_allele")
save(skca.snvs, file="./Output/ICGC_SKCA_BR_snvs.RData")


# Get colnames
mut.file <- "./ICGC_SKCA_BR/simple_somatic_mutation.open.SKCA-BR.tsv"
line <- readLines(mut.file, n=1)
line.split <- strsplit(line,"\\\t")
col.names <- line.split[[1]]


#--------- MUTATIONAL SIGNATURES ANALYSIS -------------------------------------------------------------------------

# Add "chr" to chromosome number in skca.snvs
skca.snvs$chromosome <- sapply(skca.snvs$chromosome, function(x) paste0("chr",x))

# Convert to deconstructSigs input:
skca.sigs.input <- mut.to.sigs.input(mut.ref = skca.snvs, 
                                     sample.id = "icgc_sample_id", 
                                     chr = "chromosome", 
                                     pos = "chromosome_start", 
                                     ref = "mutated_from_allele", 
                                     alt = "mutated_to_allele",
                                     bsg = BSgenome.Hsapiens.UCSC.hg19)

save(skca.sigs.input, file="./Output/ICGC_SKCA_BR_mutsig_input.RData")


# Signatures
skca.weights <- as.data.frame(t(sapply(rownames(skca.sigs.input), 
                                       function(x) whichSignatures(tumor.ref = skca.sigs.input,
                                                                   signatures.ref = signatures.cosmic,
                                                                   sample.id = x, 
                                                                   contexts.needed = TRUE,
                                                                   signature.cutoff = 0.00,              # default = 0.06
                                                                   tri.counts.method = 'default')$weights)), 
                              row.names=rownames(skca.sigs.input))

skca.weight <- as.data.frame(sapply(skca.weights, function(x) unlist(x)), row.names = rownames(skca.weights))

# Change Signature.# to S#
colnames(skca.weight) <- sapply(colnames(skca.weight), function(x) paste0("S",strsplit(x, "\\.")[[1]][2]))

# Add unknown contribution column
skca.weight$unk <- (1-rowSums(skca.weight))

save(skca.weight, file="./Output/ICGC_SKCA_BR_weights_cut0.00.RData")

# Sort by S1 (descending)
skca.sorted <- skca.weight[order(-skca.weight$S1),]

# Remove non-contributing signatures
skca.sorted <- as.matrix(skca.sorted)
skca.trim <- skca.sorted[,which(colSums(skca.sorted)>0)]
skca.melt <- melt(skca.sorted)
colnames(skca.melt) <- c("sample_id", "Signature", "weight")
skca.melt$weight <- sapply(skca.melt$weight, function(x) 100*x)

# Plot: boxplot signatures
pdf("./Figures/ICGC_SKCA_BR_mutsigs_boxplot.pdf", w=10, h=6)
ggplot(skca.melt, aes(x=Signature, y=weight, fill=Signature)) +
    geom_boxplot() +
    scale_fill_manual(values=colors2)+
    ggtitle("Mutational signature analysis ICGC SKCA-BR") +
    xlab("Signatures") +
    ylab("Contribution") +
    coord_cartesian(ylim=c(0,100), expand=FALSE)
dev.off()

# Anova
summary(mutsigs.anova <- aov(weight~Signature,  data=skca.melt))
pairwise <- as.data.frame(TukeyHSD(mutsigs.anova)[[1]])
pairwise$sign <- NA
colnames(pairwise)[4] <- "padj"
pairwise$sign <- ifelse(pairwise$padj < 0.05, "*","")

write.csv(pairwise, file="./Output/ICGC_SKCA_BR_mutsig_anova_pairwise.csv")

#--------- ~ Mean contributions per signature ---------

skca.means <- data.frame(Signatures=c(colnames(skca.trim)), value=c(colMeans(skca.trim)*100))

pdf("./Figures/ICGC_SKCA_BR_mutsigs_avg.pdf", w=13, h=6)
ggplot(skca.means, aes(x=reorder(Signatures, -value), y=value)) + 
    geom_bar(stat="identity", fill="dark blue") +
    ggtitle("Mean mutational signature contributions ICGC SKCA-BR") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,5)) +
    coord_cartesian(ylim=c(0,40), expand=FALSE)
dev.off()

save(skca.means, file="./Output/ICGC_SKCA_BR_mutsigs_means.RData")


skca.means <- data.frame(Signatures=c(colnames(skca.sorted)), value=c(colMeans(skca.sorted)*100))
skca.means$Signatures <- factor(skca.means$Signatures, levels=order)

skca.means <- data.frame(Signatures=c(colnames(trim.tm)), value=c(colMeans(trim.tm)*100))
skca.means.2.5 <- skca.means[which(skca.means$value>2.5),]
skca.means.1 <- skca.means[which(skca.means$value>1),]
skca.means.oth <- setdiff(skca.means$Signatures, skca.means.1$Signatures)
sum.oth <- sum(skca.means$value[which(skca.means$Signatures %in% skca.means.oth)])
oth <- data.frame(Signatures="oth", value=sum.oth)
means.oth.skca <- rbind(skca.means.1, oth)

means.oth.skca$Signatures <- factor(means.oth.tm$Signatures, levels=order)
means.oth.tm$sample <- "SKCA-BR"

skca.mutsigs <- ggplot(means.oth.skca, aes(x=Signatures, y=value)) + 
    geom_bar(stat="identity", fill=colors10[5]) +
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

skca.mutsigs1 <- ggplot(means.oth.skca, aes(x=Signatures, y=value)) + 
    geom_bar(stat="identity", fill=colors10[5]) +
    ggtitle("SKCA-BR") +
    xlab("Signatures") +
    ylab("Mean contribution (%)") +
    scale_y_continuous(breaks=seq(0,100,2)) +
    coord_cartesian(ylim=c(35,40), expand=FALSE) +
    theme_bw() +
    theme(legend.position=c(0.87,0.7), 
          legend.title = element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=13))

save(skca.means, file="./Output/TCGA_SKCM_TM_mutsigs_means.RData")


#--------- SURVIVAL ANALYSIS -------------------------------------------------------------------------

# Set up dataframe with all data
temp <- skca.sample[,c(1,6)]
skca.mutsigs <- skca.weight
skca.mutsigs$icgc_sample_id <- rownames(skca.mutsigs)
skca.mutsigs <- merge(skca.mutsigs, temp, by="icgc_sample_id")
skca.mutsigs$icgc_sample_id <- skca.mutsigs$icgc_donor_id
colnames(skca.mutsigs)[1] <- "icgc_donor_id"
skca.mutsigs <- skca.mutsigs[, -33]

skca.clinall <- merge(skca.clin[,c(1,4:21)], skca.ther[,c(1,4:15)], by="icgc_donor_id")
skca.clinall <- merge(skca.clinall, skca.mutsigs, by="icgc_donor_id")

# Recode stage data
skca.clinall$stage <- skca.clinall$donor_tumour_stage_at_diagnosis

for (i in 1:nrow(skca.clinall)) {
    if (skca.clinall$stage[i] == "1") {
        skca.clinall$stage[i] <- "I"
    }
    else if (skca.clinall$stage[i] == "2"){
        skca.clinall$stage[i] <- "II"
    }
    else if (skca.clinall$stage[i] == "3") {
        skca.clinall$stage[i] <- "III"
    }
    else if (skca.clinall$stage[i] == "4") {
        skca.clinall$stage[i] <- "IV"
    }
}
for (i in 1:nrow(skca.clinall)) {
    if (grepl("I", skca.clinall$stage[i], ignore.case = TRUE)) {
        skca.clinall$stage[i] <- skca.clinall$stage[i]
    }
    else {
        skca.clinall$stage[i] <- NA
    }
}

# Recode vital status into binary
skca.clinall$donor_vital_status <- ifelse(skca.clinall$donor_vital_status == "deceased", 1, 0)
save(skca.clinall, file="./Output/ICGC_SKCA_BR_clin_ther_mutsig.RData")

# Survival model (clinical, therapy, mutsigs)
summary(skca.clinall.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~donor_sex+donor_age_at_diagnosis+stage+
                              first_therapy_type+second_therapy_type+S7+S3+S12+S5+S9+S1,
                          data=skca.clinall))

summary(skca.trimmed.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~stage+
                              first_therapy_type+second_therapy_type,
                          data=skca.clinall))

summary(skca.fin.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~stage+
                              second_therapy_type,
                          data=skca.clinall))


# Survival model (mutsigs)
summary(skca.mutsigs.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~S1+S2+S3+S4+S5+S6+S7+S8+S9+S10+
                                      S11+S12+S13+S14+S15+S16+S17+S18+S19+S20+
                                      S21+S22+S23+S24+S25+S26+S27+S28+S29+S30+unk,
                              data=skca.clinall))

summary(skca.mutsig.trim.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~S1+S6+S10+
                                      S13+S15+S18+S21+S28,
                                  data=skca.clinall))

summary(skca.mutsig.fin.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~S6+S10+S13,
                                      data=skca.clinall))

summary(coxph(Surv(donor_survival_time, donor_vital_status)~second_therapy_type+stage,
                                  data=skca.clinall))

skca.cox <- summary(coxph(Surv(donor_survival_time, donor_vital_status)~stage+S6+second_therapy_type,
              data=skca.clinall))

e <- as.data.frame(skca.cox[[7]])
f <- as.data.frame(skca.cox[[8]])
skca.cox.mod <- cbind(rownames(e), e, f[,c(3:4)])
colnames(skca.cox.mod) <- c("var","coef", "HR", "SE", "z", "p", "CI_low", "CI_high")

skca <- ggplot(skca.cox.mod, aes(x=factor(var, skca.cox.mod$var), y=coef)) + 
    scale_color_manual(values = c("(-Inf,0]" = "blue",
                                  "(0, Inf]" = "red")) +
    geom_errorbar(aes(ymin=coef-SE, ymax=coef+SE), width=0.05) +
    geom_point(aes(colour = cut(coef, c(-Inf, 0, Inf))),
               size = 3) +
    ggtitle("SKCA-BR (n=95)") +
    ylab("log(HR)") + 
    scale_y_continuous(breaks=seq(-10,50,10)) +
    coord_cartesian(ylim=c(-5,50), expand=TRUE) + 
    theme_bw() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title = element_text(size=13),
          axis.text = element_text(size=12)) +
    scale_x_discrete(labels=c("Stage\nII", "Stage\nIII", "Stage\nIV", "S6",
        "Immuno-\ntherapy", "No\ntreatment","Radiation", "Surgery")) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    annotate("text", x=1, y=-5, label=paste0("p=",round(skca.cox.mod$p[1], digits=3))) +
    annotate("text", x=2, y=-5, label=paste0("p=",round(skca.cox.mod$p[2], digits=3))) +
    annotate("text", x=3, y=-5, label=paste0("p=",round(skca.cox.mod$p[3], digits=3))) +
    annotate("text", x=4, y=-5, label=paste0("p=",round(skca.cox.mod$p[4], digits=3))) +
    annotate("text", x=5, y=-5, label=paste0("p=",round(skca.cox.mod$p[5], digits=3))) +
    annotate("text", x=6, y=-5, label=paste0("p=",round(skca.cox.mod$p[6], digits=3))) +
    annotate("text", x=7, y=-5, label=paste0("p=",round(skca.cox.mod$p[7], digits=3))) +
    annotate("text", x=8, y=-5, label=paste0("p=",round(skca.cox.mod$p[8], digits=3)))

# KM plot stage
skca.stage.fitKM <- survfit(Surv(donor_survival_time, donor_vital_status)~stage,
                            data=skca.clinall)
summary(tp.fitKM)

skca.stage.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~stage,
                       data=skca.clinall)
cox <- summary(skca.stage.cox)

hr <- round(cox$coefficients[2], digits=2)
CI1 <- round(cox$conf.int[3], digits=2)
CI2 <- round(cox$conf.int[4], digits=2)
pval <- round(cox$sctest[3], digits=4)

pdf("./Figures/ICGC_SKCA_BR_survival_stage_KM_annot.pdf", w=8, h=6)
plot(skca.stage.fitKM, col=col[2:5],
     mark.time = TRUE,
     main="SKCA-BR overall survival by cancer stage",
     xlab="Survival time (days)", ylab="% Overall survival")
legend("topright", 
       legend = c("Stage I","Stage II", "Stage III", "Stage IV"), 
       fill=col[2:5], 
       bty="n")
#y.intersp=1, x.intersp=2, text.width=0.8)
legend("bottomleft", 
       legend = c("HR: ", paste0(hr, ", 95% CI: ", CI1, "-", CI2), paste0("p = ", pval)), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1)
dev.off()




#--------- UV EXPOSURE ANALYSIS -------------------------------------------------------------------------

skca.sun <- merge(skca.expo[,c(1,5)], skca.mutsigs[,c(1,8)], by="icgc_donor_id")
skca.sun.melt <- melt(skca.sun)
skca.sun$S7 <- skca.sun$S7*100
skca.sun$intensity <- skca.sun$exposure_intensity
skca.sun$intensity[which(skca.sun$intensity==0)] <- "0"
skca.sun$intensity[which(skca.sun$intensity==1)] <- "1"
skca.sun$intensity[which(skca.sun$intensity==2)] <- "2"

skca.sun <- skca.sun[!is.na(skca.sun$intensity),]
skca.sun$intensity <- factor(skca.sun$intensity, levels=c("None","Low","High"))

skca.sun.trim <- skca.sun[which(skca.sun$S7 > 0.8 & skca.sun$S7 < quantile(skca.sun$S7)[4][[1]]),]
skca.sun.trim <- skca.sun.trim[!is.na(skca.sun.trim$intensity),]
skca.sun.trim$intensity <- factor(skca.sun.trim$intensity, levels=c("None","Low","High"))

# Boxplot by sun exposure
col <- brewer.pal(5, "Blues")
pdf("./Figures/ICGC_SKCA_BR_expo_s7.pdf", w=8, h=6)
ggplot(skca.sun, aes(x=intensity, y=S7, fill=intensity)) +
    geom_boxplot() +
    scale_fill_manual(values=col[2:4]) +
    theme(legend.position = "none") +
    ggtitle("S7 contribution by sun exposure ICGC SKCA-BR") +
    xlab("Exposure intensity") +
    ylab("S7 contribution") +
    coord_cartesian(ylim=c(0,100), expand=FALSE)
dev.off()

# Anova
summary(sun.anova <- aov(S7~intensity,  data=skca.sun))
TukeyHSD(sun.anova)[[1]]


summary(coxph(Surv(donor_survival_time, donor_vital_status)~first_therapy_type,
                                  data=skca.clinall))



#--------- EXPOSURE VS SURVIVAL -------------------------------------------------------------------------
skca.clinex <- merge(skca.clinall, skca.expo[,c(1,5:9)], by="icgc_donor_id")

# Survival model (clinical, therapy, mutsigs, expo)
summary(skca.clinex.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~tobacco_smoking_intensity+
                                     alcohol_history+alcohol_history_intensity+second_therapy_type,
                                  data=skca.clinex))

summary(skca.clinex.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~tobacco_smoking_intensity+
                                     alcohol_history+alcohol_history_intensity,
                                 data=skca.clinex))


# KM plot sun
skca.sun.fitKM <- survfit(Surv(donor_survival_time, donor_vital_status)~exposure_intensity,
                            data=skca.clinex)
summary(tp.fitKM)

skca.sun.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~exposure_intensity,
                        data=skca.clinex)
cox <- summary(skca.sun.cox)

hr <- round(cox$coefficients[2], digits=2)
CI1 <- round(cox$conf.int[3], digits=2)
CI2 <- round(cox$conf.int[4], digits=2)
pval <- round(cox$sctest[3], digits=4)


col1 <- brewer.pal(9, "Blues")
pdf("./Figures/ICGC_SKCA_BR_survival_sun_KM_annot.pdf", w=8, h=6)
plot(skca.sun.fitKM, col=col1[c(5,7,9)],
     mark.time = TRUE,
     main="SKCA-BR overall survival by sun exposure intensity",
     xlab="Survival time (days)", ylab="% Overall survival")
legend("topright", 
       legend = c(0, 1, 2), 
       fill=col1[c(5,7,9)], 
       bty="n")
#y.intersp=1, x.intersp=2, text.width=0.8)
legend("bottomleft", 
       legend = c("HR: ", paste0(hr, ", 95% CI: ", CI1, "-", CI2), paste0("p = ", pval)), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1)
dev.off()



fitCoxbystage <- coxph(Surv(donor_survival_time, donor_vital_status)~stage,
                       data=skca.clinall)
predictions.stage <- predict(fitCoxbystage,
                             newdata = skca.clinall,
                             type="expected")

skca.clinall$probEventstage <- predictions.stage
boxplot(probEventstage~stage,data=skca.clinall)





skca.clinall2 <- merge(skca.clin[,c(1,4:21)], skca.ther[,c(1,4:15)], by="icgc_donor_id")
# Recode stage data
skca.clinall2$stage <- skca.clinall$donor_tumour_stage_at_diagnosis

for (i in 1:nrow(skca.clinall2)) {
    if (skca.clinall2$stage[i] == "1") {
        skca.clinall2$stage[i] <- "I"
    }
    else if (skca.clinall2$stage[i] == "2"){
        skca.clinall2$stage[i] <- "II"
    }
    else if (skca.clinall2$stage[i] == "3") {
        skca.clinall2$stage[i] <- "III"
    }
    else if (skca.clinall2$stage[i] == "4") {
        skca.clinall$stage[i] <- "IV"
    }
}
for (i in 1:nrow(skca.clinall2)) {
    if (grepl("I", skca.clinall2$stage[i], ignore.case = TRUE)) {
        skca.clinall2$stage[i] <- skca.clinall2$stage[i]
    }
    else {
        skca.clinall2$stage[i] <- NA
    }
}

# Recode vital status into binary
skca.clinall2$donor_vital_status <- ifelse(skca.clinall2$donor_vital_status == "deceased", 1, 0)

# Survival model (clinical, therapy, mutsigs)
summary(skca.clinall2.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~stage+
                                      first_therapy_type,
                                  data=skca.clinall2))