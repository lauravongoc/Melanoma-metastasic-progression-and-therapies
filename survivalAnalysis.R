# Survival Analysis
# DT Laura Vo Ngoc
# Start: 27/06/2018

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

# Stage
load("./Output/TCGA_SKCM_TP_stage.RData")       # SKCM TP
load("./Output/TCGA_SKCM_TM_stage.RData")       # SKCM TM

# Mutsigs, metscores, and stage
load("./Output/TCGA_SKCM_TP_metastatic_score_stage.RData")  # SKCM TP
load("./Output/TCGA_SKCM_TM_metastatic_score_stage.RData")  # SKCM TM

# Survival
load("./Output/TCGA_SKCM_TP_survival.RData")        # SKCM TP
load("./Output/TCGA_SKCM_TM_survival.RData")        # SKCM TM

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


#--------- SKCM COXPH -------------------------------------------------------------------------------------------------

#---- ~SKCM TP ----

# Remove not in tp.skcm.metscore.stage
a <- setdiff(tp.clin$tcga_participant_barcode, tp.skcm.metscore.stage$tcga_participant_barcode)
b <- as.vector(as.numeric(rownames(tp.clin[which(tp.clin$tcga_participant_barcode %in% a),])))
b <- b*(-1)

# Merge to for clin, stage, metscore, metpotential
tp.clin1 <- tp.clin[b,] 
tp.surv <- cbind(tp.clin1,tp.skcm.metscore.stage[,c(35,36,38,39)])

# Remove NA survival data (days_to_death) -- 147 removed
tp.surv <- tp.surv[!is.na(tp.surv$days_to_death),]

# Change vital status to binary, 1=dead, 0=alive (no NA)
tp.surv$vital_status <- ifelse(tp.surv$vital_status=="dead",1,0)

save(tp.surv, file="./Output/TCGA_SKCM_TP_survival.RData")

#---- ~SKCM TM ----

# Remove not in tm.skcm.metscore.stage
c <- setdiff(tm.clin$tcga_participant_barcode, tm.skcm.metscore.stage$tcga_participant_barcode)
d <- as.vector(as.numeric(rownames(tm.clin[which(tm.clin$tcga_participant_barcode %in% c),])))
d <- d*(-1)

# Merge to for clin, stage, metscore, metpotential
tm.clin1 <- tm.clin[d,] 
tm.surv <- cbind(tm.clin1,tm.skcm.metscore.stage[,c(35,36,38,39)])

# Remove NA survival data (days_to_death) -- 147 removed
tm.surv <- tm.surv[!is.na(tm.surv$days_to_death),]

# Generate column encompassing all treatment data
tm.surv$treatment <- tm.surv$prior_systemic_therapy_type

for (i in 1:nrow(tm.surv)) {
    if (!is.na(tm.surv$radiation_therapy[i]) && tm.surv$radiation_therapy[i]=="yes") {
        if (is.na(tm.surv$treatment[i])) {
            tm.surv$treatment[i] <- "radiation"
        }
        else {
            tm.surv$treatment[i] <- paste0(tm.surv$treatment[i], " + radiation")
        }
    }
}

# Change vital status to binary, 1=dead, 0=alive (no NA)
tm.surv$vital_status <- ifelse(tm.surv$vital_status=="dead",1,0)

# Change radiation therapy to binary, 1=yes, 0=no (no NA)
tm.surv$radiation_therapy <- ifelse(tm.surv$radiation_therapy=="yes",1,0)

# Change gender to binary, 1=male, 0=female (no NA)
tm.surv$gender <- ifelse(tm.surv$gender=="male",1,0)

# Create binary interferon column, 1=yes, 0=no (no NA)
tm.surv$radiation_therapy <- ifelse(tm.surv$radiation_therapy=="yes",1,0)

# Change met_potential to binary, 1=high, 0=low (no NA)
tm.surv$met_potential <- ifelse(tm.surv$met_potential=="high",1,0)

save(tm.surv, file="./Output/TCGA_SKCM_TM_survival.RData")


# ---- ~MODELS ---- 
# Outcome = overall survival = days_to_death

# Model TM -- Radiation only
tm.cox <- coxph(Surv(days_to_death, vital_status)~radiation_therapy,
                data=tm.surv)
summary(tm.cox)


# Model TM -- Treatments only
tm.cox <- coxph(Surv(days_to_death, vital_status)~treatment,
                data=tm.surv)
summary(tm.cox)


# Model -- Stage only
tp.cox <- coxph(Surv(days_to_death, vital_status)~stage,
                data=tp.surv)
summary(tp.cox)

tm.cox <- coxph(Surv(days_to_death, vital_status)~stage,
                data=tm.surv)
summary(tm.cox)


# Model -- Main mutsigs
tp.cox <- coxph(Surv(days_to_death, vital_status)~S7+S1+S11+S23+S6,
                data=tp.surv)
summary(tp.cox)

tm.cox <- coxph(Surv(days_to_death, vital_status)~S7+S1+S11+S23+S6,
                data=tm.surv)
summary(tm.cox)


# Model -- Mutsigs, one by one
tp.cox <- coxph(Surv(days_to_death, vital_status)~S7,
                data=tp.surv)
summary(tp.cox)
tp.cox <- coxph(Surv(days_to_death, vital_status)~S1,
                data=tp.surv)
summary(tp.cox)
tp.cox <- coxph(Surv(days_to_death, vital_status)~S11,
                data=tp.surv)
summary(tp.cox)
tp.cox <- coxph(Surv(days_to_death, vital_status)~S23,
                data=tp.surv)
summary(tp.cox)
tp.cox <- coxph(Surv(days_to_death, vital_status)~S6,
                data=tp.surv)
summary(tp.cox)

tm.cox <- coxph(Surv(days_to_death, vital_status)~S7,
                data=tm.surv)
summary(tm.cox)
tm.cox <- coxph(Surv(days_to_death, vital_status)~S1,
                data=tm.surv)
summary(tm.cox)
tm.cox <- coxph(Surv(days_to_death, vital_status)~S11,
                data=tm.surv)
summary(tm.cox)
tm.cox <- coxph(Surv(days_to_death, vital_status)~S23,
                data=tm.surv)
summary(tm.cox)
tm.cox <- coxph(Surv(days_to_death, vital_status)~S6,
                data=tm.surv)
summary(tm.cox)




# Model TM -- Treatment, main sigs, gender, stage, metscore
tm.cox <- coxph(Surv(days_to_death, vital_status)~treatment+
                    S7+S1+S11+S23+S6+
                    gender+
                    stage+earlyLate+
                    metscore+met_potential,
                data=tm.surv)
summary(tm.cox)
# Remove non-significant variables
tm.cox <- coxph(Surv(days_to_death, vital_status)~treatment+
                    S7+S6,
                data=tm.surv)
summary(tm.cox)


# Model TP -- Main sigs, gender, stage, metscore
tp.cox <- coxph(Surv(days_to_death, vital_status)~S7+S1+S11+S23+S6+
                    gender+
                    stage+earlyLate+
                    metscore+met_potential,
                data=tp.surv)
summary(tp.cox)
# Remove non-significant variables
tp.cox <- coxph(Surv(days_to_death, vital_status)~S11+S23+
                    stage,
                data=tp.surv)
summary(tp.cox)


# Model TM -- Treatment, all sigs
tm.cox <- coxph(Surv(days_to_death, vital_status)~treatment+
                    S1+S2+S3+S4+S5+S6+S7+S8+S9+S10+S11+S12+S13+S14+S15+S16+S17+S18+S19+S20+
                    S21+S22+S23+S24+S25+S26+S27+S28+S29+S30+unknown,
                    data=tm.surv)
summary(tm.cox)
# Remove sigs with NA results
tm.cox <- coxph(Surv(days_to_death, vital_status)~treatment+
                    S1+S2+S3+S4+S5+S6+S7+S10+S11+S12+S13+S16+S17+S18+S19+S20+
                    S21+S23+S26+S29,
                data=tm.surv)
summary(tm.cox)
# S19 is significant
tm.cox <- coxph(Surv(days_to_death, vital_status)~treatment+S29,
                data=tm.surv)
summary(tm.cox)



#--------- KAPLAN-MEIER: SURV VS METPOT ---------------------------------------------------------------------------------

#---- ~SKCM TP ----

survdiff(Surv(days_to_death, vital_status)~met_potential, data=tp.surv, rho=1)


# Removes samples with NA survival data (36 samples)
tp.fitKM <- survfit(Surv(days_to_death, vital_status)~met_potential,
                    data=tp.surv)
summary(tp.fitKM)

pdf("./Figures/TCGA_SKCM_TP_metpotential_survival_KM.pdf", w=8, h=6)
plot(tp.fitKM, col=c("red","blue"),
     mark.time = TRUE,
     main="SKCM TP overall survival by metastatic potential",
     xlab="Survival time (days)", ylab="% Overall survival")
legend("topright", legend = c("High metastatic potential", "Low metastatic potential"), fill=c("red","blue"))
dev.off()

tp.metpot.cox <- coxph(Surv(days_to_death, vital_status)~met_potential,
                data=tp.surv)
summary(tp.metpot.cox)


# Replace NA with largest survival value = 1070
tp.surv$surv_cens <- NA
tp.surv$surv_cens <- ifelse(is.na(tp.surv$days_to_death), 1070, tp.surv$days_to_death)

tp.cens.fitKM <- survfit(Surv(surv_cens, vital_status)~met_potential,
                    data=tp.surv)
summary(tp.cens.fitKM)

pdf("./Figures/TCGA_SKCM_TP_metpotential_survival_KM_cens.pdf", w=8, h=6)
plot(tp.cens.fitKM, col=c("red","blue"),
     mark.time = TRUE,
     main="SKCM TP overall survival by metastatic potential",
     xlab="Survival time (days)", ylab="% Overall survival")
legend("bottomright", legend = c("High metastatic potential", "Low metastatic potential"), fill=c("red","blue"))
dev.off()

survdiff(Surv(surv_cens, vital_status)~met_potential, data=tp.surv, rho=1)

tp.metpot.cox <- coxph(Surv(days_to_death, vital_status)~met_potential,
                       data=tp.surv)
summary(tp.metpot.cox)


# High vs. Low by quantiles
tp.surv$met_potential_quant <- NA
for (i in 1:nrow(tp.surv)) {
    if (tp.surv$metscore[i] < quantile(tp.surv$metscore)[2]) {
        tp.surv$met_potential_quant[i] <- "low"
    }
    else if (tp.surv$metscore[i] > quantile(tp.surv$metscore)[4]) {
        tp.surv$met_potential_quant[i] <- "high"
    }
    else {
        tp.surv$met_potential_quant[i] <- "intermediate"
    }
    
}

# Keep only high and low potentials by quantile
tp.surv.quant <- tp.surv[which(tp.surv$met_potential_quant=="low" | tp.surv$met_potential_quant=="high"),]


tp.quant.fitKM <- survfit(Surv(days_to_death, vital_status)~met_potential_quant,
                          data=tp.surv.quant)
summary(tp.quant.fitKM)

pdf("./Figures/TCGA_SKCM_TP_metpotential_survival_KM_quant.pdf", w=8, h=6)
plot(tp.quant.fitKM, col=c("red","blue"),
     mark.time = TRUE,
     main="SKCM TP overall survival by metastatic potential",
     xlab="Survival time (days)", ylab="% Overall survival")
legend("topright", legend = c("High metastatic potential", "Low metastatic potential"), fill=c("red","blue"))
dev.off()


tp.surv.quant$surv_cens <- NA
tp.surv.quant$surv_cens <- ifelse(is.na(tp.surv.quant$days_to_death), 1070, tp.surv.quant$days_to_death)

tp.quant.cens.fitKM <- survfit(Surv(surv_cens, vital_status)~met_potential_quant,
                          data=tp.surv.quant)
summary(tp.quant.cens.fitKM)

pdf("./Figures/TCGA_SKCM_TP_metpotential_survival_KM_quant_cens.pdf", w=8, h=6)
plot(tp.quant.cens.fitKM, col=c("red","blue"),
     mark.time = TRUE,
     main="SKCM TP overall survival by metastatic potential",
     xlab="Survival time (days)", ylab="% Overall survival")
legend("bottomright", legend = c("High metastatic potential", "Low metastatic potential"), fill=c("red","blue"))
dev.off()

survdiff(Surv(surv_cens, vital_status)~met_potential_quant, data=tp.surv.quant, rho=1)

tp.metpot.quant.cox <- coxph(Surv(days_to_death, vital_status)~met_potential_quant,
                       data=tp.surv.quant)
summary(tp.metpot.quant.cox)



#---- ~SKCM TM ----

survdiff(Surv(days_to_death, vital_status)~met_potential, data=tm.surv, rho=1)

tm.fitKM <- survfit(Surv(days_to_death, vital_status)~met_potential,
                 data=tm.surv)
summary(tm.fitKM)

pdf("./Figures/TCGA_SKCM_TM_metpotential_survival_KM.pdf", w=8, h=6)
plot(tm.fitKM, col=c("red","blue"),
     mark.time = TRUE,
     main="SKCM TM overall survival by metastatic potential",
     xlab="Survival time (days)", ylab="% Overall survival")
legend("topright", legend = c("High metastatic potential", "Low metastatic potential"), fill=c("red","blue"))
dev.off()

tm.metpot.cox <- coxph(Surv(days_to_death, vital_status)~met_potential,
                       data=tm.surv)
summary(tm.metpot.cox)