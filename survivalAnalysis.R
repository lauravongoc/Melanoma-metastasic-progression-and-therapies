# Survival Analysis
# DT Laura Vo Ngoc
# Start: 27/06/2018

library(deconstructSigs)
library(ggplot2)
library(ggpubr)
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


# ICGC data
load("./Output/ICGC_SKCA_BR_clin_ther_mutsig.RData")    # SKCA-BR
load("./Output/ICGC_MELA_AU_clin_ther_mutsig.RData")    # MELA-AU

# All cutaneous melanoma data survival
load("./Output/ICGC_TCGA_survival.RData")

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


#---- ~UVM ----

uvm.surv <- merge(clin, uvm.metscore.stage[,c(1,35,36,38,39)], by="tcga_participant_barcode")

# Remove NA survival data (days_to_death) -- 57 removed
uvm.surv <- uvm.surv[!is.na(uvm.surv$days_to_death),]

# Change vital status to binary, 1=dead, 0=alive (no NA)
uvm.surv$vital_status <- ifelse(uvm.surv$vital_status=="dead",1,0)

save(uvm.surv, file="./Output/TCGA_UVM_TP_survival.RData")




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


# Model TM -- Treatment, main sigs, gender, stage, metscore
tm.cox <- coxph(Surv(days_to_death, vital_status)~treatment+
                    S7+S1+S11+S23+S19+
                    gender+
                    stage+earlyLate+
                    metscore+met_potential,
                data=tm.surv)
summary(tm.cox)



#### ... tm final ####
tm.cox <- summary(tm.c <- coxph(Surv(days_to_death, vital_status)~treatment+
                    S7+S6,
                data=tm.surv))
summary(tm.cox)


# Model TP -- Main sigs, gender, stage, metscore
tp.cox <- coxph(Surv(days_to_death, vital_status)~S7+S1+S11+S23+S6+
                    gender+
                    stage+earlyLate+
                    metscore+met_potential,
                data=tp.surv)
summary(tp.cox)



#### ... tp final ####
tp.cox <- summary(tp.c <- coxph(Surv(days_to_death, vital_status)~stage+S11+S23,
                data=tp.surv))





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
tm.cox <- summary(coxph(Surv(days_to_death, vital_status)~treatment+S29,
                data=tm.surv))
summary(tm.cox)


#### ooooooooo ####

#--------- ALIVE VS DEAD SIGS ---------------------------------------------------------------------------------

# Variables with only patients annotated with survival data (and therapy data)
tp.alive.dead <- tp.surv[!is.na(tp.surv$days_to_death),]        # 27 patients
tm.alive.dead <- tm.surv[!is.na(tm.surv$days_to_death),]        # 162 patients
tm.alive.dead <- tm.alive.dead[!is.na(tm.alive.dead$treatment),]        # 21 patients

# Change vital status back to "alive" and "deceased
tp.alive.dead$vital_status <- ifelse(tp.alive.dead$vital_status==1, "Deceased", "Alive")
tm.alive.dead$vital_status <- ifelse(tm.alive.dead$vital_status==1, "Deceased", "Alive")

# Plot SKCM TP S11 alive vs. dead
ggplot(tp.alive.dead, aes(x=vital_status, y=S11, fill=vital_status)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("") +
    xlab("Vital status") +
    ylab("Signature contribution (%)") +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    #coord_cartesian( ylim=c(0.4,1.8), expand = TRUE) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=13))
    #scale_x_discrete(labels=c("Alive\nn=0","Deceased\nn=27"))

# Plot SKCM TP S23 alive vs. dead
ggplot(tp.alive.dead, aes(x=vital_status, y=S23, fill=vital_status)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("") +
    xlab("Vital status") +
    ylab("Signature contribution (%)") +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    #coord_cartesian( ylim=c(0.4,1.8), expand = TRUE) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=13))
#scale_x_discrete(labels=c("Alive\nn=0","Deceased\nn=27"))

# Plot SKCM TP S6 alive vs. dead
ggplot(tm.alive.dead, aes(x=vital_status, y=S6, fill=vital_status)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("") +
    xlab("Vital status") +
    ylab("Signature contribution (%)") +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    #coord_cartesian( ylim=c(0.4,1.8), expand = TRUE) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=13))
#scale_x_discrete(labels=c("Alive\nn=0","Deceased\nn=27"))

# Plot SKCM TP S7 alive vs. dead
ggplot(tp.alive.dead, aes(x=vital_status, y=S7, fill=vital_status)) + 
    scale_fill_manual(values=colors4) +
    geom_boxplot() +
    ggtitle("") +
    xlab("Vital status") +
    ylab("Signature contribution (%)") +
    scale_y_continuous(breaks=seq(0,10,0.2)) +
    #coord_cartesian( ylim=c(0.4,1.8), expand = TRUE) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=13))
#scale_x_discrete(labels=c("Alive\nn=0","Deceased\nn=27"))


#--------- KAPLAN-MEIER: SURV VS METPOT ---------------------------------------------------------------------------------

#---- ~SKCM TP ----

survdiff(Surv(days_to_death, vital_status)~met_potential, data=tp.surv, rho=0)


# Survival by metpotential
tp.fitKM <- survfit(Surv(days_to_death, vital_status)~met_potential,
                    data=tp.surv)
summary(tp.fitKM)

tp.metpotential.cox <- coxph(Surv(days_to_death, vital_status)~met_potential,
                data=tp.surv)
cox <- summary(tp.metpotential.cox)

hr <- round(cox$coefficients[2], digits=2)
CI1 <- round(cox$conf.int[3], digits=2)
CI2 <- round(cox$conf.int[4], digits=2)
pval <- round(cox$sctest[3], digits=4)

pdf("./Figures/TCGA_SKCM_TP_metpotential_survival_KM_annot.pdf", w=7, h=6)
plot(tp.fitKM, col=c("red","blue"),
     mark.time = TRUE,
     xlab="Survival time (days)", ylab="% Overall survival",
     lwd = 3,
     cex.axis = 1.3,
     cex.lab = 1.5)
legend("topright", 
       legend = c("High potential (n=6)", "Low potential (n=21)"), 
       fill=c("red","blue"), 
       bty="n",
       cex = 1.5)
legend("bottomleft", 
       legend = c("HR: ", paste0(hr, ", 95% CI: ", CI1, "-", CI2), paste0("p = ", pval)), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1,
       cex = 1.5)
dev.off()



# Survival by stage
tp.stage.fitKM <- survfit(Surv(days_to_death, vital_status)~stage,
                    data=tp.surv)
summary(tp.stage.fitKM)

tp.stage.cox <- coxph(Surv(days_to_death, vital_status)~stage,
                       data=tp.surv)
cox <- summary(tp.stage.cox)

hr <- round(cox$coefficients[2], digits=2)
CI1 <- round(cox$conf.int[3], digits=2)
CI2 <- round(cox$conf.int[4], digits=2)
pval <- round(cox$sctest[3][[1]], digits=7)


# All stage I and IV were NA in days_to_death
pdf("./Figures/TCGA_SKCM_TP_stage_survival_KM_annot.pdf", w=7, h=6)
plot(tp.stage.fitKM, col=c("blue","forest green"),
     mark.time = TRUE,
     xlab="Survival time (days)", ylab="% Overall survival",
     lwd = 3,
     cex.axis = 1.3,
     cex.lab = 1.5)
legend("topright", 
       legend = c("Stage II (n=21)", "Stage III (n=6)"), 
       fill=c("blue","forest green"), 
       bty="n",
       cex = 1.5)
#y.intersp=1, x.intersp=2, text.width=0.8)
legend("bottomleft", 
       legend = c("HR: ", paste0(hr, ", 95% CI: ", CI1, "-", CI2), paste0("p = ", pval)), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1,
       cex=1.5)
dev.off()


# Survival by gender
tp.fitKM <- survfit(Surv(days_to_death, vital_status)~gender,
                    data=tp.surv)
summary(tp.fitKM)

tp.gender.cox <- coxph(Surv(days_to_death, vital_status)~gender,
                             data=tp.surv)
cox <- summary(tp.gender.cox)

hr <- round(cox$coefficients[2], digits=2)
CI1 <- round(cox$conf.int[3], digits=2)
CI2 <- round(cox$conf.int[4], digits=2)
pval <- round(cox$sctest[3], digits=4)

pdf("./Figures/TCGA_SKCM_TP_metpotential_survival_KM_annot.pdf", w=8, h=6)
plot(tp.fitKM, col=c("red","blue"),
     mark.time = TRUE,
     main="SKCM TP overall survival by gender",
     xlab="Survival time (days)", ylab="% Overall survival")
legend("topright", 
       legend = c("Female (n=11)", "Male (n=16)"), 
       fill=c("red","blue"), 
       bty="n")
#y.intersp=1, x.intersp=2, text.width=0.8)
legend("bottomleft", 
       legend = c("HR: ", paste0(hr, ", 95% CI: ", CI1, "-", CI2), paste0("p = ", pval)), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1)
dev.off()

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


# Survival by stage
tm.stage.fitKM <- survfit(Surv(days_to_death, vital_status)~stage,
                          data=tm.surv)
summary(tm.stage.fitKM)

tm.stage.cox <- coxph(Surv(days_to_death, vital_status)~stage,
                      data=tm.surv)
cox <- summary(tm.stage.cox)

hr <- round(cox$coefficients[2], digits=2)
CI1 <- round(cox$conf.int[3], digits=2)
CI2 <- round(cox$conf.int[4], digits=2)
pval <- round(cox$sctest[3][[1]], digits=4)

pdf("./Figures/TCGA_SKCM_TM_stage_survival_KM_annot.pdf", w=7, h=6)
plot(tm.stage.fitKM, col=c("red","blue","forest green","orange"),
     mark.time = TRUE,
     xlab="Survival time (days)", ylab="% Overall survival",
     lwd = 3,
     cex.axis = 1.3,
     cex.lab =1.5)
legend("topright", 
       legend = c("Stage I (n=37)","Stage II (n=36)", "Stage III (n=77)", "Stage IV (n=12)"), 
       fill=c("red","blue","forest green","orange"), 
       bty="n",
       cex = 1.5)
#y.intersp=1, x.intersp=2, text.width=0.8)
legend("bottomleft", 
       legend = c(paste0("p = ", pval)), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1,
       cex = 1.5)
dev.off()


# Survival by treatment, radiation vs others
tm.rad <- tm.surv

for (i in 1:nrow(tm.rad)) {
    if (!is.na(tm.rad$treatment[i])) {
        if (tm.rad$treatment[i]!="radiation") {
            tm.rad$treatment[i] <- "other"
        }
    }
}
    
tm.rad.fitKM <- survfit(Surv(days_to_death, vital_status)~treatment,
                          data=tm.rad)    
summary(tm.ther.fitKM)

tm.rad.cox <- coxph(Surv(days_to_death, vital_status)~treatment,
                      data=tm.rad)
cox <- summary(tm.rad.cox)

hr <- round(cox$coefficients[2], digits=2)
CI1 <- round(cox$conf.int[3], digits=2)
CI2 <- round(cox$conf.int[4], digits=2)
pval <- round(cox$sctest[3][[1]], digits=4)

pdf("./Figures/TCGA_SKCM_TM_rad_survival_KM_annot.pdf", w=7, h=6)
plot(tm.rad.fitKM, col=c("blue","red"),
     mark.time = TRUE,
     xlab="Survival time (days)", ylab="% Overall survival",
     lwd = 3,
     cex.axis = 1.3,
     cex.lab = 1.5)
legend("topright", 
       legend = c("radiation (n=7)","other (n=14)"), 
       fill=c("red","blue"), 
       bty="n",
       cex=1.5)
#y.intersp=1, x.intersp=2, text.width=0.8)
legend("bottomleft", 
       legend = c("HR: ", paste0(hr, ", 95% CI: ", CI1, "-", CI2), paste0("p = ", pval)), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1,
       cex=1.5)
dev.off()



# Model censored data
tp.cens.cox <- coxph(Surv(surv_cens, vital_status)~S1+S7+S23+S11+S6+stage+earlyLate+metscore+met_potential,
                         data=tp.surv)
summary(tp.cens.cox)


tp.cens.cox <- coxph(Surv(surv_cens, vital_status)~S23+met_potential,
                     data=tp.surv)
summary(tp.cens.cox)



#---- ~UVM ----

survdiff(Surv(days_to_death, vital_status)~met_potential, data=uvm.surv, rho=1)

uvm.fitKM <- survfit(Surv(days_to_death, vital_status)~met_potential,
                    data=uvm.surv)
summary(uvm.fitKM)

cox <- summary(uvm.metpot.cox <- coxph(Surv(days_to_death, vital_status)~met_potential,
                                data=uvm.surv))

hr <- round(cox$coefficients[2], digits=2)
CI1 <- round(cox$conf.int[3], digits=2)
CI2 <- round(cox$conf.int[4], digits=2)
pval <- round(cox$sctest[3], digits=4)

pdf("./Figures/TCGA_UVM_TP_metpotential_survival_KM.pdf", w=7, h=6)
plot(uvm.fitKM, col=c("red","blue"),
     mark.time = TRUE,
     xlab="Survival time (days)", ylab="% Overall survival",
     lwd = 3,
     cex.axis = 1.3,
     cex.lab = 1.5)
legend("topright", 
       legend = c("High potential (n=19)", "Low potential (n=4)"), 
       fill = c("red","blue"),
       bty = "n",
       cex = 1.5)
legend("bottomleft", 
       legend = c("HR: ", paste0(hr, ", 95% CI: ", CI1, "-", CI2), paste0("p = ", pval)), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1,
       cex = 1.5)
dev.off()




# Censoring

# Replace NA with largest survival value = 1581
uvm.surv$surv_cens <- NA
uvm.surv$surv_cens <- ifelse(is.na(uvm.surv$days_to_death), 1581, uvm.surv$days_to_death)


survdiff(Surv(surv_cens, vital_status)~met_potential, data=uvm.surv, rho=1)

uvm.cens.fitKM <- survfit(Surv(surv_cens, vital_status)~met_potential,
                     data=uvm.surv)
summary(uvm.fitKM)

cox.cens <- summary(uvm.metpot.cens.cox <- coxph(Surv(surv_cens, vital_status)~met_potential,
                                     data=uvm.surv))

hr <- round(cox.cens$coefficients[2], digits=2)
CI1 <- round(cox.cens$conf.int[3], digits=2)
CI2 <- round(cox.cens$conf.int[4], digits=2)
pval <- round(cox.cens$sctest[3], digits=9)

pdf("./Figures/TCGA_UVM_TP_metpotential_survival_KM_censored.pdf", w=8, h=6)
plot(uvm.cens.fitKM, col=c("red","blue"),
     mark.time = TRUE,
     main="UVM TP censored overall survival by metastatic potential",
     xlab="Survival time (days)", ylab="% Overall survival")
legend("bottomright", 
       legend = c("High metastatic potential (n=38)", "Low metastatic potential (n=42)"), 
       fill = c("red","blue"),
       bty = "n")
legend("bottomleft", 
       legend = c("HR: ", paste0(hr, ", 95% CI: ", CI1, "-", CI2), paste0("p = ", pval)), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1)
dev.off()



#--------- INDIV MODEL PLOTS  ---------------------------------------------------------------------------------

#---- ~SKCM TP ----
#         coef 	exp(coef) 	se(coef) 	z 	Pr(>|z|) 
# S11 	6.545e+00 	6.957e+02 	3.240e+00 	2.020 	0.043377 * 
# S23 	1.086e+01 	5.219e+04 	4.802e+00 	2.262 	0.023690 * 
# stageIII 	2.519e+00 	1.241e+01 	6.546e-01 	3.848 	0.000119 *** 

a <- as.data.frame(tp.cox[[7]])
b <- as.data.frame(tp.cox[[8]])
tp.cox.mod <- cbind(rownames(a), a, b[,c(3:4)])
colnames(tp.cox.mod) <- c("var","coef", "HR", "SE", "z", "p", "CI_low", "CI_high")

tp <- ggplot(tp.cox.mod, aes(x=var, y=coef)) + 
    scale_color_manual(values = c("(-Inf,0]" = "blue",
                                  "(0, Inf]" = "red")) +
    geom_errorbar(aes(ymin=coef-SE, ymax=coef+SE), width=0.05) +
    geom_point(aes(colour = cut(coef, c(-Inf, 0, Inf))),
               size = 3) +
    ggtitle("SKCM TP (n=27)") +
    ylab("log(HR)") + 
    scale_y_continuous(breaks=seq(-10,20,5)) +
    coord_cartesian(ylim=c(-8,18), expand=TRUE) + 
    theme_bw() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title = element_text(size=13),
          axis.text = element_text(size=12)) +
    scale_x_discrete(labels=c("S11", "S23", "Stage\nIII")) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    annotate("text", x=1, y=-8, label=paste0("p=",round(tp.cox.mod$p[1], digits=3))) +
    annotate("text", x=2, y=-8, label=paste0("p=",round(tp.cox.mod$p[2], digits=3))) +
    annotate("text", x=3, y=-8, label=paste0("p=",round(tp.cox.mod$p[3], digits=3)))
    

#---- ~SKCM TM ----
c <- as.data.frame(tm.cox[[7]])
d <- as.data.frame(tm.cox[[8]])
tm.cox.mod <- cbind(rownames(c), c, d[,c(3:4)])
tm.cox.mod$prognostic <- ifelse(tm.cox.mod$coef<0, "neg", "pos")
colnames(tm.cox.mod) <- c("var","coef", "HR", "SE", "z", "p", "CI_low", "CI_high")

tm <- ggplot(tm.cox.mod, aes(x=var, y=coef)) + 
    scale_color_manual(values = c("(-Inf,0]" = "blue",
                                  "(0, Inf]" = "red")) +
    geom_errorbar(aes(ymin=coef-SE, ymax=coef+SE), width=0.05) +
    geom_point(aes(colour = cut(coef, c(-Inf, 0, Inf))),
               size = 3) +
    ggtitle("SKCM TM (n=21)") +
    ylab("log(HR)") + 
    scale_y_continuous(breaks=seq(-10,20,5)) +
    coord_cartesian(ylim=c(-8,18), expand=TRUE) + 
    theme_bw() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title = element_text(size=13),
          axis.text = element_text(size=12)) +
    scale_x_discrete(labels=c("S6", "S7", "Interferon", "Interferon\n+ radiation", "Radiation")) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    annotate("text", x=1, y=-8, label=paste0("p=",round(tm.cox.mod$p[1], digits=3))) +
    annotate("text", x=2, y=-8, label=paste0("p=",round(tm.cox.mod$p[2], digits=3))) +
    annotate("text", x=3, y=-8, label=paste0("p=",round(tm.cox.mod$p[3], digits=3))) +
    annotate("text", x=4, y=-8, label=paste0("p=",round(tm.cox.mod$p[4], digits=3))) +
    annotate("text", x=5, y=-8, label=paste0("p=",round(tm.cox.mod$p[5], digits=3)))

# Combined plot, get skca and mel from respective cohort mutational analysis R files
pdf("./Figures/TCGA_ICGC_survival.pdf", w=8, h=10)
ggarrange(ggarrange(tp, tm,
                    widths = c(.8, 1.5), ncol=2),
          skca, melaau,
          nrow=3)
dev.off()


#--------- LOW HIGH SPLIT  ---------------------------------------------------------------------------------

# Split cohort into low and high potential groups
tp.surv.low <- tp.surv[which(tp.surv$met_potential=="low"),]        # dim 45 106  -- 24 alive, 21 dead
tp.surv.high <- tp.surv[which(tp.surv$met_potential=="high"),]      # dim 53 106  -- 47 alive, 6 dead

summary(tp.low.cox <- coxph(Surv(surv_cens, vital_status)~S1+S11+S23,
                data=tp.surv.low))
summary(tp.high.cox <- coxph(Surv(surv_cens, vital_status)~S11+S23,
                            data=tp.surv.high))


a<-summary(tp.low.cox <- coxph(Surv(surv_cens, vital_status)~stage+S7+S11+S23,
                            data=tp.surv.low))




#--------- TCGA + ICGC  ---------------------------------------------------------------------------------

# ---- ~ Combine data ----

# Add earlyLate disease progression to MELA-AU and SKCA-BR
mela.clinall$earlyLate <- ifelse(mela.clinall$stage=="I"|mela.clinall$stage=="II", "early", "late")
skca.clinall$earlyLate <- ifelse(skca.clinall$stage=="I"|skca.clinall$stage=="II", "early", "late")

# Extract relevant data (ID, age, survival, gender, vital status, mutsigs, stage, earlyLate, treatment)
tcga.tp <- tp.surv[,c(2, 3, 17, 27, 65, 71:101, 102, 103, 106)]
tcga.tp$surv_cens <- NA         # Make an empty column due to no treatment data
colnames(tcga.tp)[39] <- "treatment"

tcga.tm <- tm.surv[,c(2, 3, 17, 27, 65, 71:101, 102, 103, 106)]
icgc.br <- skca.clinall[,c(1, 7, 15, 3, 4, 32:62, 63, 64, 25)]
icgc.au <- mela.clinall[,c(1, 7, 15, 3, 4, 20:50, 65, 66, 58)]

# Change colnames to match
sigs <- colnames(icgc.br)[6:36]
colnames(tcga.tp) <- c("patient_id", "age", "survival_time", "gender", "vital_status", sigs, "stage", "progression", "treatment")
colnames(tcga.tm) <- c("patient_id", "age", "survival_time", "gender", "vital_status", sigs, "stage", "progression", "treatment")
colnames(icgc.br) <- c("patient_id", "age", "survival_time", "gender", "vital_status", sigs, "stage", "progression", "treatment")
colnames(icgc.au) <- c("patient_id", "age", "survival_time", "gender", "vital_status", sigs, "stage", "progression", "treatment")

# Add column with cohort
tcga.tp$cohort <- "SKCM TP"
tcga.tm$cohort <- "SKCM TM"
icgc.br$cohort <- "SKCA-BR"
icgc.au$cohort <- "MELA-AU"

# Add column with country
tcga.tp$country <- "USA"
tcga.tm$country <- "USA"
icgc.br$country <- "BRA"
icgc.au$country <- "AUS"

# Combine all cohorts together onto one variable
cut.mela <- rbind(tcga.tp, tcga.tm, icgc.br, icgc.au)

# Recode treatment for clarity
for (i in 1:nrow(cut.mela)) {
    if (!is.na(cut.mela$treatment[i])) {
        if(cut.mela$treatment[i]=="other therapy") {
            cut.mela$treatment[i] <- "immunotherapy + radiation"
        }
        else if (cut.mela$treatment[i]=="radiation therapy") {
            cut.mela$treatment[i] <- "radiation"
        }
    }
}

# Add combined NER and MMR signatures
cut.mela$NER <- cut.mela$S4 + cut.mela$S7 + cut.mela$S11 + cut.mela$S22 + cut.mela$S24 + cut.mela$S29
cut.mela$MMR <- cut.mela$S6 + cut.mela$S15 + cut.mela$S20 + cut.mela$S26


save(cut.mela, file="./Output/ICGC_TCGA_survival.RData")



# ---- ~ Univar survival ----
summary(mela.clin.cox <- coxph(Surv(survival_time, vital_status)~stage+age+cohort,
                            data=cut.mela))

summary(mela.geo.cox <- coxph(Surv(survival_time, vital_status)~cohort,
                               data=cut.mela))

summary(mela.sigs.cox <- coxph(Surv(survival_time, vital_status)~S1+S6+S11+S23+cohort,
                               data=cut.mela))

summary(mela.ther.cox <- coxph(Surv(survival_time, vital_status)~treatment+cohort,
                               data=cut.mela))


# ---- ~ KM plots univar ----

# Survival by disease progression
cut.prog.fitKM <- survfit(Surv(survival_time, vital_status)~progression,
                          data=cut.mela)

cut.prog.cox <- coxph(Surv(survival_time, vital_status)~progression,
                      data=cut.mela)
prog.cox <- summary(cut.prog.cox)

prog.hr <- round(prog.cox$coefficients[2], digits=2)
prog.CI1 <- round(prog.cox$conf.int[3], digits=2)
prog.CI2 <- round(prog.cox$conf.int[4], digits=2)
prog.pval <- round(prog.cox$sctest[3], digits=4)

pdf("./Figures/TCGA_ICGC_survival_progression_KM_annot.pdf", w=7, h=6)
plot(cut.prog.fitKM, col=c("blue","red"),
     mark.time = TRUE,
     xlab="Survival time (days)", ylab="% Overall survival",
     lwd = 3,
     cex.axis = 1.3,
     cex.lab = 1.5)
legend("topright", 
       legend = c("Early (n=194)", "Late (n=184)"), 
       fill=c("blue","red"), 
       bty="n",
       cex = 1.5)
legend("bottomleft", 
       legend = c("HR: ", paste0(prog.hr, ", 95% CI: ", prog.CI1, "-", prog.CI2), paste0("p = ", prog.pval)), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1,
       cex = 1.5)
dev.off()

# Survival by cohort
cut.coh.fitKM <- survfit(Surv(survival_time, vital_status)~cohort,
                         data=cut.mela)

cut.coh.cox <- coxph(Surv(survival_time, vital_status)~cohort,
                     data=cut.mela)
coh.cox <- summary(cut.coh.cox)

coh.hr <- round(coh.cox$coefficients[2], digits=2)
coh.CI1 <- round(coh.cox$conf.int[3], digits=2)
coh.CI2 <- round(coh.cox$conf.int[4], digits=2)
coh.pval <- round(coh.cox$sctest[3], digits=5)

pdf("./Figures/TCGA_ICGC_survival_cohort_KM_annot.pdf", w=7, h=6)
plot(cut.coh.fitKM, col=c("orange", "red","dark blue", "steelblue3"),
     mark.time = TRUE,
     xlab="Survival time (days)", ylab="% Overall survival",
     lwd = 3,
     cex.axis = 1.3,
     cex.lab = 1.5)
legend("topright", 
       legend = c("SKCM TP (n=27)", "SKCM TM (n=162)", "SKCA-BR (n=100)", "MELA-AU (n=104)"), 
       fill=c("steelblue3", "dark blue", "red", "orange"), 
       bty="n",
       cex = 1.5)
legend("bottomleft", 
       legend = "p < 2e-16", 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1,
       cex = 1.5)
dev.off()

# Survival by contort
cut.cont.fitKM <- survfit(Surv(survival_time, vital_status)~country,
                         data=cut.mela)

cut.cont.cox <- coxph(Surv(survival_time, vital_status)~country,
                     data=cut.mela)
cont.cox <- summary(cut.cont.cox)

cont.hr <- round(cont.cox$coefficients[2], digits=2)
cont.CI1 <- round(cont.cox$conf.int[3], digits=2)
cont.CI2 <- round(cont.cox$conf.int[4], digits=2)
cont.pval <- round(cont.cox$sctest[3], digits=12)

pdf("./Figures/TCGA_ICGC_survival_country_KM_annot.pdf", w=7, h=6)
plot(cut.cont.fitKM, col=c("orange", "red","dark blue"),
     mark.time = TRUE,
     xlab="Survival time (days)", ylab="% Overall survival",
     lwd = 3,
     cex.axis = 1.3,
     cex.lab = 1.5)
legend("topright", 
       legend = c("USA (n=189)", "Brazil (n=100)", "Australia (n=104)"), 
       fill=c("dark blue", "red", "orange"), 
       bty="n",
       cex = 1.5)
legend("bottomleft", 
       legend = paste0("p = ", cont.pval), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1,
       cex = 1.5)
dev.off()


# Survival by treatment
cut.mela1 <- cut.mela
cut.mela1$treatment1 <- ifelse(grepl("+ radiation", cut.mela1$treatment), "radiation + other", cut.mela1$treatment)
cut.mela1$treatment1 <- ifelse(!is.na(grepl(" + ", cut.mela1$treatment)), "radiation + other", cut.mela$treatment)

cut.treat.fitKM <- survfit(Surv(survival_time, vital_status)~treatment1,
                           data=cut.mela1)

cut.treat.cox <- coxph(Surv(survival_time, vital_status)~treatment1,
                       data=cut.mela1)
treat.cox <- summary(cut.treat.cox)

treat.hr <- round(treat.cox$coefficients[2], digits=2)
treat.CI1 <- round(treat.cox$conf.int[3], digits=2)
treat.CI2 <- round(treat.cox$conf.int[4], digits=2)
treat.pval <- round(treat.cox$sctest[3], digits=4)

treat.col <- brewer.pal(7, "Paired")

pdf("./Figures/TCGA_ICGC_survival_treatment_KM_annot.pdf", w=7, h=6)
plot(cut.treat.fitKM, col=treat.col,
     mark.time = TRUE,
     xlab="Survival time (days)", ylab="% Overall survival",
     lwd = 3,
     cex.axis = 1.3,
     cex.lab = 1.5)
legend("bottomleft", 
       legend = paste0("p = ", treat.pval), 
       bty = "n",
       y.intersp=0.8,x.intersp=0.5,text.width=0.1,
       cex = 1.5)
dev.off()

# Legend alone 
leg.fitKM <- survfit(Surv(survival_time, vital_status)~cohort,
                     data=cut.mela[which(cut.mela$cohort=="SKCM TP"),])

pdf("./Figures/TCGA_ICGC_survival_treatment_KM_leg.pdf", w=7, h=6)
plot(leg.fitKM, col=treat.col,
     mark.time = TRUE,
     xlab="Survival time (days)", ylab="% Overall survival",
     lwd = 3,
     xlim=c(0,10000),
     ylim=c(0,2),
     cex.axis = 1.3,
     cex.lab = 1.5)
legend("topright", 
       legend = c("chemotherapy (n=18)", "immunotherapy (n=3)", 
                  "interferon (n=10)", "no treatment (n=68)",
                  "radiation (n=19)", "radiation + other (n=4)",
                  "surgery (n=57)"), 
       fill=treat.col, 
       bty="n",
       cex = 1.5)
dev.off()


# ---- ~ Multivar survival ----

# Multivariate no therapy
cut <- summary(mela.cox <- coxph(Surv(survival_time, vital_status)~stage+age+
                              S1+MMR,
                          data=cut.mela))

cut.coh <- summary(mela.cox <- coxph(Surv(survival_time, vital_status)~stage+age+
                              NER+cohort,
                               data=cut.mela))

# Multivariate with therapy
cut.ther <- summary(mela.cox.ther <- coxph(Surv(survival_time, vital_status)~stage+treatment+
                              S6+cohort,
                          data=cut.mela))
cut.ther <- summary(mela.cox.ther <- coxph(Surv(survival_time, vital_status)~stage+treatment+
                                               S6,
                                           data=cut.mela))
# Multivariate with therapy with radiation + immuno and radiation + interferon combined
cut.ther1 <- summary(mela.cox.ther1 <- coxph(Surv(survival_time, vital_status)~stage+treatment1+
                                   S6+cohort,
                               data=cut.mela1))
cut.ther1 <- summary(mela.cox.ther1 <- coxph(Surv(survival_time, vital_status)~stage+treatment1+
                                                 S6,
                                             data=cut.mela1))





summary(tp.cox <- coxph(Surv(days_to_death, vital_status)~S7+S11+S22+S23+S9,
                               data=tp.surv))

summary(tm.cox <- coxph(Surv(days_to_death, vital_status)~S7+S11+S22+S23+S9,
                         data=tm.surv))

summary(mela.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~S7+S11+S22+S23+S9,
                         data=mela.clinall))

summary(skca.cox <- coxph(Surv(donor_survival_time, donor_vital_status)~S7+S11+S22+S23+S9,
                          data=skca.clinall))

summary(cut.cox <- coxph(Surv(survival_time, vital_status)~cohort+S7+S11+S22+S23+S9,
                          data=cut.mela))


summary(cut.cox <- coxph(Surv(survival_time, vital_status)~S6+S15,
                         data=cut.mela))




#--------- FINAL MODEL PLOTS  ---------------------------------------------------------------------------------

# No therapy, no cohort
i <- as.data.frame(cut[[7]])
j <- as.data.frame(cut[[8]])
cut.mod <- cbind(rownames(i), i, j[,c(3:4)])
colnames(cut.mod) <- c("var","coef", "HR", "SE", "z", "p", "CI_low", "CI_high")

ggplot(cut.mod[3:6,], aes(x=factor(var, cut.mod$var), y=coef)) + 
    scale_color_manual(values = c("(-Inf,0]" = "blue",
                                  "(0, Inf]" = "red")) +
    geom_errorbar(aes(ymin=coef-SE, ymax=coef+SE), width=0.05) +
    geom_point(aes(colour = cut(coef, c(-Inf, 0, Inf))),
               size = 3) +
    ggtitle("") +
    ylab("log(HR)") + 
    scale_y_continuous(breaks=seq(-10,20,1)) +
    coord_cartesian(ylim=c(-1,4), expand=TRUE) + 
    theme_bw() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title = element_text(size=13),
          axis.text = element_text(size=12)) +
    #scale_x_discrete(labels=c("Age", "MMR\nsigs", "S1","Stage\nIII")) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    annotate("text", x=1, y=-1, label=paste0("p=",round(cut.mod$p[3], digits=3))) +
    annotate("text", x=2, y=-1, label=paste0("p=",round(cut.mod$p[4], digits=4))) +
    annotate("text", x=3, y=-1, label=paste0("p=",round(cut.mod$p[5], digits=5))) +
    annotate("text", x=4, y=-1, label=paste0("p=",round(cut.mod$p[6], digits=3)))


# No therapy, cohort
k <- as.data.frame(cut.coh[[7]])
l <- as.data.frame(cut.coh[[8]])
cut.coh.mod <- cbind(rownames(k), k, l[,c(3:4)])
colnames(cut.coh.mod) <- c("var","coef", "HR", "SE", "z", "p", "CI_low", "CI_high")

var.ord <- cut.coh.mod$var

ggplot(cut.coh.mod[3:8,], aes(x=factor(var, levels=var.ord), y=coef))
wother <- ggplot(cut.coh.mod, aes(x=factor(var, levels=var.ord), y=coef)) + 
    scale_color_manual(values = c("(-Inf,0]" = "blue",
                                  "(0, Inf]" = "red")) +
    geom_errorbar(aes(ymin=coef-SE, ymax=coef+SE), width=0.05) +
    geom_point(aes(colour = cut(coef, c(-Inf, 0, Inf))),
               size = 3) +
    ggtitle("") +
    ylab("log(HR)") + 
    scale_y_continuous(breaks=seq(-10,20,1)) +
    coord_cartesian(ylim=c(-1,3), expand=TRUE) + 
    theme_bw() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title = element_text(size=15),
          axis.text = element_text(size=14)) +
    scale_x_discrete(labels=c("Stage\nII","Stage\nIII","Stage\nIV", "Age", "NER\nsigs","SKCA\nBR", "SKCM\nTM", "SKCM\nTP")) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    annotate("text", x=1, y=-1, label=paste0("p=",round(cut.coh.mod$p[1], digits=3))) +
    annotate("text", x=2, y=-1, label=paste0("p=",round(cut.coh.mod$p[2], digits=3))) +
    annotate("text", x=3, y=-1, label=paste0("p=",round(cut.coh.mod$p[3], digits=3))) +
    annotate("text", x=4, y=-1, label=paste0("p=",round(cut.coh.mod$p[4], digits=7))) +
    annotate("text", x=5, y=-1, label=paste0("p=",round(cut.coh.mod$p[5], digits=3))) +
    annotate("text", x=6, y=-1, label=paste0("p=",round(cut.coh.mod$p[6], digits=3))) +
    annotate("text", x=7, y=-1, label=paste0("p=",round(cut.coh.mod$p[7], digits=8))) +
    annotate("text", x=8, y=-1, label=paste0("p<2e-16"))


# Therapy
m <- as.data.frame(cut.ther1[[7]])
n <- as.data.frame(cut.ther1[[8]])
cut.ther.mod <- cbind(rownames(m), m, n[,c(3:4)])
colnames(cut.ther.mod) <- c("var","coef", "HR", "SE", "z", "p", "CI_low", "CI_high")

var.ord2 <- cut.ther.mod$var

wther <- ggplot(cut.ther.mod, aes(x=factor(var, levels=var.ord2), y=coef)) + 
    scale_color_manual(values = c("(-Inf,0]" = "blue",
                                  "(0, Inf]" = "red")) +
    geom_errorbar(aes(ymin=coef-SE, ymax=coef+SE), width=0.05) +
    geom_point(aes(colour = cut(coef, c(-Inf, 0, Inf))),
               size = 3) +
    ggtitle("") +
    ylab("log(HR)") + 
    scale_y_continuous(breaks=seq(-10,100,2)) +
    coord_cartesian(ylim=c(-3.5,12), expand=TRUE) + 
    theme_bw() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title = element_text(size=13),
          axis.text = element_text(size=12)) +
    scale_x_discrete(labels=c("Stage\nII", "Stage\nIII", "Stage\nIV", "Immuno-\ntherapy", "Interferon",
                              "No\ntreatment", "Radiation", "Radiation\n+ other", "Surgery","S6\n(MMR)",
                              "SKCA\nBR", "SKCM\nTM")) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    annotate("text", x=1, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[1], digits=3))) +
    annotate("text", x=2, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[2], digits=3))) +
    annotate("text", x=3, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[3], digits=3))) +
    annotate("text", x=4, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[4], digits=3))) +
    annotate("text", x=5, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[5], digits=3))) +
    annotate("text", x=6, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[6], digits=3))) +
    annotate("text", x=7, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[7], digits=3))) +
    annotate("text", x=8, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[8], digits=3))) +
    annotate("text", x=9, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[9], digits=4))) +
    annotate("text", x=10, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[10], digits=3)))
    #annotate("text", x=11, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[11], digits=3))) +
    #annotate("text", x=12, y=-3.5, label=paste0("p=",round(cut.ther.mod$p[12], digits=3)))

pdf("./Figures/TCGA_ICGC_comb_survival2.pdf", w=10, h=10)
ggarrange(wother, wther,
          nrow=2,
          labels = c("a)", "b)"))
dev.off()