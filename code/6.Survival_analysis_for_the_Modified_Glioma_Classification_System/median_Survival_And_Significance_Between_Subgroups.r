setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")

library(survival)
source("code/6.Survival_analysis_for_the_Modified_Glioma_Classification_System/Survival_Analysis_Functions.r")

load("data/GliomaMolecularSubtype.RData")
clonal.subtype <- glioma.molecular.subtype[, 3, drop=F]
load("data/PanGlioma.clinical.outcome.endpoints.RData")

subtype <- unique(clonal.subtype)[,1]
############# OS
subtype.index <- which(rownames(PanGlioma.OS.data) %in% rownames(clonal.subtype))
clinical.data <- PanGlioma.OS.data[subtype.index,,drop=F]
survival.type <- "OS"
survival.time.info <- NULL
p.value <- NULL
d=1
for(i in 1:3){
  for(j in (i+1):4){
    ## extracte samples containing two subtypes
    cancer.subtype <- clonal.subtype[which(clonal.subtype[,1] %in% subtype[c(i,j)]), ,drop=F]
    commonSamples <- intersect(rownames(cancer.subtype), rownames(clinical.data))
    clinicalData <- cbind(clinical.data[commonSamples, ,drop = FALSE],cancer.subtype[commonSamples, ,drop = FALSE])
    clinical.outcome <- ClinicalOutcome(clinicalData, survival.type)
    # The significance between subgroups
    km <- survdiff(Surv(clinical.outcome$time/365, clinical.outcome$event) ~ final.subtype, data = clinicalData)
    pair.p.value <- 1 - pchisq(km$chisq, length(km$n) - 1)
    p.value <- c(p.value,list(pair.p.value))
    names(p.value)[d] <- paste(subtype[c(i,j)],collapse = ".vs.")
    # Median lifetime of each subtype
    fit <- survfit(Surv(clinical.outcome$time/365, clinical.outcome$event) ~ final.subtype, data = clinicalData)
    survival.time.info <- rbind(survival.time.info, summary(fit)$table)
    d <- d+1
  }
}
median.survival <- survival.time.info[!duplicated(survival.time.info[,7:9]),7:9]
save(median.survival,p.value,file="result/6.Survival_analysis_for_the_Modified_Glioma_Classification_System/median.survival.P.OS.RData")

#########################DSS
subtype.index <- which(rownames(PanGlioma.DSS.data) %in% rownames(clonal.subtype))
clinical.data <- PanGlioma.DSS.data[subtype.index,,drop=F]
survival.type <- "DSS"
survival.time.info <- NULL
p.value <- NULL
d=1
for(i in 1:3){
  for(j in (i+1):4){
    ## extracte samples containing two subtypes
    cancer.subtype <- clonal.subtype[which(clonal.subtype[,1] %in% subtype[c(i,j)]), ,drop=F]
    commonSamples <- intersect(rownames(cancer.subtype), rownames(clinical.data))
    clinicalData <- cbind(clinical.data[commonSamples, ,drop = FALSE],cancer.subtype[commonSamples, ,drop = FALSE])
    clinical.outcome <- ClinicalOutcome(clinicalData, survival.type)
    # The significance between subgroups
    km <- survdiff(Surv(clinical.outcome$time/365, clinical.outcome$event) ~ final.subtype, data = clinicalData)
    pair.p.value <- 1 - pchisq(km$chisq, length(km$n) - 1)
    p.value <- c(p.value,list(pair.p.value))
    names(p.value)[d] <- paste(subtype[c(i,j)],collapse = ".vs.")
    # Median lifetime of each subtype
    fit <- survfit(Surv(clinical.outcome$time/365, clinical.outcome$event) ~ final.subtype, data = clinicalData)
    survival.time.info <- rbind(survival.time.info, summary(fit)$table)
    d <- d+1
  }
}
median.survival <- survival.time.info[!duplicated(survival.time.info[,7:9]),7:9]
save(median.survival,p.value,file="result/6.Survival_analysis_for_the_Modified_Glioma_Classification_System/median.survival.P.DSS.RData")

