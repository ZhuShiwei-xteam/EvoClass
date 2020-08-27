setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")
library(dplyr)
library(survival)
library(survminer)

load("data/GliomaMolecularSubtype.RData")
load("data/gliomaClinicalData.RData")
rownames(pan.glio.mer.cli.data) <- pan.glio.mer.cli.data$bcr_patient_barcode
load("data/PanGlioma.clinical.outcome.endpoints.RData")
### OS
clonal.subtype <- glioma.molecular.subtype[, 3, drop=F]
samples <- PanGlioma.OS.data$bcr_patient_barcode
Glioma.OS.data <- cbind(cbind(pan.glio.mer.cli.data[samples, , drop =FALSE], PanGlioma.OS.data[samples, c("OS", "OS.time") , drop =FALSE]), clonal.subtype[samples, , drop =FALSE])

# for each subtypes
subtypes <- unique(clonal.subtype[,1])
lapply(subtypes, function(x){
	per.subtype.clinical.data <- subset(Glioma.OS.data, final.subtype == x)
	
	pdf(paste("result/6.Survival_analysis_for_the_Modified_Glioma_Classification_System/Statistic_distr_of_clinical_factors/",x, ".pdf", sep=""))
	# km plot
	per.subtype.clinical.data$OS.time <- per.subtype.clinical.data$OS.time/365
	fit <- survfit(Surv(OS.time, OS)~final.subtype, data = per.subtype.clinical.data)
	ggsurvplot(fit,risk.table=FALSE, conf.int=FALSE, censor = FALSE)
	
	# Age distribution curve
	y <- density(per.subtype.clinical.data$age_at_initial_pathologic_diagnosis) 
	plot(y, main = "", xlab ="", ylab ="")
	
	# Histological type histogram
	histological.statistic <- as.data.frame(table(per.subtype.clinical.data$molecular_histological_type))
	barplot(histological.statistic$Freq, names.arg=histological.statistic$Var1, xlab="", ylab="", main="")
	
	# who grade histogram
	grade.statistic <- as.data.frame(table(per.subtype.clinical.data$neoplasm_histologic_grade))
	barplot(grade.statistic$Freq, names.arg=grade.statistic$Var1, xlab="", ylab="", main="")
	
	dev.off()
})