setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")
library(survival)
source("code/6.Survival_analysis_for_the_Modified_Glioma_Classification_System/Survival_Analysis_Functions.r")

# load clinical outcome data
load("/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Resources/PanGlioma.clinical.outcome.endpoints.RData")

load("/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Resources/GliomaMolecularSubtype.RData")
clonal.subtype <- glioma.molecular.subtype[, 3, drop=F]
clonal.subtype$final.subtype <- factor(clonal.subtype$final.subtype, levels = c("IDHmut.codel","IDHcolmut.non.codel", "IDHsubmut.non.codel", "IDHwt")) 	
	
file.path <- "result/6.Survival_analysis_for_the_Modified_Glioma_Classification_System/"
########## OS

# clonal.subtype <- subset(clonal.subtype, final.subtype %in% c("IDHcolmut.non.codel", "IDHsubmut.non.codel"))
subtype.index <- which(rownames(PanGlioma.OS.data) %in% rownames(clonal.subtype))
clinical.outcome <- ClinicalOutcome(PanGlioma.OS.data[subtype.index,,drop=F], "OS")
clinical.outcome$time <- clinical.outcome$time/365
# 生存曲线的绘制
colors <- c("#CF802D", "#BF5354", "#3878A6","#266479")
clonalSubtype.os.KM.plot <- KMplot(clinical.outcome, clonal.subtype,colors)
########## DSS
subtype.index <- which(rownames(PanGlioma.DSS.data) %in% rownames(clonal.subtype))
clinical.outcome <- ClinicalOutcome(PanGlioma.DSS.data[subtype.index,,drop=F], "DSS")
clinical.outcome$time <- clinical.outcome$time/365
# 生存曲线的绘制
clonalSubtype.dss.KM.plot <- KMplot(clinical.outcome, clonal.subtype, colors)

pdf(paste0(c(file.path,"clonal_subtype_KM", ".pdf"), collapse = ""))
clonalSubtype.os.KM.plot
clonalSubtype.dss.KM.plot
dev.off()

###################################################################################################
################### Univariate and Multivariate Cox regression

######### Subtype information, clinical characteristics information and outcome were combined into one data set
load("data/gliomaClinicalData.RData")
rownames(pan.glio.mer.cli.data) <- pan.glio.mer.cli.data$bcr_patient_barcode
pan.glio.mer.cli.data$Extent_of_surgery_resection[which(pan.glio.mer.cli.data$Extent_of_surgery_resection == "Biopsy")] <- NA
load("data/PanGlioma.clinical.outcome.endpoints.RData")

############################# OS
samples <- PanGlioma.OS.data$bcr_patient_barcode
Glioma.OS.data <- cbind(cbind(pan.glio.mer.cli.data[samples, , drop =FALSE], PanGlioma.OS.data[samples, c("OS", "OS.time") , drop =FALSE]), clonal.subtype[samples, , drop =FALSE])
survival.type <- "OS"

clinical.factors <- c("age_at_initial_pathologic_diagnosis","neoplasm_histologic_grade", "final.subtype")
clinical.data <- RmnaBeforeMulCox(Glioma.OS.data, clinical.factors)# 873 samples
clinical.data$final.subtype <- factor(clinical.data$final.subtype, levels = c("IDHmut.codel", "IDHcolmut.non.codel", "IDHsubmut.non.codel", "IDHwt"))

# Univariate Cox regression
cox.univ <- lapply(clinical.factors, function(x) CoxUniv(clinical.data, interData = T, clinical.factor = x, survival.type = survival.type))
cox.univ.result <- lapply(cox.univ, function(x) CoxphResult(x))
names(cox.univ.result) <- clinical.factors
# Multivariate Cox regression
clinical.factors <- c("age_at_initial_pathologic_diagnosis","neoplasm_histologic_grade", "final.subtype")
clonal.a.g.coxmulti <- CoxMulti(clinical.data, clinical.factors = clinical.factors, survival.type = survival.type)
cox.multi.result <- CoxphResult(clonal.a.g.coxmulti)
a.g.coxmulti <- CoxMulti(clinical.data, clinical.factors = c("age_at_initial_pathologic_diagnosis","neoplasm_histologic_grade"), survival.type = survival.type)
clonal.LRT <- anova(a.g.coxmulti,clonal.a.g.coxmulti,test = 'Chisq')$`P(>|Chi|)`[2] # 5.739975e-12
save(cox.univ.result, cox.multi.result, file = paste(file.path,"COX.result.OS.RData", sep=""))

### C-index 
clinical.outcome <- ClinicalOutcome(clinical.data, "OS")
# c-index of each clinical factor
models <- list(clonal="final.subtype",
               #original="IDH.codel.subtype", 
               age="age_at_initial_pathologic_diagnosis",
               grade="neoplasm_histologic_grade")
sig.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, models,diff=F)
# c-index of multivariate cox model
model <- list(clonal=c("age_at_initial_pathologic_diagnosis", "neoplasm_histologic_grade", "final.subtype"))
multi.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, model,diff=F)
save(sig.cindex, multi.cindex, file = paste(file.path, "cindex.OS.RData", sep=""))

############################## DSS
samples <- PanGlioma.DSS.data$bcr_patient_barcode
Glioma.DSS.data <- cbind(cbind(pan.glio.mer.cli.data[samples, , drop =FALSE], PanGlioma.DSS.data[samples, c("DSS", "DSS.time") , drop =FALSE]), clonal.subtype[samples, , drop =FALSE])
survival.type <- "DSS"

######### Univariate and Multivariate Cox regression
clinical.factors <- c("age_at_initial_pathologic_diagnosis","neoplasm_histologic_grade", "final.subtype")
clinical.data <- RmnaBeforeMulCox(Glioma.DSS.data, clinical.factors)# 832个样本
clinical.data$final.subtype <- factor(clinical.data$final.subtype, levels = c("IDHmut.codel", "IDHcolmut.non.codel", "IDHsubmut.non.codel", "IDHwt"))
# Univariate Cox regression
cox.univ <- lapply(clinical.factors, function(x) CoxUniv(clinical.data, interData = T, clinical.factor = x, survival.type = survival.type))
cox.univ.result <- lapply(cox.univ, function(x) CoxphResult(x))
names(cox.univ.result) <- clinical.factors
# Multivariate Cox regression
clinical.factors <- c("age_at_initial_pathologic_diagnosis","neoplasm_histologic_grade", "final.subtype")
clonal.a.g.coxmulti <- CoxMulti(clinical.data, clinical.factors = clinical.factors, survival.type = survival.type)
cox.multi.result <- CoxphResult(clonal.a.g.coxmulti)
a.g.coxmulti <- CoxMulti(clinical.data, clinical.factors = c("age_at_initial_pathologic_diagnosis","neoplasm_histologic_grade"), survival.type = survival.type)
clonal.LRT <- anova(a.g.coxmulti,clonal.a.g.coxmulti,test = 'Chisq')$`P(>|Chi|)`[2] # 7.049076e-13
save(cox.univ.result, cox.multi.result, file = paste(file.path, "COX.result.DSS.RData", sep=""))

### C-index 
clinical.outcome <- ClinicalOutcome(clinical.data, "DSS")
#  c-index of each clinical factor
models <- list(clonal="final.subtype",
               #original="IDH.codel.subtype", 
               age="age_at_initial_pathologic_diagnosis",
               grade="neoplasm_histologic_grade")
sig.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, models,diff=F)
# c-index of multivariate cox model
model <- list(clonal=c("age_at_initial_pathologic_diagnosis", "neoplasm_histologic_grade", "final.subtype"))
multi.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, model,diff=F)
save(sig.cindex, multi.cindex, file = paste(file.path,"cindex.DSS.RData", sep=""))

