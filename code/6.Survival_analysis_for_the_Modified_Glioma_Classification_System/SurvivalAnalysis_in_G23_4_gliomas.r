setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")

source("code/6.Survival_analysis_for_the_Modified_Glioma_Classification_System/Survival_Analysis_Functions.r")


load("data/PanGlioma.clinical.outcome.endpoints.RData")
load("data/GliomaMolecularSubtype.RData")
glioma.molecular.subtype$final.subtype <- factor(glioma.molecular.subtype$final.subtype, levels = c("IDHmut.codel","IDHcolmut.non.codel", "IDHsubmut.non.codel", "IDHwt")) 	
file.path <- "result/6.Survival_analysis_for_the_Modified_Glioma_Classification_System/"

G23.clonal.subtype <- subset(glioma.molecular.subtype, grade == "grade_II-III", select = final.subtype) # 503
G4.clonal.subtype <- subset(glioma.molecular.subtype, grade == "grade_IV", select = final.subtype) # 372
G4.clonal.subtype[,] <- factor(G4.clonal.subtype[,], levels = c("IDHcolmut.non.codel", "IDHsubmut.non.codel", "IDHwt"))
########## OS
# grade 2/3
subtype.index <- which(rownames(PanGlioma.OS.data) %in% rownames(G23.clonal.subtype))
clinical.outcome <- ClinicalOutcome(PanGlioma.OS.data[subtype.index,,drop=F], "OS")
clinical.outcome$time <- clinical.outcome$time/365
# km plot
colors <- c("#CF802D", "#BF5354", "#3878A6","#266479")
G23.clonalSubtype.os.KM.plot <- KMplot(clinical.outcome, G23.clonal.subtype,colors)
# grade 4
subtype.index <- which(rownames(PanGlioma.OS.data) %in% rownames(G4.clonal.subtype))
clinical.outcome <- ClinicalOutcome(PanGlioma.OS.data[subtype.index,,drop=F], "OS")
clinical.outcome$time <- clinical.outcome$time/365
# km plot
colors <- c("#BF5354", "#3878A6","#266479")
G4.clonalSubtype.os.KM.plot <- KMplot(clinical.outcome, G4.clonal.subtype,colors)

########## DSS
# grade 2/3
subtype.index <- which(rownames(PanGlioma.DSS.data) %in% rownames(G23.clonal.subtype))
clinical.outcome <- ClinicalOutcome(PanGlioma.DSS.data[subtype.index,,drop=F], "DSS")
clinical.outcome$time <- clinical.outcome$time/365
# km plot
colors <- c("#CF802D", "#BF5354", "#3878A6","#266479")
G23.clonalSubtype.dss.KM.plot <- KMplot(clinical.outcome, G23.clonal.subtype, colors)
# grade 4
subtype.index <- which(rownames(PanGlioma.DSS.data) %in% rownames(G4.clonal.subtype))
clinical.outcome <- ClinicalOutcome(PanGlioma.DSS.data[subtype.index,,drop=F], "DSS")
clinical.outcome$time <- clinical.outcome$time/365
# km plot
colors <- c("#BF5354", "#3878A6","#266479")
G4.clonalSubtype.dss.KM.plot <- KMplot(clinical.outcome, G4.clonal.subtype, colors)
	
pdf(paste0(c(file.path,"Strify_Grade_clonal_subtype_KM", ".pdf"), collapse = ""))
	G23.clonalSubtype.os.KM.plot
	G23.clonalSubtype.dss.KM.plot
	G4.clonalSubtype.os.KM.plot
	G4.clonalSubtype.dss.KM.plot 
dev.off()

### Univariate and multivariate Cox regression analysis
load("data/gliomaClinicalData.RData")
rownames(pan.glio.mer.cli.data) <- pan.glio.mer.cli.data$bcr_patient_barcode

## OS
# grage 2/3 G23.clonal.subtype
# data preparation 
clinical.factors <- c("age_at_initial_pathologic_diagnosis", "neoplasm_histologic_grade", "final.subtype")
samples <- intersect(PanGlioma.OS.data$bcr_patient_barcode, rownames(G23.clonal.subtype))
Glioma.OS.data <- cbind(cbind(pan.glio.mer.cli.data[samples, , drop =FALSE], PanGlioma.OS.data[samples, c("OS", "OS.time") , drop =FALSE]), G23.clonal.subtype[samples, , drop =FALSE])
survival.type <- "OS"
clinical.data <- RmnaBeforeMulCox(Glioma.OS.data, clinical.factors) # 501

# Univariate Cox regression analysis
G23.cox.univ <- lapply(clinical.factors, function(x) CoxUniv(clinical.data, interData = T, clinical.factor = x, survival.type = survival.type))
G23.cox.univ.result <- lapply(G23.cox.univ, function(x) CoxphResult(x))
names(G23.cox.univ.result) <- clinical.factors

# multivariate Cox regression analysis
G23.clonal.coxmulti <- CoxMulti(clinical.data, clinical.factors = clinical.factors, survival.type = survival.type)
G23.clonal.multi.result <- CoxphResult(G23.clonal.coxmulti)

### C-index 
clinical.outcome <- ClinicalOutcome(clinical.data, "OS")
models <- list(clonal="final.subtype",
               age="age_at_initial_pathologic_diagnosis",
               grade="neoplasm_histologic_grade")
G23.sig.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, models,diff=F)
model <- list(clonal= clinical.factors)
G23.multi.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, model,diff=F)
save(G23.cox.univ, G23.cox.univ.result, G23.clonal.coxmulti, G23.clonal.multi.result, G23.sig.cindex, G23.multi.cindex, file = paste(file.path, "grade2_3_OS_cox_result.RData", sep=""))

# grade 4 G4.clonal.subtype
# data preparation 
clinical.factors <- c("age_at_initial_pathologic_diagnosis", "final.subtype")
samples <- intersect(PanGlioma.OS.data$bcr_patient_barcode, rownames(G4.clonal.subtype))
Glioma.OS.data <- cbind(cbind(pan.glio.mer.cli.data[samples, , drop =FALSE], PanGlioma.OS.data[samples, c("OS", "OS.time") , drop =FALSE]), G4.clonal.subtype[samples, , drop =FALSE])
survival.type <- "OS"
clinical.data <- RmnaBeforeMulCox(Glioma.OS.data, clinical.factors)# 372
# Univariate Cox regression analysis
G4.cox.univ <- lapply(clinical.factors, function(x) CoxUniv(clinical.data, interData = T, clinical.factor = x, survival.type = survival.type))
G4.cox.univ.result <- lapply(G4.cox.univ, function(x) CoxphResult(x))
names(G4.cox.univ.result) <- clinical.factors

# multivariate Cox regression analysis
G4.clonal.coxmulti <- CoxMulti(clinical.data, clinical.factors = clinical.factors, survival.type = survival.type)
G4.clonal.multi.result <- CoxphResult(G4.clonal.coxmulti)

### C-index 
clinical.outcome <- ClinicalOutcome(clinical.data, "OS")
models <- list(clonal="final.subtype",
               age="age_at_initial_pathologic_diagnosis")
G4.sig.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, models,diff=F)
model <- list(clonal= clinical.factors)
G4.multi.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, model,diff=F)
save(G4.cox.univ, G4.cox.univ.result, G4.clonal.coxmulti, G4.clonal.multi.result, G4.sig.cindex, G4.multi.cindex, file = paste(file.path,"grade4_OS_cox_result.RData", sep=""))


# DSS
# grade 2/3 G23.clonal.subtype
# data preparation 
clinical.factors <- c("age_at_initial_pathologic_diagnosis", "neoplasm_histologic_grade", "final.subtype")
samples <- intersect(PanGlioma.DSS.data$bcr_patient_barcode, rownames(G23.clonal.subtype))
Glioma.DSS.data <- cbind(cbind(pan.glio.mer.cli.data[samples, , drop =FALSE], PanGlioma.DSS.data[samples, c("DSS", "DSS.time") , drop =FALSE]), G23.clonal.subtype[samples, , drop =FALSE])
survival.type <- "DSS"
clinical.data <- RmnaBeforeMulCox(Glioma.DSS.data, clinical.factors)# 493

# Univariate Cox regression analysis
G23.cox.univ <- lapply(clinical.factors, function(x) CoxUniv(clinical.data, interData = T, clinical.factor = x, survival.type = survival.type))
G23.cox.univ.result <- lapply(G23.cox.univ, function(x) CoxphResult(x))
names(G23.cox.univ.result) <- clinical.factors

# multivariate Cox regression analysis
G23.clonal.coxmulti <- CoxMulti(clinical.data, clinical.factors = clinical.factors, survival.type = survival.type)
G23.clonal.multi.result <- CoxphResult(G23.clonal.coxmulti)

### C-index 
clinical.outcome <- ClinicalOutcome(clinical.data, "DSS")
models <- list(clonal="final.subtype",
               age="age_at_initial_pathologic_diagnosis",
               grade="neoplasm_histologic_grade")
G23.sig.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, models,diff=F)
model <- list(clonal= clinical.factors)
G23.multi.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, model,diff=F)
save(G23.cox.univ, G23.cox.univ.result, G23.clonal.coxmulti, G23.clonal.multi.result, G23.sig.cindex, G23.multi.cindex, file = paste(file.path,"grade2_3_DSS_cox_result.RData", sep=""))

# grade 4 G4.clonal.subtype
# data preparation
clinical.factors <- c("age_at_initial_pathologic_diagnosis", "final.subtype")
samples <- intersect(PanGlioma.DSS.data$bcr_patient_barcode, rownames(G4.clonal.subtype))
Glioma.DSS.data <- cbind(cbind(pan.glio.mer.cli.data[samples, , drop =FALSE], PanGlioma.DSS.data[samples, c("DSS", "DSS.time") , drop =FALSE]), G4.clonal.subtype[samples, , drop =FALSE])
survival.type <- "DSS"
clinical.data <- RmnaBeforeMulCox(Glioma.DSS.data, clinical.factors)# 339
# Univariate Cox regression analysis
G4.cox.univ <- lapply(clinical.factors, function(x) CoxUniv(clinical.data, interData = T, clinical.factor = x, survival.type = survival.type))
G4.cox.univ.result <- lapply(G4.cox.univ, function(x) CoxphResult(x))
names(G4.cox.univ.result) <- clinical.factors

#  multivariate Cox regression analysis
G4.clonal.coxmulti <- CoxMulti(clinical.data, clinical.factors = clinical.factors, survival.type = survival.type)
G4.clonal.multi.result <- CoxphResult(G4.clonal.coxmulti)

### C-index 
clinical.outcome <- ClinicalOutcome(clinical.data, "DSS")
models <- list(clonal="final.subtype",
               age="age_at_initial_pathologic_diagnosis")
G4.sig.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, models,diff=F)
model <- list(clonal= clinical.factors)
G4.multi.cindex <- cIndex(clinical.data,time=clinical.outcome$time, event=clinical.outcome$event, model,diff=F)
save(G4.cox.univ, G4.cox.univ.result, G4.clonal.coxmulti, G4.clonal.multi.result, G4.sig.cindex, G4.multi.cindex, file = paste(file.path,"grade4_DSS_cox_result.RData", sep=""))

