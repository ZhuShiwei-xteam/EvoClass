setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')

library('powerSurvEpi')

source('code/3.Glioma_Maunally_defined_subtypes/Cox.function.R')
load('data/GliomaMolecularSubtype.RData')
load('data/PanGlioma.clinical.outcome.endpoints.RData')

PanGlioma.OS.data <- merge(glioma.molecular.subtype, PanGlioma.OS.data, by = 'row.names')
PanGlioma.OS.data <- subset(PanGlioma.OS.data, explore.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))

# univariate cox proportional hazards analysis
Cox.function(time = PanGlioma.OS.data$OS.time, event = PanGlioma.OS.data$OS, 
 clinical.data = PanGlioma.OS.data, clinical.variate = 3)

# IDHcolmut.non.codel IDHsubmut.non.codel 
                # 235                  34 

  # univ HR (95% CI for HR) univ p value
# 1   2.2083 (1.0661-4.574)        0.033
 
PanGlioma.DSS.data <- merge(glioma.molecular.subtype, PanGlioma.DSS.data, by = 'row.names')
PanGlioma.DSS.data <- subset(PanGlioma.DSS.data, explore.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))

Cox.function(time = PanGlioma.DSS.data$DSS.time, event = PanGlioma.DSS.data$DSS, 
 clinical.data = PanGlioma.DSS.data, clinical.variate = 3)

# IDHcolmut.non.codel IDHsubmut.non.codel 
                # 228                  33 

  # univ HR (95% CI for HR) univ p value
# 1  2.3592 (1.0878-5.1165)       0.0298


# Sample size calculation for the Comparison of Survival Curves 
# Between Two Groups under the Cox Proportional-Hazards Model for clinical trials

PanGlioma.OS.data$explore.subtype <- ifelse(PanGlioma.OS.data$explore.subtype == 'IDHsubmut.non.codel', 'E', 'C')
ssizeCT(formula = Surv(OS.time, OS) ~ explore.subtype, dat = PanGlioma.OS.data, power = 0.80, k = 0.145, RR = 2.21, alpha = 0.05)$ssize
 # nE  nC 
 # 34 231 

PanGlioma.DSS.data$explore.subtype <- ifelse(PanGlioma.DSS.data$explore.subtype == 'IDHsubmut.non.codel', 'E', 'C')
ssizeCT(formula = Surv(DSS.time, DSS) ~ explore.subtype, dat = PanGlioma.DSS.data, power = 0.80, k = 0.145, RR = 2.36, alpha = 0.05)$ssize
 # nE  nC 
 # 32 217 


# Random perturbation strategy
load('data/GliomaMolecularSubtype.RData')

glioma.molecular.subtype <- subset(glioma.molecular.subtype, 
 explore.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'), explore.subtype)

set.seed(1024)
glioma.random.subtype <- lapply(seq(1000), function(times){
 random.subtype <- glioma.molecular.subtype
 random.subtype$explore.subtype <- sample(random.subtype$explore.subtype, length(random.subtype$explore.subtype), replace=FALSE) 
 
 return(random.subtype)
})

save(glioma.random.subtype, file = 'result/3.Glioma_Maunally_defined_subtypes/nonCodelGliomaRandomSubtype.RData')

# random log-rank p value
load('data/PanGlioma.clinical.outcome.endpoints.RData')

random.logrank.pvalue <- lapply(seq(1000), function(times){
 require(survival)
 random.subtype <- glioma.random.subtype[[times]]
 
 os.data <- merge(random.subtype, PanGlioma.OS.data, by = 'row.names')
 dss.data <- merge(random.subtype, PanGlioma.DSS.data, by = 'row.names')
 
 os.pvalue <- pchisq(survdiff(Surv(OS.time, OS)~explore.subtype, data = os.data)$chisq, 1, lower.tail = FALSE)
 dss.pvalue <- pchisq(survdiff(Surv(DSS.time, DSS)~explore.subtype, data = dss.data)$chisq, 1, lower.tail = FALSE)
 
 return(c(os.pvalue, dss.pvalue))
})
random.logrank.pvalue <- do.call(rbind, random.logrank.pvalue)


# p value 
PanGlioma.OS.data <- merge(glioma.molecular.subtype, PanGlioma.OS.data, by = 'row.names')
os.logrank.pvalue <- pchisq(survdiff(Surv(OS.time, OS)~explore.subtype, data=PanGlioma.OS.data)$chisq, 1, lower.tail = FALSE)
sum(random.logrank.pvalue[, 1] < os.logrank.pvalue)/nrow(random.logrank.pvalue) # 0.041


PanGlioma.DSS.data <- merge(glioma.molecular.subtype, PanGlioma.DSS.data, by = 'row.names')
dss.logrank.pvalue <- pchisq(survdiff(Surv(DSS.time, DSS)~explore.subtype, data=PanGlioma.DSS.data)$chisq, 1, lower.tail = FALSE)
sum(random.logrank.pvalue[, 2] < dss.logrank.pvalue)/nrow(random.logrank.pvalue) # 0.031



