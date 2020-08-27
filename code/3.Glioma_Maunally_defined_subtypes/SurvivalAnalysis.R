setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")

library('mclust')
source('code/3.Glioma_Maunally_defined_subtypes/SystemSurvivalPlot.R')
source('code/3.Glioma_Maunally_defined_subtypes/Cox.function.R')

# load clinical data
load('data/gliomaClinicalData.RData')
rownames(pan.glio.mer.cli.data) <- pan.glio.mer.cli.data$bcr_patient_barcode

# load oncosign molecular subtype
load(file = 'result/3.Glioma_Maunally_defined_subtypes/OncosignMolecularClassification/OncosignMolecularSubtype.RData')
load(file = 'result/3.Glioma_Maunally_defined_subtypes/OncosignMolecularClassification/OncosignResults.RData')
rownames(molecular.subtype) <- molecular.subtype$patient

# merge data
clinical.annotation <- merge(molecular.subtype, pan.glio.mer.cli.data, by = 'row.names')

# adjusted Rand index
adjustedRandIndex(clinical.annotation$level.1.subtype, clinical.annotation$IDH.codel.subtype) # [1] 0.6348182


table(subset(clinical.annotation, level.3.subtype == 'OSC.3.5.3')[, c('IDH.codel.subtype', 'level.3.subtype')]) # OSC.3.5.3(26/26, 100%)
table(subset(clinical.annotation, level.3.subtype == 'OSC.3.5.13')[, c('IDH.codel.subtype', 'level.3.subtype')])# OSC.3.5.3(38/49, 78%)

# survival analysis
load('data/GliomaMolecularSubtype.RData')
load('data/PanGlioma.clinical.outcome.endpoints.RData')


# IDH clonal VS subclonal mutation
setwd(paste(getwd(), "/result/3.Glioma_Maunally_defined_subtypes", sep=""))
	
SystemSurvivalPlot(PanGlioma.OS.data, PanGlioma.DSS.data, PanGlioma.PFI.data, PanGlioma.DFI.data, 
 glioma.molecular.subtype[,'IDH.mut.status', drop = FALSE], risk.table.index = TRUE, 'IDHmut.subtype')

SystemSurvivalPlot(PanGlioma.OS.data, PanGlioma.DSS.data, PanGlioma.PFI.data, PanGlioma.DFI.data, 
 subset(glioma.molecular.subtype, grade == 'grade_II-III', 'IDH.mut.status'), risk.table.index = TRUE, 'LGGIDHmut.subtype')

SystemSurvivalPlot(PanGlioma.OS.data, PanGlioma.DSS.data, PanGlioma.PFI.data, PanGlioma.DFI.data, 
 subset(glioma.molecular.subtype, grade == 'grade_IV', 'IDH.mut.status'), risk.table.index = TRUE, 'GBMIDHmut.subtype')


SystemSurvivalPlot(PanGlioma.OS.data, PanGlioma.DSS.data, PanGlioma.PFI.data, PanGlioma.DFI.data, 
 subset(glioma.molecular.subtype, explore.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'), explore.subtype), 
 risk.table.index = TRUE, 'IDHmut.non.codel.subtype')

SystemSurvivalPlot(PanGlioma.OS.data, PanGlioma.DSS.data, PanGlioma.PFI.data, PanGlioma.DFI.data, 
 subset(glioma.molecular.subtype, 
 explore.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel') & grade == 'grade_II-III', explore.subtype), 
 risk.table.index = TRUE, 'LGGIDHmut.non.codel.subtype')

SystemSurvivalPlot(PanGlioma.OS.data, PanGlioma.DSS.data, PanGlioma.PFI.data, PanGlioma.DFI.data, 
 subset(glioma.molecular.subtype, 
 explore.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel') & grade == 'grade_IV', explore.subtype), 
 risk.table.index = TRUE, 'GBMIDHmut.non.codel.subtype')


SystemSurvivalPlot(PanGlioma.OS.data, PanGlioma.DSS.data, PanGlioma.PFI.data, PanGlioma.DFI.data, 
 subset(glioma.molecular.subtype, explore.subtype %in% c('IDHcolmut.codel', 'IDHsubmut.codel'), explore.subtype), 
 risk.table.index = TRUE, 'IDHmut.codel.subtype')



# Multivariate cox proportional hazards analysis

# overall survival
PanGlioma.OS.data <- merge(glioma.molecular.subtype, PanGlioma.OS.data, by = 'row.names')
PanGlioma.OS.data <- PanGlioma.OS.data[, c('bcr_patient_barcode', 'OS.time', 'OS', 
 'age_at_initial_pathologic_diagnosis', 'grade', 'explore.subtype'), drop = FALSE] 


PanGlioma.OS.data <- subset(PanGlioma.OS.data, explore.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))
PanGlioma.OS.data$explore.subtype <- factor(PanGlioma.OS.data$explore.subtype, 
 levels = c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))


PanGlioma.OS.data$age_at_initial_pathologic_diagnosis <- as.numeric(PanGlioma.OS.data$age_at_initial_pathologic_diagnosis)
PanGlioma.OS.data$grade <- factor(PanGlioma.OS.data$grade, levels = c('grade_II-III', 'grade_IV')) 

os.cox.results <- Cox.function(time = PanGlioma.OS.data$OS.time, event = PanGlioma.OS.data$OS, 
 clinical.data = PanGlioma.OS.data, clinical.variate = NULL)


# disease specific survival
PanGlioma.DSS.data <- merge(glioma.molecular.subtype, PanGlioma.DSS.data, by = 'row.names')
PanGlioma.DSS.data <- PanGlioma.DSS.data[, c('bcr_patient_barcode', 'DSS.time', 'DSS', 
 'age_at_initial_pathologic_diagnosis', 'grade', 'explore.subtype'), drop = FALSE] 


PanGlioma.DSS.data <- subset(PanGlioma.DSS.data, explore.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))
PanGlioma.DSS.data$explore.subtype <- factor(PanGlioma.DSS.data$explore.subtype, 
 levels = c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))


PanGlioma.DSS.data$age_at_initial_pathologic_diagnosis <- as.numeric(PanGlioma.DSS.data$age_at_initial_pathologic_diagnosis)
PanGlioma.DSS.data$grade <- factor(PanGlioma.DSS.data$grade, levels = c('grade_II-III', 'grade_IV')) 

dss.cox.results <- Cox.function(time = PanGlioma.DSS.data$DSS.time, event = PanGlioma.DSS.data$DSS, 
 clinical.data = PanGlioma.DSS.data, clinical.variate = NULL)

save(os.cox.results, dss.cox.results, file = "cox.result.RData")