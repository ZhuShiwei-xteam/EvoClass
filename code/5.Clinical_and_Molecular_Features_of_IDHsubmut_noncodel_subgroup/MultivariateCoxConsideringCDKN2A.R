setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')

source('code/3.Glioma_Maunally_defined_subtypes/Cox.function.R')
# CDKN2A Deletion
gene.cnv.alt <- read.table(file = 'data/OriginalData/all_thresholded.by_genes.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(gene.cnv.alt) <- gene.cnv.alt$Gene.Symbol
gene.cnv.alt <- gene.cnv.alt[, -c(1:3)]
colnames(gene.cnv.alt) <- gsub('\\.', '-', substr(colnames(gene.cnv.alt), 1, 12))

# Multivariate cox proportional hazards analysis
load('data/GliomaMolecularSubtype.RData')
load('data/PanGlioma.clinical.outcome.endpoints.RData')

# overall survival
PanGlioma.OS.data <- merge(glioma.molecular.subtype, PanGlioma.OS.data, by = 'row.names')
PanGlioma.OS.data$CDKN2A <- as.character(t(gene.cnv.alt['CDKN2A', PanGlioma.OS.data$Row.names]))

PanGlioma.OS.data <- PanGlioma.OS.data[, c('bcr_patient_barcode', 'OS.time', 'OS', 
 'age_at_initial_pathologic_diagnosis', 'grade', 'CDKN2A', 'explore.subtype'), drop = FALSE] 


PanGlioma.OS.data <- subset(PanGlioma.OS.data, explore.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))
PanGlioma.OS.data$explore.subtype <- factor(PanGlioma.OS.data$explore.subtype, 
 levels = c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))


PanGlioma.OS.data$age_at_initial_pathologic_diagnosis <- as.numeric(PanGlioma.OS.data$age_at_initial_pathologic_diagnosis)
PanGlioma.OS.data$grade <- factor(PanGlioma.OS.data$grade, levels = c('grade_II-III', 'grade_IV')) 
PanGlioma.OS.data$CDKN2A <- ifelse(PanGlioma.OS.data$CDKN2A %in% c('-2'), TRUE, FALSE)


os.cox.results <- Cox.function(time = PanGlioma.OS.data$OS.time, event = PanGlioma.OS.data$OS, 
 clinical.data = PanGlioma.OS.data, clinical.variate = NULL)


# disease specific survival
PanGlioma.DSS.data <- merge(glioma.molecular.subtype, PanGlioma.DSS.data, by = 'row.names')
PanGlioma.DSS.data$CDKN2A <- as.character(t(gene.cnv.alt['CDKN2A', PanGlioma.DSS.data$Row.names]))

PanGlioma.DSS.data <- PanGlioma.DSS.data[, c('bcr_patient_barcode', 'DSS.time', 'DSS', 
 'age_at_initial_pathologic_diagnosis', 'grade', 'CDKN2A', 'explore.subtype'), drop = FALSE] 


PanGlioma.DSS.data <- subset(PanGlioma.DSS.data, explore.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))
PanGlioma.DSS.data$explore.subtype <- factor(PanGlioma.DSS.data$explore.subtype, 
 levels = c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'))


PanGlioma.DSS.data$age_at_initial_pathologic_diagnosis <- as.numeric(PanGlioma.DSS.data$age_at_initial_pathologic_diagnosis)
PanGlioma.DSS.data$grade <- factor(PanGlioma.DSS.data$grade, levels = c('grade_II-III', 'grade_IV')) 
PanGlioma.DSS.data$CDKN2A <- ifelse(PanGlioma.DSS.data$CDKN2A %in% c('-2'), TRUE, FALSE)


dss.cox.results <- Cox.function(time = PanGlioma.DSS.data$DSS.time, event = PanGlioma.DSS.data$DSS, 
 clinical.data = PanGlioma.DSS.data, clinical.variate = NULL)

