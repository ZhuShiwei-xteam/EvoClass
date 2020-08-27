setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')

# Statistics of clinical characteristics
load('data/gliomaClinicalData.RData')
rownames(pan.glio.mer.cli.data) <- pan.glio.mer.cli.data$bcr_patient_barcode

load(file = 'data/GliomaMolecularSubtype.RData')
pan.glio.mer.cli.data <- merge(pan.glio.mer.cli.data, glioma.molecular.subtype, by = 'row.names')

pan.glio.mer.cli.data$final.subtype <- factor(pan.glio.mer.cli.data$final.subtype, 
 levels = c('IDHmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel', 'IDHwt'))

# age_at_initial_pathologic_diagnosis

subtype.age.dis <- sapply(c('IDHmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel', 'IDHwt'), function(subtype){
 subtype.cli.data <- subset(pan.glio.mer.cli.data, final.subtype == subtype)
 age.median <- median(subtype.cli.data$age_at_initial_pathologic_diagnosis)
 age.mean <- round(mean(subtype.cli.data$age_at_initial_pathologic_diagnosis), 1)
 age.sd <- round(sd(subtype.cli.data$age_at_initial_pathologic_diagnosis), 1)
 age.range <- range(subtype.cli.data$age_at_initial_pathologic_diagnosis)
 
 return(c(age.median, age.mean, age.sd, age.range))
})

median(pan.glio.mer.cli.data$age_at_initial_pathologic_diagnosis)
round(mean(pan.glio.mer.cli.data$age_at_initial_pathologic_diagnosis), 1)
round(sd(pan.glio.mer.cli.data$age_at_initial_pathologic_diagnosis), 1)
range(pan.glio.mer.cli.data$age_at_initial_pathologic_diagnosis)

kruskal.test(age_at_initial_pathologic_diagnosis~final.subtype, data=pan.glio.mer.cli.data)$p.value # 1.794662e-73
wilcox.test(age_at_initial_pathologic_diagnosis~final.subtype, 
 data = subset(pan.glio.mer.cli.data, final.subtype %in% c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')))$p.value # 0.9924686


# gender, race, neoplasm_histologic_grade, tumor_location, Laterality
pan.glio.mer.cli.data$first_presenting_symptom[pan.glio.mer.cli.data$first_presenting_symptom %in% 
 c('Sensory Changes', 'Visual Changes')] <- 'Sensory or visual change'
 
pan.glio.mer.cli.data$tumor_location[pan.glio.mer.cli.data$tumor_locatio %in% 
 c('Posterior Fossa, Brain Stem', 'Posterior Fossa, Cerebellum',
 'Supratentorial, Not Otherwise Specified', 'Supratentorial, Occipital Lobe')] <- 'Other'      

pan.glio.mer.cli.data$race[pan.glio.mer.cli.data$race != 'WHITE'] <- 'Other'

clinical.fea <- c('gender', 'race', 'neoplasm_histologic_grade', 'tumor_location', 'laterality', 'first_presenting_symptom')
 
clinical.sig <- sapply(clinical.fea, function(feature){

 clinical.sta <- table(pan.glio.mer.cli.data[, c(feature, 'final.subtype')])
 
 across.subtype.pvalue <- chisq.test(clinical.sta)$p.value
 within.subtype.pvalue <- fisher.test(clinical.sta[, c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')])$p.value

 return(c(across.subtype.pvalue, within.subtype.pvalue))
}) 
 
# stastic
clinical.sta <- sapply(clinical.fea, function(feature){

 cli.sta <- table(pan.glio.mer.cli.data[, c(feature, 'final.subtype')])
 cli.sta <- cbind(cli.sta, all.patient = rowSums(cli.sta))
 
 cli.feq <- round(prop.table(cli.sta, 2), 2)
 
 cli.sta.feq <- sapply(rownames(cli.sta), function(index){
  sta.feq <- paste0(cli.sta[index, ], '(', cli.feq[index, ], ')')
  
  return(sta.feq)
 })
 rownames(cli.sta.feq) <- colnames(cli.sta)

 return(cli.sta.feq)
}) 


# Family history of cancer
# table(pan.glio.mer.cli.data[,c('family_history_of_primary_brain_tumor', 'family_history_of_cancer', 'final.subtype')])
# chisq.test(matrix(c(69,3,35,100, 8, 40, 13, 0, 8,26,3,26), 3))$p.value
# fisher.test(matrix(c(100, 8, 40, 13, 0, 8), 3))$p.value

