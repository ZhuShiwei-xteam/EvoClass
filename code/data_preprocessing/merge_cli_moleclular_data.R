setwd("/pub5/xiaoyun/Jobs/J22/temp/Github/")


suppressPackageStartupMessages(library(tidyverse))

# load clinical characteristic data
pan.cli.fea.data <- read.csv("data/OriginalData/clinical_PANCAN_patient_with_followup.tsv", header = TRUE, stringsAsFactors = FALSE, sep = '\t', na.strings = c("", "[Not Available]", "[Not Applicable]", '[Unknown]', '[Discrepancy]', '[Not Evaluated]', 'Not listed in Medical Record'))
pan.glio.cli.fea.data <- pan.cli.fea.data[pan.cli.fea.data$acronym %in% c('LGG', 'GBM'), ] # 1111

# 36 common clinical characteristics
analyzed.cli.fea <- c('age_at_initial_pathologic_diagnosis', 'neoplasm_histologic_grade', 'karnofsky_performance_score', 'gender', 'histological_type', 
                      'new_tumor_event_after_initial_treatment', 'radiation_therapy', 'race', 'person_neoplasm_cancer_status', 'tissue_source_site', 'primary_therapy_outcome_success', 
                      'laterality', 'targeted_molecular_therapy', 'additional_pharmaceutical_therapy', 'additional_radiation_therapy', 'days_to_new_tumor_event_after_initial_treatment',
                      'family_history_of_cancer', 'supratentorial_localization', 'mold_or_dust_allergy_history', 'motor_movement_changes', 
                      'asthma_history','tumor_location','ldh1_mutation_test_method','first_presenting_symptom','animal_insect_allergy_history',
                      'first_presenting_symptom_longest_duration','visual_changes','sensory_changes','food_allergy_history','mental_status_changes','hay_fever_history',
                      'preoperative_corticosteroids','preoperative_antiseizure_meds','family_history_of_primary_brain_tumor', 'headache_history', 'eczema_history')

pan.glio.cli.fea.data <- pan.glio.cli.fea.data[, c("bcr_patient_barcode", analyzed.cli.fea)]

# load molecular features data
pan.glio.cli.data <- read.csv("data/OriginalData/lgggbm_tcga_pub_clinical_data.tsv", header = TRUE, stringsAsFactors = FALSE, sep = '\t', na.strings = "")

# intersection of clinical characteristic and molecular features data

cli.mol.inter.sam <- intersect(pan.glio.cli.fea.data$bcr_patient_barcode, pan.glio.cli.data$Patient.ID) # 1111
pan.glio.cli.data <- pan.glio.cli.data[match(pan.glio.cli.fea.data$bcr_patient_barcode, pan.glio.cli.data$Patient.ID), ]

# 13 salient molecular features
pan.glio.cli.data <- pan.glio.cli.data[match(pan.glio.cli.fea.data$bcr_patient_barcode, pan.glio.cli.data$Patient.ID), c('IDH.status', 
                                                                                                                         'IDH.codel.subtype', 'MGMT.promoter.status', 'ATRX.status', 'TERT.promoter.status', 'Transcriptome.Subtype', 
                                                                                                                         'Chr.19.20.co.gain', 'Chr.7.gain.Chr.10.loss', 'IDH.1P10Q.Subtype', 'Supervised.DNA.Methylation.Cluster', 'TERT.expression..log2.', 
                                                                                                                         'TERT.expression.status', 'Telomere.Maintenance')]

## merge

pan.glio.mer.cli.data <- cbind(pan.glio.cli.fea.data, pan.glio.cli.data)

# load mutation data
load("data/PureGliomaMutData.RData")

# revise IDH status and IDH codel subtype based on mutation data

IDH_from_mut <- pure.glioma.mut.data %>%
  mutate(IDH_logic = ifelse(Hugo_Symbol %in% c("IDH1", "IDH2") & !Variant_Classification %in% c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "Intron", "RNA", "Silent"), 1, 0)) %>%
  group_by(Patient) %>%
  summarise(IDH_fromMut = sum(IDH_logic)) %>%
  mutate(IDH_fromMut = ifelse(IDH_fromMut == 1, "Mutant", "WT"), bcr_patient_barcode = str_sub(Patient, 1, 12))

pan.glio.mer.cli.data <- pan.glio.mer.cli.data %>%
  left_join(IDH_from_mut, by = "bcr_patient_barcode") %>%
  mutate(IDH.status = ifelse(is.na(IDH_fromMut), IDH.status, IDH_fromMut), 
         IDH.temp = ifelse(IDH.status == "Mutant", "IDHmut", "IDHwt"), 
         IDH.codel.subtype = str_c(IDH.temp, IDH.1P10Q.Subtype, sep = "-"), 
         IDH.codel.subtype = ifelse(IDH.codel.subtype %in% "IDHwt-non-codel", "IDHwt", IDH.codel.subtype), 
         Grade.rough = ifelse(str_detect(histological_type, "GBM"), "grade_IV", "grade_II-III"), 
         WHO_subtype = str_c(IDH.codel.subtype, Grade.rough, sep = "_")) %>%
  select(-(Patient:IDH.temp))

save(pan.glio.mer.cli.data, file = "data/gliomaClinicalData.RData")








