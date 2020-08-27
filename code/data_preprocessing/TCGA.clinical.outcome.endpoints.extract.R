setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")

	
# load pancancer clinical outcome data
PanCancer.clinical.outcome.endpoints <- read.csv(file = 'data/OriginalData/TCGA-CDR-SupplementalTableS1.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# extract diffuse gliomas' clinial outcome
PanGlioma.clinical.outcome.endpoints <- PanCancer.clinical.outcome.endpoints[PanCancer.clinical.outcome.endpoints$type %in% c('LGG', 'GBM'), ]
rownames(PanGlioma.clinical.outcome.endpoints) <- PanGlioma.clinical.outcome.endpoints$bcr_patient_barcode

# OS remove NA
PanGlioma.OS.data <- PanGlioma.clinical.outcome.endpoints[PanGlioma.clinical.outcome.endpoints$OS != '#N/A', ]
PanGlioma.OS.data <- PanGlioma.OS.data[PanGlioma.OS.data$OS.time != '#N/A', ]
PanGlioma.OS.data$OS <- as.numeric(PanGlioma.OS.data$OS)
PanGlioma.OS.data$OS.time <- as.numeric(PanGlioma.OS.data$OS.time)
lgg.OS.data <- PanGlioma.OS.data[PanGlioma.OS.data$type == 'LGG', ]
gbm.OS.data <- PanGlioma.OS.data[PanGlioma.OS.data$type == 'GBM', ]

# DSS remove NA
PanGlioma.DSS.data <- PanGlioma.clinical.outcome.endpoints[PanGlioma.clinical.outcome.endpoints$DSS != '#N/A', ]
PanGlioma.DSS.data <- PanGlioma.DSS.data[PanGlioma.DSS.data$DSS.time != '#N/A', ]
PanGlioma.DSS.data$DSS <- as.numeric(PanGlioma.DSS.data$DSS)
PanGlioma.DSS.data$DSS.time <- as.numeric(PanGlioma.DSS.data$DSS.time)
lgg.DSS.data <- PanGlioma.DSS.data[PanGlioma.DSS.data$type == 'LGG', ]
gbm.DSS.data <- PanGlioma.DSS.data[PanGlioma.DSS.data$type == 'GBM', ]

# DFI remove NA
PanGlioma.DFI.data <- PanGlioma.clinical.outcome.endpoints[PanGlioma.clinical.outcome.endpoints$DFI != '#N/A', ]
PanGlioma.DFI.data <- PanGlioma.DFI.data[PanGlioma.DFI.data$DFI.time != '#N/A', ]
PanGlioma.DFI.data$DFI <- as.numeric(PanGlioma.DFI.data$DFI)
PanGlioma.DFI.data$DFI.time <- as.numeric(PanGlioma.DFI.data$DFI.time)
lgg.DFI.data <- PanGlioma.DFI.data[PanGlioma.DFI.data$type == 'LGG', ]
gbm.DFI.data <- PanGlioma.DFI.data[PanGlioma.DFI.data$type == 'GBM', ]

# PFI remove NA
PanGlioma.PFI.data <- PanGlioma.clinical.outcome.endpoints[PanGlioma.clinical.outcome.endpoints$PFI != '#N/A', ]
PanGlioma.PFI.data <- PanGlioma.PFI.data[PanGlioma.PFI.data$PFI.time != '#N/A', ]
PanGlioma.PFI.data$PFI <- as.numeric(PanGlioma.PFI.data$PFI)
PanGlioma.PFI.data$PFI.time <- as.numeric(PanGlioma.PFI.data$PFI.time)
lgg.PFI.data <- PanGlioma.PFI.data[PanGlioma.PFI.data$type == 'LGG', ]
gbm.PFI.data <- PanGlioma.PFI.data[PanGlioma.PFI.data$type == 'GBM', ]


save(PanGlioma.OS.data, lgg.OS.data, gbm.OS.data, PanGlioma.DSS.data, lgg.DSS.data, gbm.DSS.data,
PanGlioma.DFI.data, lgg.DFI.data, gbm.DFI.data, PanGlioma.PFI.data, lgg.PFI.data, gbm.PFI.data, file = 'data/PanGlioma.clinical.outcome.endpoints.RData')






