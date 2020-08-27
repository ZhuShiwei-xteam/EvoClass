setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")

suppressPackageStartupMessages(library(data.table))

# load mutation data
ori.pan.mut <- fread("data/OriginalData/mc3.v0.2.8.PUBLIC.maf", sep = '\t', header = TRUE, stringsAsFactors = FALSE, 
                         select = c('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 't_alt_count', 't_ref_count', 'Hugo_Symbol', 'Variant_Classification', 'HGVSp_Short', 'FILTER', 'NCALLERS'))
colnames(ori.pan.mut) <- c('Patient', 'Chr', 'Start_position', 'End_position', 'Reference', 'Alternate', 'Variant_freq', 'Ref_freq', 'Hugo_Symbol', 'Variant_Classification', 'HGVSp_Short', 'FILTER', 'NCALLERS')
# load clinical data
pan.cli.data <- read.csv("data/OriginalData/TCGA-CDR-SupplementalTableS1.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE)


ori.pan.mut$Patient <- substr(ori.pan.mut$Patient, 1, 15)
ori.pan.mut$cancer.type <- pan.cli.data[match(substr(ori.pan.mut$Patient, 1, 12), pan.cli.data$bcr_patient_barcode), 'type']
ori.pan.mut$sample.type <- substr(ori.pan.mut$Patient, 14, 15)
ori.pan.mut <- ori.pan.mut[-which(is.na(ori.pan.mut$cancer.type)), ] # 78 samples, 21945 mutations
ori.pan.mut$mutation.id <- paste(ori.pan.mut$Tumor_Sample_Barcode, ori.pan.mut$Chromosome, ori.pan.mut$Start_Position, ori.pan.mut$Reference_Allele, ori.pan.mut$Tumor_Seq_Allele2, sep = ':')

# retain mutations with specific mutation calling labels
ori.pan.mut.filter <- ori.pan.mut[ori.pan.mut$FILTER %in% c('PASS', 'wga', 'native_wga_mix'), ]

# identify hypermutated samples -------------------------------------------

ori.pan.pri.mut <- ori.pan.mut.filter
cancer.name <- unique(ori.pan.pri.mut$cancer.type)
hyper.mut.sample <- sapply(cancer.name, function(sin.cancer){

  if(sin.cancer %in% 'LUAD'){
    hyper.mut.threshold <- 1047
  }else if(sin.cancer %in% 'SKCM'){
    hyper.mut.threshold <- 2122
  }else if(sin.cancer %in% 'UCEC'){
    hyper.mut.threshold <- 2545
  }else{
    hyper.mut.threshold <- 1000
  }
  
  sin.cancer.mut.data <- ori.pan.pri.mut[ori.pan.pri.mut$cancer.type %in% sin.cancer, ]
  sin.cancer.sam.mutBur <- table(sin.cancer.mut.data$Patient)
  
  # hypermutator criteria 1
  sin.cancer.hyper.sam1 <- names(sin.cancer.sam.mutBur[sin.cancer.sam.mutBur > hyper.mut.threshold])

  # hypermutator criteria 2
  sin.cancer.IQR <- quantile(sin.cancer.sam.mutBur)[4] - quantile(sin.cancer.sam.mutBur)[2]
  sin.cancer.outlier.value <- (sin.cancer.IQR * 1.5) + quantile(sin.cancer.sam.mutBur)[4]
  sin.cancer.hyper.sam2 <- names(sin.cancer.sam.mutBur[sin.cancer.sam.mutBur > sin.cancer.outlier.value])
  
  # hypermutator determination
  sin.cancer.hyper.sam <- intersect(sin.cancer.hyper.sam1, sin.cancer.hyper.sam2)
  
})

# transform list to data frame
hyper.mut.sample.mat <- data.frame(Sample = NA, Cancer = NA)
for(i in 1:length(cancer.name)){
  sin.can.hyper <- hyper.mut.sample[[cancer.name[i]]]
  if(length(sin.can.hyper) > 0){
    sin.can.hyper.data <- data.frame(Sample = sin.can.hyper, Cancer = rep(cancer.name[i], length(sin.can.hyper)))
  }else{
    next
  }
  hyper.mut.sample.mat <- rbind(hyper.mut.sample.mat, sin.can.hyper.data)
}
hyper.mut.sample.mat <- hyper.mut.sample.mat[-1, ]
save(hyper.mut.sample.mat, file = "data/panCan_hypermutator_data.RData")

# produce mutation data for downstream analysis ---------------------------

# remove duplicated mutations
ori.pan.mut.filter.dedup <- ori.pan.mut.filter[!duplicated(ori.pan.mut.filter$mutation.id), ]

# only retain primary samples (except metastatic samples of SKCM) and non-hypermutated samples
pure.pan.mut.data <- ori.pan.mut.filter.dedup[ori.pan.mut.filter.dedup$sample.type %in% c('01', '03') | (ori.pan.mut.filter.dedup$sample.type %in% '06' & ori.pan.mut.filter.dedup$cancer.type %in% 'SKCM'), ]
pure.pan.mut.data <- pure.pan.mut.data[!pure.pan.mut.data$Tumor_Sample_Barcode %in% hyper.mut.sample.mat$Sample, ]
pure.glioma.mut.data <- pure.pan.mut.data[pure.pan.mut.data$cancer.type %in% c("GBM", "LGG"), ]
save(pure.glioma.mut.data, file = "data/PureGliomaMutData.RData")




