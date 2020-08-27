setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')

library(data.table)

### load Pan-Cancer mut data
col.names <- c("Tumor_Sample_Barcode",  "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
	"t_alt_count", "t_ref_count", "Hugo_Symbol","Variant_Classification", "HGVSp_Short")
panCanAtlas.maf <- fread(input = 'data/OriginalData/mc3.v0.2.8.PUBLIC.maf.gz', sep = '\t', header = TRUE, select = col.names)
panCanAtlas.maf <- as.data.frame(panCanAtlas.maf)

panCanAtlas.maf$Tumor_Sample_Barcode <- substr(x = panCanAtlas.maf$Tumor_Sample_Barcode, start = 1, stop = 15)
panCanAtlas.maf$Chromosome[panCanAtlas.maf$Chromosome == 'X'] <- 23
panCanAtlas.maf$Chromosome[panCanAtlas.maf$Chromosome == 'Y'] <- 24
colnames(panCanAtlas.maf) <- c('Patient', 'Chr', 'Start_position', 'End_position', 'Reference', 
									'Alternate', 'Variant_freq', 'Ref_freq', 'Hugo_Symbol', 'Variant_Classification', "HGVSp_Short")

### load Pan-Cancer absolute copy number data
col.names <- c('Sample', 'Chromosome', 'Start', 'End', 'Num_Probes', 'Modal_Total_CN', 'Modal_HSCN_1', 'Modal_HSCN_2')
panCanAtlas.abs.seg <- fread(input = 'data/OriginalData/TCGA_mastercalls.abs_segtabs.fixed.txt', sep = '\t', header = TRUE, select = col.names)
panCanAtlas.abs.seg <- as.data.frame(panCanAtlas.abs.seg)
colnames(panCanAtlas.abs.seg) <- c('SampleID', 'Chr', 'Start', 'End', 'nProbes', 'cn', 'nA', 'nB')

### load Pan-Cancer absolute tumor purity data
col.names <- c('array', 'ploidy', 'purity')
panCanAtlas.purity <- fread(input = 'data/OriginalData/TCGA_mastercalls.abs_tables_JSedit.fixed.txt', sep = '\t', header = TRUE, select = col.names)
panCanAtlas.purity <- as.data.frame(panCanAtlas.purity)
colnames(panCanAtlas.purity) <- c('SampleID', 'Ploidy', 'Aberrant Cell Fraction')

### Merge absolute copy number and purity data
# panCanAtlas.abs.seg <- merge(x = panCanAtlas.abs.seg, y = panCanAtlas.purity, by = 'SampleID')
patients <- intersect(unique(panCanAtlas.abs.seg$SampleID), panCanAtlas.purity$SampleID)
panCanAtlas.abs.seg <- panCanAtlas.abs.seg[panCanAtlas.abs.seg$SampleID%in%patients, ]
panCanAtlas.abs.seg <- cbind(panCanAtlas.abs.seg, panCanAtlas.purity[match(panCanAtlas.abs.seg$SampleID, panCanAtlas.purity$SampleID), 
c('Ploidy', 'Aberrant Cell Fraction')])

### load Pan-Cancer clinical data
panCanAtlas.clinical.data <- read.csv(file = 'data/OriginalData/TCGA-CDR-SupplementalTableS1.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
panCanAtlas.clinical.data <- panCanAtlas.clinical.data[, c('bcr_patient_barcode', 'gender', 'type')]
panCanAtlas.clinical.data$gender <- tolower(panCanAtlas.clinical.data$gender)


### sex file
patients <- unique(panCanAtlas.maf$Patient)
bcr_patient_barcode <- substr(x = patients, start = 1, stop = 12)
patients <- data.frame(Patient = patients, bcr_patient_barcode = bcr_patient_barcode, stringsAsFactors = FALSE)
panCanAtlas.sex <- merge(x = patients, y = panCanAtlas.clinical.data, by = 'bcr_patient_barcode')
colnames(panCanAtlas.sex) <- c('barcode', 'sample', 'gender', 'type')

save(panCanAtlas.maf, file = 'data/panCanAtlas.maf.RData')
save(panCanAtlas.abs.seg, file = 'data/panCanAtlas.abs.seg.RData')
save(panCanAtlas.sex, file = 'data/panCanAtlas.sex.RData')