setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")

# load mutation clonaltiy data
load('data/glioma.mut.ccf.RData')
mut.clonality.data <- com.data

silent.mut.label <- c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "Intron", "RNA", "Silent")
mut.clonality.data <- subset(mut.clonality.data, !(Variant_Classification %in% silent.mut.label))
mut.clonality.data$patient <- substr(mut.clonality.data$patient, 1, 12)

# Combination of IDH1 and IDH2
mut.clonality.data$Hugo_Symbol[mut.clonality.data$Hugo_Symbol %in% c('IDH1', 'IDH2')] <- 'IDH'

IDHmut.subclonal <- subset(mut.clonality.data, (Hugo_Symbol == 'IDH')&(CI95.timing == 'Subclonal'))$patient
IDHmut.clonal <- subset(mut.clonality.data, (Hugo_Symbol == 'IDH')&(CI95.timing == 'Clonal'))$patient


# A core sample set of 876 primary diffuse gliomas
# whole exon sequencing samples
load('data/PureGliomaMutData.RData')
wes.patients <- unique(substr(pure.glioma.mut.data$Patient, 1, 12))

load('data/panCanAtlas.abs.seg.RData')
panCanAtlas.abs.seg <- panCanAtlas.abs.seg[substr(panCanAtlas.abs.seg$SampleID, 14, 15) %in% '01', ]
cnv.patients <- unique(substr(panCanAtlas.abs.seg$SampleID, 1, 12))

absolute.purity <- read.table(file = 'data/OriginalData/TCGA_mastercalls.abs_tables_JSedit.fixed.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
absolute.purity <- absolute.purity[substr(absolute.purity$array, 14, 15) %in% '01', ]
purity.patients <- substr(absolute.purity$array, 1, 12)

core.patients <- Reduce(intersect, list(wes.patients, cnv.patients, purity.patients))

# load clinical data
load('data/gliomaClinicalData.RData')
rownames(pan.glio.mer.cli.data) <- pan.glio.mer.cli.data$bcr_patient_barcode
pan.glio.mer.cli.data <- pan.glio.mer.cli.data[core.patients, ]

IDHmut.codel <- subset(pan.glio.mer.cli.data, IDH.codel.subtype == 'IDHmut-codel')$bcr_patient_barcode
IDHmut.non.codel <- subset(pan.glio.mer.cli.data, IDH.codel.subtype == 'IDHmut-non-codel')$bcr_patient_barcode
IDHwt <- subset(pan.glio.mer.cli.data, IDH.codel.subtype == 'IDHwt')$bcr_patient_barcode


# manually defined molecular subtype
IDHcolmut.codel <- intersect(IDHmut.clonal, IDHmut.codel)
IDHsubmut.codel <- intersect(IDHmut.subclonal, IDHmut.codel)

IDHcolmut.non.codel <- intersect(IDHmut.clonal, IDHmut.non.codel)
IDHsubmut.non.codel <- intersect(IDHmut.subclonal, IDHmut.non.codel)

glioma.molecular.subtype <- data.frame(IDH.mut.status = rep.int(x = c('IDHcolmut', 'IDHsubmut', 'IDHcolmut', 'IDHsubmut', 'IDHwt'), 
 times = lapply(list(IDHcolmut.codel, IDHsubmut.codel, IDHcolmut.non.codel, IDHsubmut.non.codel, IDHwt), length)), 
 stringsAsFactors = FALSE)
 
rownames(glioma.molecular.subtype) <- c(IDHcolmut.codel, IDHsubmut.codel, IDHcolmut.non.codel, IDHsubmut.non.codel, IDHwt)
 

glioma.molecular.subtype$explore.subtype <- rep.int(
 x = c('IDHcolmut.codel', 'IDHsubmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel', 'IDHwt'), 
 times = lapply(list(IDHcolmut.codel, IDHsubmut.codel, IDHcolmut.non.codel, IDHsubmut.non.codel, IDHwt), length))


glioma.molecular.subtype$final.subtype <- rep.int(
 x = c('IDHmut.codel', 'IDHmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel', 'IDHwt'), 
 times = lapply(list(IDHcolmut.codel, IDHsubmut.codel, IDHcolmut.non.codel, IDHsubmut.non.codel, IDHwt), length))

glioma.molecular.subtype$grade <- pan.glio.mer.cli.data[rownames(glioma.molecular.subtype), 'Grade.rough']

save(glioma.molecular.subtype, file = 'data/GliomaMolecularSubtype.RData')
