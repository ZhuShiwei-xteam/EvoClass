setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")

library(oncosign)
source('code/function/OncosignInputAndSubtypeExtract.R')
source('code/function/Oncosign.R')

# molecular classification of glioma with oncosign software

# load mutational clonality data
load('data/glioma.mut.ccf.RData')
mut.clonality.data <- com.data

silent.mut.label <- c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "Intron", "RNA", "Silent")
mut.clonality.data <- subset(mut.clonality.data, !(Variant_Classification %in% silent.mut.label))

# Combination of IDH1 and IDH2
mut.clonality.data$Hugo_Symbol[mut.clonality.data$Hugo_Symbol %in% c('IDH1', 'IDH2')] <- 'IDH'


# load driver genes with a mutation frequency of at least 2% of the sample
load('data/pan.glioma.cancer.gene.RData')
glioma.driver.gene <- pan.glioma.cancer.gene


# load('/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Resources/PureGliomaMutData.RData')
# # length(unique(pure.glioma.mut.data$Patient)) 896

# silent.mut.label <- c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "Intron", "RNA", "Silent")
# pure.glioma.mut.data <- subset(pure.glioma.mut.data, !(Variant_Classification %in% silent.mut.label))
# gene.mut.freq <- table(pure.glioma.mut.data$Hugo_Symbol)
# glioma.driver.gene <- intersect(glioma.driver.gene, names(gene.mut.freq)[gene.mut.freq > 896 *0.02])

# pan.glioma.cancer.gene <- glioma.driver.gene
# save(pan.glioma.cancer.gene, file = '/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Resources/pan.glioma.cancer.gene.RData')

glioma.driver.gene[glioma.driver.gene %in% c('IDH1', 'IDH2')] <- 'IDH'


# prepare input data for oncosign
input.file.path <- paste(getwd(), '/result/3.Glioma_Maunally_defined_subtypes/OncosignMolecularClassification',sep="")
OncosignInputFile(mut.clonality.data, glioma.driver.gene, input.file.path)


RunOncosign <- function(data.file.path, max.levels, min.alterations){

 setwd(data.file.path)
 set.seed(42)
 oncosign.res <- oncosign.sea(cases = 'sample.data.txt', gam = 'mut.data.txt', max.levels, min.alterations)
 
 return(oncosign.res)
}

# run oncosign

oncosign.results <- RunOncosign(input.file.path, max.levels = 3, min.alterations = 0.01)
save(oncosign.results, file = 'OncosignResults.RData')


molecular.subtype <- lapply(1:3, function(level) SubtypeExtract(oncosign.results$classes, level))
molecular.subtype <- Reduce(function(x, y) merge(x, y, by = 'patient'), molecular.subtype)

save(molecular.subtype, file = 'OncosignMolecularSubtype.RData')




