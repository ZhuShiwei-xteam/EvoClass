setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')


source('code/function/PycloneFunctionModule.R')
load('data/PureGliomaMutData.RData')
load('data/panCanAtlas.abs.seg.RData')
load('data/panCanAtlas.sex.RData') 

# remove indel alterations
panCanAtlas.maf <- as.data.frame(pure.glioma.mut.data)
panCanAtlas.maf <- panCanAtlas.maf[-(union(grep('-', panCanAtlas.maf$Alternate), grep('-', panCanAtlas.maf$Reference))), ]
panCanAtlas.maf <- panCanAtlas.maf[(nchar(panCanAtlas.maf$Alternate) == 1) & (nchar(panCanAtlas.maf$Reference) == 1), ]


panCanAtlas.sex <- split(panCanAtlas.sex, f = panCanAtlas.sex$type)
output.file.path <- 'result/4.Intratumoral_Heterogeneity_Analysis/PycloneITH/'

PycloneInputData('LGG', panCanAtlas.sex, panCanAtlas.maf, panCanAtlas.abs.seg, output.file.path)
PycloneInputData('GBM', panCanAtlas.sex, panCanAtlas.maf, panCanAtlas.abs.seg, output.file.path)

# determined the tumor subclonal composition with an independent algorithm Pyclone
tumor.purity.data <- panCanAtlas.abs.seg[, c('SampleID', 'Aberrant Cell Fraction')]
tumor.purity.data <- tumor.purity.data[!duplicated(tumor.purity.data), ]
rownames(tumor.purity.data ) <- tumor.purity.data $SampleID

setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')
input.file.path <- paste(getwd(), '/result/4.Intratumoral_Heterogeneity_Analysis/PycloneITH/TCGA_LGG/', sep="")
try(RunPyClone(input.file.path, tumor.purity.data, s.index = 1, e.index = 525), silent = FALSE)


input.file.path <- 'result/4.Intratumoral_Heterogeneity_Analysis/PycloneITH/TCGA_GBM/'
try(RunPyClone(input.file.path, tumor.purity.data, s.index = 301, e.index = 400), silent = FALSE)



setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')
# The number of detectable subclonal populations
file.path <- '/result/4.Intratumoral_Heterogeneity_Analysis/PycloneITH/TCGA_LGG/'
lgg.clonal.num <- PyCloneClonalPopulationNumber(file.path)
file.path <- '/result/4.Intratumoral_Heterogeneity_Analysis/PycloneITH/TCGA_GBM/'
gbm.clonal.num <- PyCloneClonalPopulationNumber(file.path)

clonal.num <- c(lgg.clonal.num, gbm.clonal.num)
clonal.num <- as.data.frame(clonal.num)
rownames(clonal.num) <- substr(rownames(clonal.num), 1, 12)

load("data/GliomaMolecularSubtype.RData")
glioma.molecular.subtype <- merge(glioma.molecular.subtype, clonal.num, by = 'row.names')

non.codel.clonal.pvalue <- wilcox.test(clonal.num~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
	c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value # p-value = 0.007716565
	
codel.clonal.pvalue <- wilcox.test(clonal.num~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
	c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value # p-value = 0.0008332387


# plot
ClonalPopulationBoxplot <- function(clonal.population.data, file.name, file.path){
 library('ggplot2')
 
 box.plot <- ggplot(data = clonal.population.data, aes(explore.subtype, clonal.num)) +
  geom_boxplot(aes(fill = explore.subtype) , outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8), limits=c(0, 8)) + xlab('Molecular subtypes') + ylab('Subclonal population') + 
  guides(fill = 'none') + theme(axis.text.x  = element_text(angle=30, vjust=0.5))
 
 ggsave(filename = file.name, plot = box.plot, path = file.path)
}	


glioma.molecular.subtype <- subset(glioma.molecular.subtype, explore.subtype != 'IDHwt')
glioma.molecular.subtype$explore.subtype <- factor(glioma.molecular.subtype$explore.subtype, 
	levels = c('IDHcolmut.codel', 'IDHsubmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel'))

ClonalPopulationBoxplot(glioma.molecular.subtype, 
	'ClonalPopulationBoxplot.pdf', 'result/4.Intratumoral_Heterogeneity_Analysis/PycloneITH/')


# Chi-square test
glioma.molecular.subtype$clonal.num[glioma.molecular.subtype$clonal.num >= 4] <- 4

# chisq.test(table(glioma.molecular.subtype[, c('explore.subtype','clonal.num')])
 # [c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'), ]) # p-value = 0.05641
 
# chisq.test(table(glioma.molecular.subtype[, c('explore.subtype','clonal.num')])
 # [c('IDHcolmut.codel', 'IDHsubmut.codel'),]) # p-value = 0.00277


glioma.molecular.subtype$clonal.num <- factor(glioma.molecular.subtype$clonal.num, levels = c(4, 3, 2, 1))
p <- ggplot(data=glioma.molecular.subtype, mapping = aes(explore.subtype))+geom_bar(aes(fill=clonal.num), position="fill") + 
 xlab('Molecular subtypes') + ylab('Percentage') + theme(axis.text.x  = element_text(angle=30, vjust=0.5))
ggsave(filename = 'ClonalPopulationStackplot.pdf', plot = p, 
 path = 'result/4.Intratumoral_Heterogeneity_Analysis/PycloneITH/')


