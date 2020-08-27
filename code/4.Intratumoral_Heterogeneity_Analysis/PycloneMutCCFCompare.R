setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')
	
 library('ggplot2')
# the cellular prevalence (proportion of cancer cells harboring a mutation) of mutations
source('code/function/PycloneFunctionModule.R')

file.path <- '/result/4.Intratumoral_Heterogeneity_Analysis/PycloneITH/TCGA_LGG/' # 这个数据是什么？
lgg.mut.ccf <- PyCloneMutCCF(file.path)
file.path <- '/result/4.Intratumoral_Heterogeneity_Analysis/PycloneITH/TCGA_GBM/'
gbm.mut.ccf <- PyCloneMutCCF(file.path)
glioma.mut.ccf <- rbind(lgg.mut.ccf, gbm.mut.ccf)

glioma.mut.ccf$sample_id <- substr(glioma.mut.ccf$sample_id , 1, 12)


load("data/GliomaMolecularSubtype.RData")
mut.clonality.data <- cbind(glioma.mut.ccf, glioma.molecular.subtype[glioma.mut.ccf$sample_id, 'explore.subtype', FALSE])

mut.clonality.data <- subset(mut.clonality.data, explore.subtype %in% 
 c('IDHcolmut.codel', 'IDHsubmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel'))
mut.clonality.data$explore.subtype <- factor(mut.clonality.data$explore.subtype, 
	levels = c('IDHcolmut.codel', 'IDHsubmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel'))


non.codel.ccf.pvalue <- wilcox.test(cellular_prevalence~explore.subtype, subset(mut.clonality.data, explore.subtype %in% 
	c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value # 0.02475847

codel.ccf.pvalue <- wilcox.test(cellular_prevalence~explore.subtype, subset(mut.clonality.data, explore.subtype %in% 
	c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value # 1.465032e-06


# Cumulative distribution of the somatic mutations by CCF
ccf.ecdf.plot <- function(ccfData, fileName, filePath){
	
 ccf.cumu <- ggplot(data = ccfData, aes(cellular_prevalence)) + 
  stat_ecdf(aes(colour = explore.subtype), geom="line", size = 1) + 
  xlab('Cancer cell fraction') + ylab('Fraction of mutations') # + guides(colour = 'none')

 ggsave(filename = fileName, plot = ccf.cumu, path = filePath)
}
	

ccf.ecdf.plot(mut.clonality.data, 'PycloneCCFCumulativeDensity.pdf', '/result/4.Intratumoral_Heterogeneity_Analysis/PycloneITH/')

