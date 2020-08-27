setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')


load("data/GliomaMolecularSubtype.RData")
load("data/glioma.mut.ccf.RData")
com.data$patient <- substr(com.data$patient, 1, 12)

mut.clonality.data <- com.data
mut.clonality.data <- cbind(mut.clonality.data, glioma.molecular.subtype[mut.clonality.data$patient, 'explore.subtype', FALSE])

mut.clonality.data <- subset(mut.clonality.data, explore.subtype %in% 
 c('IDHcolmut.codel', 'IDHsubmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel'))
mut.clonality.data$explore.subtype <- factor(mut.clonality.data$explore.subtype, 
	levels = c('IDHcolmut.codel', 'IDHsubmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel'))


non.codel.ccf.pvalue <- wilcox.test(absolute.ccf~explore.subtype, subset(mut.clonality.data, explore.subtype %in% 
	c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value # 2.182084e-21

codel.ccf.pvalue <- wilcox.test(absolute.ccf~explore.subtype, subset(mut.clonality.data, explore.subtype %in% 
	c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value # 3.887847e-33


# random p value
# load('/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Section3/Results/Glioma/unbalancedAnalysis/nonCodelGliomaRandomSubtype.RData')
# load(file = '/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Section4/Results/MutationBurden/CodelRandomSubtype.RData')

# non.codel.random.ccf.pvalue <- sapply(seq(1000), function(times){
 
 # random.subtype <- glioma.random.subtype[[times]]
 # random.subtype <- cbind(com.data, random.subtype[com.data$patient, 'explore.subtype', FALSE])

 # ccf.pvalue <- wilcox.test(absolute.ccf~explore.subtype, random.subtype, alternative = 'two.sided')$p.value
 
 # return(ccf.pvalue)
# })

# codel.random.ccf.pvalue <- sapply(seq(1000), function(times){
 
 # random.subtype <- codel.random.subtype[[times]]
 # random.subtype <- cbind(com.data, random.subtype[com.data$patient, 'explore.subtype', FALSE])

 # ccf.pvalue <- wilcox.test(absolute.ccf~explore.subtype, random.subtype, alternative = 'two.sided')$p.value
 
 # return(ccf.pvalue)
# })

# sum(non.codel.random.ccf.pvalue < non.codel.ccf.pvalue)/1000 # 0.002
# sum(codel.random.ccf.pvalue < codel.ccf.pvalue)/1000 # 0


# Cumulative distribution of the somatic mutations by CCF
ccf.ecdf.plot <- function(ccfData, fileName, filePath){
 library('ggplot2')
	
 ccf.cumu <- ggplot(data = ccfData, aes(absolute.ccf)) + 
  stat_ecdf(aes(colour = explore.subtype), geom="line", size = 1) + 
  xlab('Cancer cell fraction') + ylab('Fraction of mutations') # + guides(colour = 'none')

 ggsave(filename = fileName, plot = ccf.cumu, path = filePath)
}
	

ccf.ecdf.plot(mut.clonality.data, 'CCFCumulativeDensity.pdf', filePath=paste(getwd(), '/result/4.Intratumoral_Heterogeneity_Analysis/MutationBurden/', sep=""))


