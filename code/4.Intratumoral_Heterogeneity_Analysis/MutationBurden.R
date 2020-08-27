setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")



# mutation berden
load('data/PureGliomaMutData.RData')

glioma.mut.burden <- data.frame(mut.burden = tapply(pure.glioma.mut.data$mutation.id, factor(pure.glioma.mut.data$Patient), length))
rownames(glioma.mut.burden) <- substr(rownames(glioma.mut.burden), 1, 12)

# clonal and subclonal mutation burden
load("data/glioma.mut.ccf.RData")
mut.clonality.data <- com.data

glioma.clonality.burden <- split(mut.clonality.data, factor(mut.clonality.data$patient))
glioma.clonality.burden <- lapply(glioma.clonality.burden, function(clonality.data){
 
 clonal.num <- sum(clonality.data$CI95.timing == 'Clonal')
 subclonal.num <- sum(clonality.data$CI95.timing == 'Subclonal')

 return(c(clonal.num, subclonal.num))
})

glioma.clonality.burden <- as.data.frame(do.call(rbind, glioma.clonality.burden))
colnames(glioma.clonality.burden) <- c('cmut.burden', 'smut.burden')
rownames(glioma.clonality.burden) <- substr(rownames(glioma.clonality.burden), 1, 12)

# molecular subtype
load("data/GliomaMolecularSubtype.RData")

glioma.molecular.subtype <- cbind(glioma.molecular.subtype, glioma.clonality.burden[rownames(glioma.molecular.subtype), ], 
 glioma.mut.burden[rownames(glioma.molecular.subtype), ,drop = FALSE])


non.codel.mut.pvalue <- wilcox.test(mut.burden~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
	c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value # p-value = 0.8346
codel.mut.pvalue <- wilcox.test(mut.burden~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
	c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value # p-value = 0.8096

non.codel.cmut.pvalue <- wilcox.test(cmut.burden~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
	c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value # p-value = 0.05196
codel.cmut.pvalue <- wilcox.test(cmut.burden~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
	c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value # p-value = 0.004786

non.codel.smut.pvalue <- wilcox.test(smut.burden~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
	c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value # p-value = 0.004643
codel.smut.pvalue <- wilcox.test(smut.burden~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
	c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value # p-value = 0.01441


# non-codel random p value
# load('/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Section3/Results/Glioma/unbalancedAnalysis/nonCodelGliomaRandomSubtype.RData')

# random.mutBurden.pvalue <- lapply(seq(1000), function(times){
 
 # random.subtype <- glioma.random.subtype[[times]]
 # random.subtype <- merge(random.subtype, glioma.molecular.subtype[, c('cmut.burden', 'smut.burden', 'mut.burden')], by = 'row.names')

 # mut.pvalue <- wilcox.test(mut.burden~explore.subtype, subset(random.subtype, explore.subtype %in% 
	# c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value

 # cmut.pvalue <- wilcox.test(cmut.burden~explore.subtype, subset(random.subtype, explore.subtype %in% 
	# c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value
 
 # smut.pvalue <- wilcox.test(smut.burden~explore.subtype, subset(random.subtype, explore.subtype %in% 
	# c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value
 
 
 # return(c(mut.pvalue, cmut.pvalue, smut.pvalue))
# })

# random.mutBurden.pvalue <- do.call(rbind, random.mutBurden.pvalue)

# sum(random.mutBurden.pvalue[, 1] < non.codel.mut.pvalue)/nrow(random.mutBurden.pvalue) # 0.859
# sum(random.mutBurden.pvalue[, 2] < non.codel.cmut.pvalue)/nrow(random.mutBurden.pvalue) # 0.051
# sum(random.mutBurden.pvalue[, 3] < non.codel.smut.pvalue)/nrow(random.mutBurden.pvalue) # 0.009


# # codel random p value
# codel.subtype <- subset(glioma.molecular.subtype, 
 # explore.subtype %in% c('IDHcolmut.codel', 'IDHsubmut.codel'), explore.subtype)

# set.seed(1024)
# codel.random.subtype <- lapply(seq(1000), function(times){
 # random.subtype <- codel.subtype
 # random.subtype$explore.subtype <- sample(random.subtype$explore.subtype, length(random.subtype$explore.subtype), replace=FALSE) 
 
 # return(random.subtype)
# })
# # save(codel.random.subtype, file = '/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Section4/Results/MutationBurden/CodelRandomSubtype.RData')

# random.mutBurden.pvalue <- lapply(seq(1000), function(times){
 
 # random.subtype <- codel.random.subtype[[times]]
 # random.subtype <- merge(random.subtype, glioma.molecular.subtype[, c('cmut.burden', 'smut.burden', 'mut.burden')], by = 'row.names')

 # mut.pvalue <- wilcox.test(mut.burden~explore.subtype, subset(random.subtype, explore.subtype %in% 
	# c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value

 # cmut.pvalue <- wilcox.test(cmut.burden~explore.subtype, subset(random.subtype, explore.subtype %in% 
	# c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value
 
 # smut.pvalue <- wilcox.test(smut.burden~explore.subtype, subset(random.subtype, explore.subtype %in% 
	# c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value
 
 
 # return(c(mut.pvalue, cmut.pvalue, smut.pvalue))
# })

# random.mutBurden.pvalue <- do.call(rbind, random.mutBurden.pvalue)

# sum(random.mutBurden.pvalue[, 1] < codel.mut.pvalue)/nrow(random.mutBurden.pvalue) # 0.798
# sum(random.mutBurden.pvalue[, 2] < codel.cmut.pvalue)/nrow(random.mutBurden.pvalue) # 0.006
# sum(random.mutBurden.pvalue[, 3] < codel.smut.pvalue)/nrow(random.mutBurden.pvalue) # 0.008


# boxplot-mutation berden
MutBurdenBoxplot <- function(mutBurden, fileName, filePath){
 library('ggplot2')
	
 mut.box <- ggplot(data = mutBurden, aes(explore.subtype, mut.burden)) +
  geom_boxplot(aes(fill = explore.subtype), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_y_continuous(breaks=c(0, 10, 30, 50, 75, 100), limits=c(0, 100)) + 
  xlab('Molecular subtypes') + ylab('Overall mutation burden') + guides(fill = 'none') + 
  theme(axis.text.x  = element_text(angle=30, vjust=0.5))
		
  ggsave(filename = fileName, plot = mut.box, path = filePath)
}


glioma.molecular.subtype <- subset(glioma.molecular.subtype, explore.subtype != 'IDHwt')
glioma.molecular.subtype$explore.subtype <- factor(glioma.molecular.subtype$explore.subtype, 
	levels = c('IDHcolmut.codel', 'IDHsubmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel'))

MutBurdenBoxplot(glioma.molecular.subtype, 'MutBurdenBoxplot.pdf', 
 filePath=paste(getwd(),'/result/4.Intratumoral_Heterogeneity_Analysis/MutationBurden', sep=""))


# boxplot-clonal mutation berden
ClonalMutBurdenBoxplot <- function(mutBurden, fileName, filePath){
 library('ggplot2')
	
 mut.box <- ggplot(data = mutBurden, aes(explore.subtype, cmut.burden)) +
  geom_boxplot(aes(fill = explore.subtype), outlier.shape = NA) + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_y_continuous(breaks=c(0, 10, 30, 50, 70), limits=c(0, 70)) + 
  # scale_y_continuous(breaks=c(0, 50, 100, 200, 300), limits=c(0, 300)) + 
  xlab('Molecular subtypes') + ylab('Clonal mutation burden') + guides(fill = 'none') + 
  theme(axis.text.x  = element_text(angle=30, vjust=0.5))
		
 ggsave(filename = fileName, plot = mut.box, path = filePath)
}


ClonalMutBurdenBoxplot(glioma.molecular.subtype, 'ClonalMutBurdenBoxplot1.pdf', 
 filePath=paste(getwd(),'/result/4.Intratumoral_Heterogeneity_Analysis/MutationBurden', sep=""))


# boxplot-subclonal mutation berden
SubclonalMutBurdenBoxplot <- function(mutBurden, fileName, filePath){
 library('ggplot2')
	
 mut.box <- ggplot(data = mutBurden, aes(explore.subtype, smut.burden)) +
  geom_boxplot(aes(fill = explore.subtype), outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_y_continuous(breaks=c(0, 10, 30, 50, 70),   limits=c(0, 70))+
  # scale_y_continuous(breaks=c(0, 50, 100, 200, 300), limits=c(0, 300)) + 
  xlab('Molecular subtypes') + ylab('Subclonal mutation burden') + guides(fill = 'none') + 
  theme(axis.text.x  = element_text(angle=30, vjust=0.5))
		
 ggsave(filename = fileName, plot = mut.box, path = filePath)
}

SubclonalMutBurdenBoxplot(glioma.molecular.subtype, 'SubclonalMutBurdenBoxplot1.pdf', 
  filePath=paste(getwd(),'/result/4.Intratumoral_Heterogeneity_Analysis/MutationBurden', sep=""))



