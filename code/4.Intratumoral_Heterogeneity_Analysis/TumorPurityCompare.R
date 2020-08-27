setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')
library('ggplot2')

# tumor purity calculated by Absolute
absolute.purity <- read.table(file = 'data/OriginalData/TCGA_mastercalls.abs_tables_JSedit.fixed.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
absolute.purity <- absolute.purity[substr(absolute.purity$array,14, 15) %in% '01', ]
rownames(absolute.purity) <- substr(absolute.purity$array, 1, 12)


load("data/GliomaMolecularSubtype.RData")
glioma.molecular.subtype <- merge(glioma.molecular.subtype, absolute.purity, by = 'row.names')

glioma.molecular.subtype$explore.subtype <- factor(glioma.molecular.subtype$explore.subtype, 
	levels = c('IDHcolmut.codel', 'IDHsubmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel', 'IDHwt'))


non.codel.purity.pvalue <- wilcox.test(purity~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
 c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value # p-value = 4.424657e-07
codel.purity.pvalue <- wilcox.test(purity~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
 c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value # p-value = 0.1542



PurityBoxplot <- function(purity.data, file.name, file.path){	
	
 purity.box <- ggplot(data = purity.data, aes(explore.subtype, purity)) + geom_boxplot(aes(fill = explore.subtype), outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + xlab('Molecular subtypes') + ylab('Purity') + guides(fill = 'none') + 
  theme(axis.text.x  = element_text(angle=30, vjust=0.5)) 
 	
 ggsave(filename = file.name, plot = purity.box, path = file.path)
}


PurityBoxplot(glioma.molecular.subtype, 'AbsolutePurityBoxplot.pdf', file.path='result/4.Intratumoral_Heterogeneity_Analysis/PurityCompare')


# tumor purity calculated by CHAT
load('data/TCGA_cancers_purity.RData') 
CHAT.purity <- rbind(TCGA_cancers_purity$lgg_purity, TCGA_cancers_purity$gbm_purity)
CHAT.purity <- CHAT.purity[substr(CHAT.purity$sampleid, 14, 15) == '01', ]

CHAT.purity <- CHAT.purity[-which(duplicated(substr(CHAT.purity$sampleid, 1, 12))), ]
rownames(CHAT.purity) <- substr(CHAT.purity$sampleid, 1, 12)


load("data/GliomaMolecularSubtype.RData")
glioma.molecular.subtype <- merge(glioma.molecular.subtype, CHAT.purity, by = 'row.names')

glioma.molecular.subtype$explore.subtype <- factor(glioma.molecular.subtype$explore.subtype, 
	levels = c('IDHcolmut.codel', 'IDHsubmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel', 'IDHwt'))


non.codel.purity.pvalue <- wilcox.test(purity~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
 c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')), alternative = 'two.sided')$p.value # p-value = 0.0007102492
codel.purity.pvalue <- wilcox.test(purity~explore.subtype, subset(glioma.molecular.subtype, explore.subtype %in% 
 c('IDHcolmut.codel', 'IDHsubmut.codel')), alternative = 'two.sided')$p.value # p-value = 0.4480102


PurityBoxplot(glioma.molecular.subtype, 'CHATPurityBoxplot.pdf', file.path='result/4.Intratumoral_Heterogeneity_Analysis/PurityCompare')

