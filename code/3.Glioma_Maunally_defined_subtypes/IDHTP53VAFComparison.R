
setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')
library('ggpubr')


load('data/glioma.mut.ccf.RData')
load('data/GliomaMolecularSubtype.RData')


com.data$patient <- substr(com.data$patient, 1, 12)
com.data <- cbind(com.data, glioma.molecular.subtype[com.data$patient, ])

idh.tp53.vaf <- subset(com.data, final.subtype %in% 'IDHsubmut.non.codel' & Hugo_Symbol %in% c('IDH1', 'TP53'))

# wilcox.test(obs.VAF~Hugo_Symbol, idh.tp53.vaf, alternative = 'two.sided')$p.value
# 2.369457e-08

# wilcox.test(absolute.ccf~Hugo_Symbol, idh.tp53.vaf, alternative = 'two.sided')$p.value
# 3.518791e-09

idh.tp53.vaf$LOH_Status <- ifelse(idh.tp53.vaf$minor_cn==0 & idh.tp53.vaf$major_cn==2, 'withLOH', 'withoutLOH')
idh.tp53.vaf$Hugo_Symbol <- paste0(idh.tp53.vaf$Hugo_Symbol, idh.tp53.vaf$LOH_Status)
idh.tp53.vaf$Hugo_Symbol <- factor(idh.tp53.vaf$Hugo_Symbol, levels = c('IDH1withoutLOH', 'TP53withoutLOH', 'TP53withLOH'))

# wilcox.test(obs.VAF~Hugo_Symbol, subset(idh.tp53.vaf, Hugo_Symbol %in% 
 # c('IDH1withoutLOH', 'TP53withoutLOH')), alternative = 'two.sided')$p.value
# 1.629085e-05



vaf.boxplot <- ggboxplot(idh.tp53.vaf, "Hugo_Symbol", "obs.VAF", add = "jitter", color = "Hugo_Symbol",
 ylab = "Variant allele frequency") + stat_compare_means(method = "wilcox.test", 
 comparisons = list(c('IDH1withoutLOH', 'TP53withoutLOH'), c('TP53withoutLOH', 'TP53withLOH'),
 c('IDH1withoutLOH', 'TP53withLOH'))) 

ggsave(vaf.boxplot, file = 'result/3.Glioma_Maunally_defined_subtypes/IDH1TP53VAF.pdf')

ccf.boxplot <- ggboxplot(idh.tp53.vaf, "Hugo_Symbol", "absolute.ccf", add = "jitter", color = "Hugo_Symbol",
 ylab = "Cancer cell fraction") + stat_compare_means(method = "wilcox.test", 
 comparisons = list(c('IDH1withoutLOH', 'TP53withoutLOH'), c('TP53withoutLOH', 'TP53withLOH'),
 c('IDH1withoutLOH', 'TP53withLOH'))) 

ggsave(ccf.boxplot, file = 'result/3.Glioma_Maunally_defined_subtypes/IDH1TP53CCF.pdf')

