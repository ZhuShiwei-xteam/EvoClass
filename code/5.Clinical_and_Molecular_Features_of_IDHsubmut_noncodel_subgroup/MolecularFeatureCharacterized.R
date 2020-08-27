setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')
library('ggpubr')

# significantly amplified and deleted events
peak.region.data <- read.table(file = 'data/OriginalData/all_lesions.conf_99.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
peak.region.data <- peak.region.data[1:(nrow(peak.region.data)/2), ]

amp.peak.data <- peak.region.data[grep('Amplification Peak', peak.region.data$Unique.Name), ]
rownames(amp.peak.data) <- trimws(amp.peak.data$Descriptor, which = "right", whitespace = "[ \t\r\n]")
amp.peak.data <- amp.peak.data[, -c(1:9, ncol(amp.peak.data))]
colnames(amp.peak.data) <- gsub('\\.', '-', substr(colnames(amp.peak.data), 1, 12))
amp.peak.data[amp.peak.data == 2] <- 1


del.peak.data <- peak.region.data[grep('Deletion Peak', peak.region.data$Unique.Name), ]
rownames(del.peak.data) <- trimws(del.peak.data$Descriptor, which = "right", whitespace = "[ \t\r\n]")
del.peak.data <- del.peak.data[, -c(1:9, ncol(del.peak.data))]
colnames(del.peak.data) <- gsub('\\.', '-', substr(colnames(del.peak.data), 1, 12))


# significantly associated with IDHsubmut-non-codel subtype
load('data/GliomaMolecularSubtype.RData')
subtype.amp.peak.data <- merge(t(amp.peak.data), glioma.molecular.subtype, by  = 'row.names')
subtype.del.peak.data <- merge(t(del.peak.data), glioma.molecular.subtype, by  = 'row.names')

amp.peak.sig <- sapply(rownames(amp.peak.data), function(each.peak){
	
 pvalue <- try(fisher.test(table(subset(subtype.amp.peak.data, explore.subtype %in% 
  c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'), c(each.peak, 'explore.subtype'))))$p.value)
 
 return(pvalue)
})


del.peak.sig <- sapply(rownames(del.peak.data), function(each.peak){
	
 pvalue <- try(fisher.test(table(subset(subtype.del.peak.data, explore.subtype %in% 
  c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'), c(each.peak, 'explore.subtype'))))$p.value)
 
 return(pvalue)
})


# Barplot shows the frequencies in the indicated molecular subtypes of significant amplified peaks
amp.peak.sta <- lapply(rownames(amp.peak.data), function(each.peak){

 peak.sta <- table(subtype.amp.peak.data[, c(each.peak, 'final.subtype')])
 peak.feq <- round(prop.table(peak.sta, 2), 2)
 
 peak.sta <- data.frame(subtype = colnames(peak.feq), peak = rep(each.peak, ncol(peak.feq)), mut.feq = peak.feq['1', ])
 
 if(max(peak.sta$mut.feq) > 0.3){
  peak.sta$facet <- '>0.3'
  
 }else if(max(peak.sta$mut.feq)<0.3 & max(peak.sta$mut.feq) > 0.15){
  peak.sta$facet <- '0.15<<0.3'
  
 }else{
  
  peak.sta$facet <- '<0.15'
 }
 
 return(peak.sta)
}) 

amp.peak.sta <- do.call(rbind, amp.peak.sta)


amp.peak.sta$subtype <- factor(amp.peak.sta$subtype, levels = 
 c("IDHmut.codel", "IDHcolmut.non.codel", "IDHsubmut.non.codel", 'IDHwt'))

cutoff <- c('>0.3', '0.15<<0.3', '<0.15')

pdf('result/5.Clinical_and_Molecular_Features_of_IDHsubmut_noncodel_subgroup/SubtypeAmpPeakFreq.pdf')
for(index in seq(length(cutoff))){
 
 subtype.amp.freq <- ggbarplot(subset(amp.peak.sta, facet == cutoff[index]), 'peak', "mut.feq", color = "subtype", 
  fill = "subtype", palette = c("#DB7F1E", "#E64457", "#3E84BC","#2D6978"), xlab = 'Amplification peak', ylab = 'Frequency of cnv', 
  position = position_dodge2(width = 0.4, padding = 0))
 
 plot(subtype.amp.freq)
}
dev.off()


# Barplot shows the frequencies in the indicated molecular subtypes of significant deleted peaks
del.peak.sta <- lapply(rownames(del.peak.data), function(each.peak){

 peak.sta <- table(subtype.del.peak.data[, c(each.peak, 'final.subtype')])
 peak.feq <- round(prop.table(peak.sta, 2), 2)
 
 peak.sta <- data.frame(subtype = colnames(peak.feq), peak = rep(each.peak, ncol(peak.feq)), mut.feq = peak.feq['1', ])
 
 if(max(peak.sta$mut.feq) > 0.6){
  peak.sta$facet <- '>0.6'
  
 }else if(max(peak.sta$mut.feq)<0.6 & max(peak.sta$mut.feq) > 0.28){
  peak.sta$facet <- '0.28<<0.6'
  
 }else if(max(peak.sta$mut.feq)<0.28 & max(peak.sta$mut.feq) > 0.16){
 
  peak.sta$facet <- '0.16<<0.28'
 }else{
  
  peak.sta$facet <- '<0.16'
 }
 
 return(peak.sta)
}) 

del.peak.sta <- do.call(rbind, del.peak.sta)

library('ggpubr')
del.peak.sta$subtype <- factor(del.peak.sta$subtype, levels = 
 c("IDHmut.codel", "IDHcolmut.non.codel", "IDHsubmut.non.codel", 'IDHwt'))

cutoff <- c('>0.6', '0.28<<0.6', '0.16<<0.28', '<0.16')

pdf('result/5.Clinical_and_Molecular_Features_of_IDHsubmut_noncodel_subgroup/SubtypeDelPeakFreq.pdf')
for(index in seq(length(cutoff))){
 
 subtype.del.freq <- ggbarplot(subset(del.peak.sta, facet == cutoff[index]), 'peak', "mut.feq", color = "subtype", 
  fill = "subtype", palette = c("#DB7F1E", "#E64457", "#3E84BC","#2D6978"), xlab = 'Deletion peak', ylab = 'Frequency of cnv', 
  position = position_dodge2(width = 0.4, padding = 0))
 
 plot(subtype.del.freq)
}
dev.off()



# molecular feature
load('data/gliomaClinicalData.RData')
rownames(pan.glio.mer.cli.data) <- pan.glio.mer.cli.data$bcr_patient_barcode

load(file = '/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Resources/GliomaMolecularSubtype.RData')
pan.glio.mer.cli.data <- merge(pan.glio.mer.cli.data, glioma.molecular.subtype, by = 'row.names')

pan.glio.mer.cli.data$final.subtype <- factor(pan.glio.mer.cli.data$final.subtype, 
 levels = c('IDHmut.codel', 'IDHcolmut.non.codel', 'IDHsubmut.non.codel', 'IDHwt'))

molecular.fea <- c('MGMT.promoter.status', 'ATRX.status', 'TERT.promoter.status', 'TERT.expression.status', 
 'Telomere.Maintenance','Chr.19.20.co.gain', 'Chr.7.gain.Chr.10.loss')
 
molecular.sig <- sapply(molecular.fea, function(feature){

 molecular.sta <- table(pan.glio.mer.cli.data[, c(feature, 'final.subtype')])
 
 across.subtype.pvalue <- chisq.test(molecular.sta)$p.value
 within.subtype.pvalue <- fisher.test(molecular.sta[, c('IDHcolmut.non.codel', 'IDHsubmut.non.codel')])$p.value

 return(c(across.subtype.pvalue, within.subtype.pvalue))
}) 
 
# stastic
molecular.sta <- sapply(molecular.fea, function(feature){

 mol.sta <- table(pan.glio.mer.cli.data[, c(feature, 'final.subtype')])
 mol.sta <- cbind(mol.sta, all.patient = rowSums(mol.sta))
 
 mol.feq <- round(prop.table(mol.sta, 2), 2)
 
 mol.sta.feq <- sapply(rownames(mol.sta), function(index){
  sta.feq <- paste0(mol.sta[index, ], '(', mol.feq[index, ], ')')
  
  return(sta.feq)
 })
 rownames(mol.sta.feq) <- colnames(mol.sta)

 return(mol.sta.feq)
}) 



# recurrent mutated genes
load('data/PureGliomaMutData.RData')
load(file = 'data/GliomaMolecularSubtype.RData')
load(file = 'data/pan.glioma.cancer.gene.RData')


ConstructedMutMatrix <- function(mutation.data, gene.set){

 variant.classification <- c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "Intron", "RNA", "Silent")	
 mutation.data <- mutation.data[!(mutation.data$Variant_Classification%in%variant.classification), ]
 
 patients <- unique(mutation.data$Patient)
 
 mut.matrix <- sapply(gene.set, function(gene){

  mut.data <- mutation.data[mutation.data$Hugo_Symbol == gene, , drop = FALSE]
  mut.patients <- unique(mut.data$Patient)
  mut.patient.lable <- ifelse(patients%in%mut.patients, 1, 0)
  
  return(mut.patient.lable)
 })
 
 rownames(mut.matrix) <- patients
 colnames(mut.matrix) <- gene.set
 
 return(mut.matrix)
}

cancer.gene.mut.matrix <- ConstructedMutMatrix(pure.glioma.mut.data, pan.glioma.cancer.gene)
rownames(cancer.gene.mut.matrix) <- substr(rownames(cancer.gene.mut.matrix), 1, 12)

cancer.gene.mut.matrix <- merge(cancer.gene.mut.matrix, glioma.molecular.subtype, by  = 'row.names')

driver.gene.sig <- sapply(pan.glioma.cancer.gene, function(driver.gene){
 
 pvalue <- NA
 index <- sum(subset(cancer.gene.mut.matrix, final.subtype %in% 
  c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'), driver.gene)) > 3
 
 if(index){
  
  pvalue <- try(fisher.test(table(subset(cancer.gene.mut.matrix, final.subtype %in% 
   c('IDHcolmut.non.codel', 'IDHsubmut.non.codel'), c(driver.gene, 'final.subtype'))))$p.value)
 }

 return(pvalue)
})


# Percentages for the indicated frequently mutated gene were shown in the bar graphs
driver.gene.sta <- lapply(pan.glioma.cancer.gene, function(driver.gene){

 driver.sta <- table(cancer.gene.mut.matrix [, c(driver.gene, 'final.subtype')])
 driver.feq <- round(prop.table(driver.sta, 2), 2)
 
 driver.sta <- data.frame(subtype = colnames(driver.feq), gene = rep(driver.gene, ncol(driver.feq)), mut.feq = driver.feq['1', ])
 
 if(max(driver.sta$mut.feq) > 0.2){
  driver.sta$facet <- '>0.2'
  
 }else if(max(driver.sta$mut.feq)<0.2 & max(driver.sta$mut.feq) > 0.07){
  driver.sta$facet <- '0.07<<0.2'
  
 }else{
  
  driver.sta$facet <- '<0.07'
 }
 
 return(driver.sta)
}) 

driver.gene.sta <- do.call(rbind, driver.gene.sta)

library('ggpubr')
driver.gene.sta$subtype <- factor(driver.gene.sta$subtype, levels = 
 c("IDHmut.codel", "IDHcolmut.non.codel", "IDHsubmut.non.codel", 'IDHwt'))

cutoff <- c('>0.2', '0.07<<0.2', '<0.07')

pdf('result/5.Clinical_and_Molecular_Features_of_IDHsubmut_noncodel_subgroup/SubtypeDriverFreq.pdf')
for(index in seq(length(cutoff))){
 
 subtype.driver.freq <- ggbarplot(subset(driver.gene.sta, facet == cutoff[index]), 'gene', "mut.feq", color = "subtype", 
  fill = "subtype", palette = c("#DB7F1E", "#E64457", "#3E84BC","#2D6978"), xlab = FALSE, ylab = 'Frequency of mutation', 
  position = position_dodge2(width = 0.4, padding = 0))
 
 plot(subtype.driver.freq)
}
dev.off()

