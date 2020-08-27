setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")



# load clinical data
load('data/gliomaClinicalData.RData')
rownames(pan.glio.mer.cli.data) <- pan.glio.mer.cli.data$bcr_patient_barcode

# oncosign molecular subtype
load(file = 'result/3.Glioma_Maunally_defined_subtypes/OncosignMolecularClassification/OncosignMolecularSubtype.RData')
load(file = 'result/3.Glioma_Maunally_defined_subtypes/OncosignMolecularClassification/OncosignResults.RData')
rownames(molecular.subtype) <- molecular.subtype$patient


# The third main picture  
pan.glio.mer.cli.data <- cbind(pan.glio.mer.cli.data[rownames(molecular.subtype), ], 
 molecular.subtype[, c('level.1.subtype', 'level.2.subtype')])

pan.glio.mer.cli.data <- pan.glio.mer.cli.data[, c('level.2.subtype', 'level.1.subtype', 'IDH.codel.subtype', 
 'TERT.promoter.status', 'ATRX.status', 'MGMT.promoter.status', 'molecular_histological_type', 
 'Grade.rough', 'gender', 'age_at_initial_pathologic_diagnosis')]
 
colnames(pan.glio.mer.cli.data) <- c('Level 2 classes', 'Level 1 classes', 'Molecular Subtype', 'TERT', 
 'ATRX', 'MGMT', 'Histologic Class', 'Grade', 'Sex', 'Age')
 
pan.glio.mer.cli.data$Age <- ifelse(pan.glio.mer.cli.data$Age >= median(pan.glio.mer.cli.data$Age), '>=52 yr', '<52 yr')


molecular.subtype$level.3.subtype <- factor(molecular.subtype$level.3.subtype, levels = c('OSC.2.6', 'OSC.2.8.1', 'OSC.2.8.2', 
 'OSC.2.10.10', 'OSC.2.10.12', 'OSC.2.10.8', 'OSC.2.10.11', 'OSC.2.10.9', 'OSC.1', 'OSC.3.5.3', 'OSC.3.5.13', 'OSC.3.5.4', 
 'OSC.3.7.6', 'OSC.3.7.5', 'OSC.3.7.7', 'OSC.3.4.18', 'OSC.3.4.15', 'OSC.3.4.16', 'OSC.3.4.14', 'OSC.3.4.17', 
 'OSC.3.2', 'OSC.3.1', 'OSC.3.3', 'OSC.3.9'))


ShowSubtypeWithHeatmap <- function(mut.matrix, subtype.label, clinical.data, output.file.name){
 
 library("pheatmap")
 
 # commom samples
 common.samples <- Reduce(intersect, list(colnames(mut.matrix), rownames(subtype.label), rownames(clinical.data)))
 
 clinical.data <- clinical.data[common.samples, ]
 subtype.label <- subtype.label[common.samples, , FALSE]
 mut.matrix <- mut.matrix[, common.samples]
 
 clinical.data <- cbind(subtype.label, clinical.data)
 
 # mutually exclusive plot
 mut.matrix <- mut.matrix[order(abs(rowSums(mut.matrix)), decreasing=T), ]
 subtype.label <- split(subtype.label, factor(subtype.label[, 1]))
 
 mut.matrix <- lapply(subtype.label, function(each.subtype){
  flag <- paste(rep(letters,rep(10, length(letters))), 0:9, sep = "")	
  sample.index <- c()
  tmp.mut.matrix <- mut.matrix[, rownames(each.subtype)]
 
  for(i in 1:ncol(tmp.mut.matrix)){
   
   index <- which(tmp.mut.matrix[,i]!=0)
   sample.index <- c(sample.index, paste(flag[index], collapse = ""))
  }
  
  tmp.mut.matrix <- tmp.mut.matrix[, order(sample.index)]
  
  
 if(each.subtype[1, 1] %in% c('OSC.2.10.11')){
  tmp.mut.matrix <- rbind(tmp.mut.matrix, clinical.data[colnames(tmp.mut.matrix), 'Molecular Subtype'])
  tmp.mut.matrix <- tmp.mut.matrix[, order(factor(tmp.mut.matrix[nrow(tmp.mut.matrix), ], 
   levels = c('IDHmut-non-codel', 'IDHmut-codel', 'IDHwt')), decreasing = TRUE)]
  
  tmp.mut.matrix <- tmp.mut.matrix[1:(nrow(tmp.mut.matrix)-1), ]
 }
 
 if(each.subtype[1, 1] %in% c('OSC.3.5.13')){
  tmp.mut.matrix <- rbind(tmp.mut.matrix, clinical.data[colnames(tmp.mut.matrix), 'Molecular Subtype'])
  tmp.mut.matrix <- tmp.mut.matrix[, order(factor(tmp.mut.matrix[nrow(tmp.mut.matrix), ], 
   levels = c('IDHwt', 'IDHmut-codel', 'IDHmut-non-codel')), decreasing = TRUE)]
   
  tmp.mut.matrix <- tmp.mut.matrix[1:(nrow(tmp.mut.matrix)-1), ]
 }
 
  return(tmp.mut.matrix)
 })	
 
 mut.matrix <- Reduce(cbind, mut.matrix)
 
 # clonal mutations were encoded with 2
 clo.mut.matrix <- mut.matrix[grep('(C)', rownames(mut.matrix), fixed = TRUE), ]
 clo.mut.matrix[clo.mut.matrix == 1] <- 2
 mut.matrix[grep('(C)', rownames(mut.matrix), fixed = TRUE), ] <- clo.mut.matrix
 
 mut.matrix[] <- lapply(mut.matrix, type.convert)
 # pheatmap
 cols <- c('#FDFEFE', '#FFF176', '#EF5350')
 pdf(file = paste0(output.file.name, ".pdf"))
  pheatmap(mut.matrix, annotation_col = clinical.data, cluster_cols = FALSE, cluster_row = FALSE, show_colnames = FALSE, 
  color = cols, breaks = c(-0.5, 0.5, 1.5, 2.5), legend_breaks = c(0, 1, 2), legend_labels = c('W', 'S', 'C'))
 dev.off()
 
}

setwd(paste(getwd(), '/result/3.Glioma_Maunally_defined_subtypes/OncosignMolecularClassification/', sep=""))
ShowSubtypeWithHeatmap(oncosign.results$data, molecular.subtype[, 4, FALSE], pan.glio.mer.cli.data, 'Level3Subtype')

