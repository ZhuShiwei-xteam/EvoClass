
# constructing the mutation matrix
ConstructMutMatrix <- function(mut.data, gene.set){
 
 mut.data <- subset(mut.data, Hugo_Symbol %in% gene.set)
 
 sample.set <- unique(mut.data$patient)
 gene.set <- unique(mut.data$Hugo_Symbol)
 
 # constructing the matrix of clonal mutation
 clonal.mut.data <- subset(mut.data, CI95.timing == 'Clonal')
 
 clonal.mut.matrix <- sapply(1:length(gene.set), function(index){
  
  gene <- gene.set[index]
  clonal.gene.mut.data <- subset(clonal.mut.data, Hugo_Symbol == gene)
  
  gene.clonal.mut.sample <- unique(clonal.gene.mut.data$patient)
  return(ifelse(sample.set%in%gene.clonal.mut.sample, 1, 0))
  
 }) 
 
 rownames(clonal.mut.matrix) <- sample.set
 colnames(clonal.mut.matrix) <- paste0(gene.set, '(C)')
 
 # constructing the matrix of subclonal mutation 
 subclonal.mut.data <- subset(mut.data, CI95.timing == 'Subclonal')
 
 subclonal.mut.matrix <- sapply(1:length(gene.set), function(index){
  gene <- gene.set[index]
  subclonal.gene.mut.data <- subset(subclonal.mut.data, Hugo_Symbol == gene)
  
  gene.subclonal.mut.sample <- unique(subclonal.gene.mut.data$patient)
  return(ifelse(sample.set%in%gene.subclonal.mut.sample, 1, 0))
  
 })  
 
 rownames(subclonal.mut.matrix) <- sample.set
 colnames(subclonal.mut.matrix) <- paste0(gene.set, '(S)')
 
 # combine clonal and subclonal mutation matrix
 mut.matrix <- cbind(clonal.mut.matrix, subclonal.mut.matrix)
 return(mut.matrix)
}


# make oncosign input file
OncosignInputFile <- function(clonality.data, driver.gene.set, input.file.path){
 
 mut.matrix <- ConstructMutMatrix(clonality.data, driver.gene.set)
 
 rownames(mut.matrix) <- substr(rownames(mut.matrix), 1, 12)
 mut.matrix <- t(mut.matrix)
 mut.matrix <- cbind(Tumor = rownames(mut.matrix), mut.matrix)
 
 setwd(input.file.path)
 write.table(mut.matrix, file = 'mut.data.txt', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
 write.table(colnames(mut.matrix), file = 'sample.data.txt', sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
	
 return(NULL)
}

# extract OncoSign molecular subtype
SubtypeExtract <- function(oncosign.subtype, level){
	
 molecular.subtype <- oncosign.subtype[, paste0('V', 1+((1:level))), drop = FALSE]
	
 cancer.subtype <- apply(molecular.subtype, 1, function(subtype){
  subtype.index <- subtype[max(which(subtype != ''))]
  
  return(subtype.index)
  })
	
 cancer.subtype <- data.frame(patient = names(cancer.subtype), subtype = cancer.subtype, stringsAsFactors = FALSE)
 
 colnames(cancer.subtype) <- c('patient', paste0('level.', level, '.subtype'))
 return(cancer.subtype)
}

