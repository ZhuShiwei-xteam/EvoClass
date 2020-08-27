
# determined the tumor subclonal composition
IdentifyMutCopyNumber <- function(x, sub.mat.mut, sub.mat.copy, gender){
  
 mut <- sub.mat.mut[x,,drop=FALSE]
 ww <- which(as.numeric(sub.mat.copy$Chr)==as.numeric(mut$Chr)
  &as.numeric(sub.mat.copy$Start)<=as.numeric(mut$Start_position)
  &as.numeric(sub.mat.copy$End)>=as.numeric(mut$Start_position))
 
 copy <- sub.mat.copy[ww,,drop=FALSE]
 mutation_id <- paste(mut$Patient,mut$Chr,mut$Start_position,mut$Reference,sep=":")
 ref_counts <- mut$Ref_freq
 var_counts <- mut$Variant_freq
 normal_cn <- 2
  
 # normal copy of sex chromosomes in male
 
 if(!is.na(gender)){
  
  if((gender=='male')&((mut$Chr==23)|(mut$Chr==24))){
    
	normal_cn <- 1
   }
  }
  
  patient <- mut$Patient
  Reference_Base <- mut$Reference
  Alternate_Base <- mut$Alternate
  Hugo_Symbol <- mut$Hugo_Symbol
  Variant_Classification <- mut$Variant_Classification
  HGVSp_Short <- mut$HGVSp_Short
  
 if(nrow(copy)!=1){
  minor_cn <- NA
  major_cn <- NA
  output <-  data.frame(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, 
   patient, Reference_Base, Alternate_Base, Variant_Classification, HGVSp_Short, Hugo_Symbol, stringsAsFactors=FALSE)
  
  return(output)
  }
  
 minor_cn  <- min(c(copy$nA,copy$nB))
 major_cn <- max(c(copy$nA,copy$nB))
  
 output  <- data.frame(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, patient, 
 Reference_Base, Alternate_Base, Variant_Classification, HGVSp_Short, Hugo_Symbol, stringsAsFactors=FALSE)

  return(output)
}


ConstructPycloneInputFile <- function(mutation.table, seg.table, data.type, TCGA.barcode, 
 sex.chr  = FALSE, gender = NA, ANALYSIS.DIR, min.var.prop = 0.05, min.alt.reads = 5, min.depth  = 30){
 # The input required is as follows
 # mutation.table # mutation tables as in ../input
 # seg.table # segmented matrices in list format as in ../input
 # data type # e.g. TCGA_BRCA
 # ANALYSIS.DIR #what directory do you want to create as a base directory, e.g. example/
 
 
 # create the root directory
 data.folder <- paste(ANALYSIS.DIR, data.type, sep="")
 
 if(!file.exists(data.folder)){
  if(!dir.create(data.folder, recursive = TRUE) ){
   
  stop("Unable to create root directory.\n")
  }
 }
  
 # create a folder for the specific patient
 patient.folder <- paste(data.folder, "/", TCGA.barcode, sep="")
  
 if( !file.exists(patient.folder)){
  if( !dir.create(patient.folder, recursive = TRUE) ){
   
   stop("Unable to create tmp directory.\n")
  }
 }
  
  
 # load somatic mutations ####

 req.mut.colnames <- c('Patient', 'Chr', 'Start_position', 'End_position', 'Reference', 'Alternate', 
 'Variant_freq', 'Ref_freq', 'Hugo_Symbol', 'Variant_Classification', 'HGVSp_Short')
  
 if(length(req.mut.colnames[!which(req.mut.colnames%in%colnames(mutation.table))])!=0){
  stop(paste('Mutation table not in correct format\nThe following columns are required:\n',
  PasteVector(req.mut.colnames, sep="\n"),sep=""))
 }
 
 
 # load somatic copy numbers ####
 seg.mat.copy <- seg.table
 req.seg.colnames <- c('SampleID', 'Chr', 'Start', 'End', 'nProbes', 'cn', 'nA', 'nB', 'Ploidy', 'Aberrant Cell Fraction')
  
 if(length(req.seg.colnames[!which(req.seg.colnames%in%colnames(seg.mat.copy))])!=0){
  stop(paste('Segment copy number not in correct format\nThe following columns are required:\n',  
  PasteVector(req.seg.colnames, sep="\n"),sep=""))
 }

						
 # combine copy number and mutation data #####
 # first of all, select only samples where both available
 barcodes <- intersect(unique(seg.mat.copy$SampleID),unique(mutation.table[,1]))
  cat(paste("\n\nThere are ", length(barcodes), " patients with both copy number and mutation data. \nOnly these will be used "))
  
 if(!TCGA.barcode%in%barcodes){
  stop("TCGA barcode not found in either mutation table or seg.mat.copy")
 }
  
 # select patient specific data
 sub.mat.copy <- seg.mat.copy[seg.mat.copy$SampleID==TCGA.barcode, , drop=FALSE]
 sub.mat.mut <- mutation.table[mutation.table[,1]==TCGA.barcode, , drop=FALSE]
  
 # remove sex chromosomes
 if((!sex.chr)|(is.na(gender))){
  sub.mat.copy <- sub.mat.copy[sub.mat.copy$Chr%in%c(1:22), ]
  sub.mat.copy$Chr <- as.numeric(sub.mat.copy$Chr)
  sub.mat.mut <- sub.mat.mut[sub.mat.mut$Chr%in%c(1:22), ]
  sub.mat.mut$Chr <- as.numeric(sub.mat.mut$Chr)
  }
  
 # combine mutation and copy number
 mut.table  <- data.frame(t(sapply(1:nrow(sub.mat.mut), IdentifyMutCopyNumber, sub.mat.mut, sub.mat.copy,gender)), stringsAsFactors=FALSE)
 mut.table  <- mut.table[!is.na(mut.table$minor_cn), ]
 mut.table <- mut.table[!is.na(mut.table$ref_counts), ]
 mut.table <- mut.table[!duplicated(mut.table$mutation_id), ]
 
 
 # min.alt.reads.cut.off
 mut.table <- mut.table[as.numeric(mut.table$var_counts)>=min.alt.reads,]
 
 # min.cov.at.site
 mut.table <- mut.table[as.numeric(mut.table$var_counts)+as.numeric(mut.table$ref_counts)>=min.depth,]
 
 # min variant proportion
 mut.table <- mut.table[(as.numeric(mut.table$var_counts)/(as.numeric(mut.table$var_counts)+as.numeric(mut.table$ref_counts)))>=min.var.prop,]
 

 mut.table <- apply(mut.table, 2, as.character)
 
 # write 
 output.file   <- paste(patient.folder, "/",TCGA.barcode, ".tsv", sep="")  
 write.table(mut.table, sep = "\t", quote=FALSE, col.names=TRUE, row.names=FALSE, file = output.file)
}


PycloneInputData <- function(cancer.type, panCan.sex.data, panCan.mut.data, panCan.seg.data, output.file.path){
 
 sex.data <- panCan.sex.data[cancer.type][[1]]
 
 lapply(seq(nrow(sex.data)), function(index){
  sample.sex.data <- sex.data[index, , drop=FALSE]
  
  try(ConstructPycloneInputFile(mutation.table  = panCan.mut.data, seg.table = panCan.seg.data, 
   data.type = paste0('TCGA_', cancer.type), TCGA.barcode = sample.sex.data$sample, sex.chr = FALSE, gender = sample.sex.data$gender, 
   ANALYSIS.DIR = output.file.path, min.var.prop = 0 , min.alt.reads = 0, min.depth = 0), silent = TRUE)
  })

}


# Run PyClone
RunPyClone <- function(input.file.path, tumor.purity.data, s.index, e.index){
 
 setwd(input.file.path)	
 file.names <- dir()	

 for(sample.name in file.names[s.index:e.index]){
 
  setwd(paste0(input.file.path, sample.name))
  file.name <- dir()
 
  if(length(file.name)==0){
   paste('Did not find mutation MAF file:\n', file.name, sep='')
   
   next
  }
 
  if(!(sample.name %in% rownames(tumor.purity.data))){
   paste('Lack of tumor purity:\n', file.name, sep='')
   
   next
  }
 
 command <- paste(c('PyClone', 'run_analysis_pipeline', '--in_files', paste0(input.file.path, sample.name, '/', sample.name, '.tsv'), 
 '--working_dir', paste0(input.file.path, sample.name), '--tumour_contents', 
 tumor.purity.data[sample.name, ]$'Aberrant Cell Fraction', '--density pyclone_binomial', 
 '--num_iters 20000', '--prior major_copy_number', '--burnin 5000', '--thin 10', '--seed 41'), collapse = ' ')
 
 system(command)
 }
 
}


PyCloneMutCCF <- function(file.path){
 setwd(file.path)
 patients <- dir()
 mut.ccf.data <- NULL
 
 for(patient in patients){
  
  setwd(paste0(file.path, patient))
  file.names <- dir()
  
  if('tables' %in% file.names){
   
   setwd(paste0(file.path, patient, '/', 'tables'))
   mut.ccf <- read.table(file = 'loci.tsv', sep = '\t', header =  TRUE, stringsAsFactors = FALSE)
  
  }else{
   
   next
  }
  
  mut.ccf.data <- rbind(mut.ccf.data, mut.ccf)
 }
 
 return(mut.ccf.data)
}


PyCloneClonalPopulationNumber <- function(file.path){
 setwd(file.path)
 patients <- dir()
 
 clonal.population.num <- NULL
	
 for(patient in patients){
 
  setwd(paste0(file.path, patient))
  file.names <- dir()
 
  if('tables' %in% file.names){
		
   mut.data <- read.csv(file = paste0(patient, '.tsv'), sep = '\t', header =  TRUE, stringsAsFactors = FALSE)
 
   setwd(paste0(file.path, patient, '/', 'tables'))
   clonal.population.data <- read.csv(file = 'loci.tsv', sep = '\t', header =  TRUE, stringsAsFactors = FALSE)
   clonal.population.data <- merge(mut.data, clonal.population.data, by = 'mutation_id')
   
   clonals <- split(clonal.population.data, factor(clonal.population.data$cluster_id))
   index <- sapply(clonals, function(clonal){
    s1 <- nrow(clonal) >= 2
	s2 <- nrow(subset(clonal, !(Variant_Classification %in% c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "Intron", "RNA", "Silent")))) >= 1
	
	return(all(c(s1, s2)))
   })
  
  population.num <- sum(index)
  names(population.num) <- patient
  clonal.population.num <- c(clonal.population.num, population.num)
 
  }else{
  
    next
  }
 
 }
 
 return(clonal.population.num)
}

