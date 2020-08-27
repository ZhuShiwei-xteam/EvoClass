setwd("/pub5/xiaoyun/Jobs/J22/temp/Github/data")

source("code/function/ccf_wrapper_function.R")

# load corresponding data
load("centromere.RData")
load("PureGliomaMutData.RData")
load("panCanAtlas.abs.seg.RData")

# calculate ccf

sapply(unique(pure.glioma.mut.data$Patient), function(sin.sam){
  try(clonality.estimation(mutation.table = pure.glioma.mut.data
                           ,seg.table          = panCanAtlas.abs.seg
                           ,data.type          = "TCGA_glioma"
                           ,TCGA.barcode       = sin.sam
                           ,ANALYSIS.DIR       = getwd()
                           ,sub.clonal.cut.off = 1), silent = TRUE)
})

# merge ccf files for each sample
ccf.merge <- function(cancer.type, input.path, out.path){

  file.path <- paste0(input.path, paste0(paste0('TCGA_', cancer.type), '/'))
  setwd(file.path)
  file.names <- dir()
  com.data <- c()
  for (file.name in file.names) {
    print(file.name)
    setwd(paste0(file.path, file.name))
    f <- dir()[grep('.tsv', dir())]
    if (length(f) < 1) {
      next
    }
    ccf.data <- read.table(file = f, header = TRUE, stringsAsFactors = F, sep = '\t', quote = '')
    if(ncol(ccf.data)==1){
		ccf.data <- t(ccf.data)
		colnames(ccf.data) <- colnames(com.data) 
	}	
	com.data <- rbind(com.data, ccf.data)
  }
  # save data
  setwd(out.path)
  save(com.data, file = paste0(cancer.type,'.mut.ccf.RData'))
}

ccf.merge("glioma", getwd(), getwd())



