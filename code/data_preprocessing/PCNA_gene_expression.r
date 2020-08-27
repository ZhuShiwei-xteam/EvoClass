setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')
library(clusterProfiler)
library(org.Hs.eg.db)

setwd(paste(getwd(), "/data/OriginalData/GDC_FPKM", sep=""))
# untar
utils::untar("gdc_download_20200722_023950.634550.tar.gz",  exdir = ".")
file_names <- list.files()
file_names <- file_names[1:(length(file_names)-3)]
# sample tag
samples.tag <- read.table(file = "/pub5/xiaoyun/Jobs/J22/EvoClass2.0/Revise/de_expression/FPKM/gdc_sample_sheet.2020-07-22.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# load PCNA genes
PCNA.signature <- read.table(file = "data/PCNA_Signature.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)	
PCNA.signature.ensembl <- bitr(PCNA.signature[,1], fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Hs.eg.db") 

# get PCNA gene expression
PCNA.FPKM <- lapply(file_names, function(per.sample){
	sample.ID <- samples.tag[which(samples.tag$File.ID == per.sample), "Sample.ID"]
	exp.file.name <- list.files(paste("./", per.sample, sep=""))[grep("*.gz", list.files(paste("./", per.sample, sep="")))]
	FPKM.data <- read.table(paste("./", per.sample, "/", exp.file.name, sep=""))
	
	FPKM.ENSEMBL <- substr(FPKM.data[,1], 1, 15)
	# PCNA gene index
	PCNA.index <- which(FPKM.ENSEMBL %in% PCNA.signature.ensembl[,2])
	PCNA.FPKM.data <- FPKM.data[PCNA.index, 2, drop = FALSE]
	colnames(PCNA.FPKM.data) <- sample.ID
	
	PCNA.ensembl <- FPKM.ENSEMBL[PCNA.index]
	PCNA.entrz <- bitr(PCNA.ensembl, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[,2]
	rownames(PCNA.FPKM.data) <- PCNA.entrz
	PCNA.FPKM.data
})
glioma.PCNA.FPKM <- do.call(cbind, PCNA.FPKM)

setwd('/pub5/xiaoyun/Jobs/J22/temp/Github')
save(glioma.PCNA.FPKM, file = "data/glioma.PCNA.FPKM.RData")
