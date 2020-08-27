setwd("/pub5/xiaoyun/Jobs/J22/temp/Github")


CliSamCountStat <- function(cli.features, cli.data){
  
  sam.count <- dim(cli.data)[1]
  cli.stat.result <- lapply(cli.features, function(sin.cli.fea){
	
	# Distinguish the variable as continuous or categorical
	cli.varibale <- cli.data[, sin.cli.fea]
    if(all(is.na(cli.varibale))){
      value <- NA
    }else{
      cli.varibale <- na.omit(cli.varibale)
	  cli.varibale <- as.character(cli.varibale)
      value <- all(!is.na(as.numeric(cli.varibale))) 
    }
	
	if(is.na(value)){
	cli.df <- as.matrix(data.frame(cli.feature = sin.cli.fea, cli.feature.stat = 'Unknown'))
	
	}else if(value == TRUE){
	# continuous, especially for age
		cli.varibale <- as.numeric(cli.varibale)
		cli.df <- as.matrix(data.frame(cli.feature = "Age", cli.feature.stat = paste(median(cli.varibale), '(',min(cli.varibale), '-', max(cli.varibale), ')', sep='')))
	
	}else if(value == FALSE){
	# categorical
	cli.table <- table(cli.data[, sin.cli.fea])
    cli.table.frac <- cli.table/sam.count
    unknown.count <- sam.count - sum(cli.table)
    unknown.frac <- 1-sum(cli.table.frac)
    cli.feature.count <- c(as.numeric(cli.table), unknown.count)
    cli.feature.frac <- c(round(cli.table.frac, 2)*100, round(unknown.frac, 2)*100)
    
    cli.fea.count.frac <- paste(cli.feature.count, '(', cli.feature.frac, ')', sep = '')
    
    cli.df <- as.matrix(data.frame(cli.feature = c(sin.cli.fea, names(cli.table), 'Unknown'), cli.feature.stat = c(sin.cli.fea, cli.fea.count.frac), stringsAsFactors = FALSE))

	}
	t(cli.df)
  })
  
  cli.stat.result
}


load("data/gliomaClinicalData.RData") # load clinical data, pan.glio.mer.cli.data
load("data/glioma_core_sam.Rdata")# load core sample labels, glioma_core_sam
	
cli.features <- c("age_at_initial_pathologic_diagnosis", "neoplasm_histologic_grade", "gender", "molecular_histological_type","Extent_of_surgery_resection","laterality", "family_history_of_cancer", "tumor_location", "first_presenting_symptom",
	"MGMT.promoter.status", "ATRX.status", "TERT.promoter.status", "TERT.expression.status", "Chr.7.gain.Chr.10.loss", "Telomere.Maintenance")
reserch.sample.cli.data <- pan.glio.mer.cli.data[which(pan.glio.mer.cli.data[,1] %in% substr(glioma_core_sam, 1, 12)), ]
# Statistical analysis of clinical molecular characteristics on all samples

all.statistic <- CliSamCountStat(cli.features, reserch.sample.cli.data)
# Statistical analysis of clinical molecular characteristics on each subgroups
subtype.name <- unique(reserch.sample.cli.data$IDH.codel.subtype)
subtype.statistic <- lapply(subtype.name, function(per.type){
	per.subtype.characters <- subset(reserch.sample.cli.data, IDH.codel.subtype == per.type)
	per.subtype.statistic <- CliSamCountStat(cli.features, per.subtype.characters)
	per.subtype.statistic
})
names(subtype.statistic) <- subtype.name
save(all.statistic, subtype.statistic, file = "result/Patient_Cohort_Characteristics/clinical.statistic.RData")

	