
 # overall survival
OSSurPlot <- function(os.data, group.label, risk.table.index = TRUE){
 
 library("survival") 
 library("survminer")
 
 rownames(os.data) <- os.data$bcr_patient_barcode
 common.samples <- intersect(rownames(group.label), rownames(os.data))
 
 os.data <- os.data[common.samples, ,drop = FALSE] 
 group.label <- group.label[common.samples, , drop = FALSE]
 
 
 if(length(unique(group.label[, 1]))<=1){
  
  return(NULL)
 }else{
 
  os.data$subtype <- group.label[, 1]
 
  cols <- c("#E7B800", "#2E9FDF", "#D1751E","#572875","#D8908F","#E0B665",
  "#ABC47A","#BA2021","#C7B3CD","#ACC7D8","#1E893A","#136295","#A45223", "#F5B041", "#AF7AC5", "#FF00FF", "#800080") 
  cols <- cols[seq(unique(os.data$subtype))]
 
  fit <- survfit(Surv(OS.time, OS) ~ subtype, data = os.data)
  
  os.plot <- ggsurvplot(data = os.data, fit, pval = TRUE, conf.int = FALSE, risk.table = risk.table.index, 
  risk.table.col = "strata", linetype = "strata", surv.median.line = "hv", xscale = 365.25, break.time.by = 365.25, 
  ggtheme = theme_bw(), palette = cols, xlab = 'Time(years)', ylab = 'Overall Survival(%)')
  
  return(os.plot)
 }
}


# disease specific survival
DSSSurPlot <- function(dss.data, group.label, risk.table.index = TRUE){

 library("survival")
 library("survminer")
 
 rownames(dss.data) <- dss.data$bcr_patient_barcode
 common.samples <- intersect(rownames(group.label), rownames(dss.data))
 
 dss.data <- dss.data[common.samples, , drop = FALSE] 
 group.label <- group.label[common.samples, , drop = FALSE]
 
 
 if(length(unique(group.label[, 1]))<=1){
  
  return(NULL)
 }else{
  
  dss.data$subtype <- group.label[, 1]
  
  cols <- c("#E7B800", "#2E9FDF", "#D1751E","#572875","#D8908F","#E0B665",
  "#ABC47A","#BA2021","#C7B3CD","#ACC7D8","#1E893A","#136295","#A45223", "#F5B041", "#AF7AC5", "#FF00FF", "#800080")
  cols <- cols[seq(unique(dss.data$subtype))]
  
  fit <- survfit(Surv(DSS.time, DSS) ~ subtype, data = dss.data)
  
  dss.plot <- ggsurvplot(data = dss.data, fit, pval = TRUE, conf.int = FALSE, risk.table = risk.table.index, 
  risk.table.col = "strata", linetype = "strata", surv.median.line = "hv", xscale = 365.25, break.time.by = 365.25,
  ggtheme = theme_bw(), palette = cols, xlab = 'Time(years)', ylab = 'Disease Specific Survival(%)')
  
  return(dss.plot)
 }
}


# progression free interval
PFISurPlot <- function(pfi.data, group.label, risk.table.index = TRUE){

 library("survival") 
 library("survminer")
 
 rownames(pfi.data) <- pfi.data$bcr_patient_barcode
 common.samples <- intersect(rownames(group.label), rownames(pfi.data))
 
 pfi.data <- pfi.data[common.samples, , drop = FALSE] 
 group.label <- group.label[common.samples, , drop = FALSE]
 

 if(length(unique(group.label[, 1]))<=1){
  
  return(NULL)
 }else{
   
  pfi.data$subtype <- group.label[, 1]
  
  cols <- c("#E7B800", "#2E9FDF", "#D1751E","#572875","#D8908F","#E0B665",
  "#ABC47A","#BA2021","#C7B3CD","#ACC7D8","#1E893A","#136295","#A45223", "#F5B041", "#AF7AC5", "#FF00FF", "#800080")
  cols <- cols[seq(unique(pfi.data$subtype))]
 
  fit <- survfit(Surv(PFI.time, PFI) ~ subtype, data = pfi.data)
  
  pfi.plot <- ggsurvplot(data = pfi.data, fit, pval = TRUE, conf.int = FALSE, risk.table = risk.table.index, 
  risk.table.col = "strata", linetype = "strata", surv.median.line = "hv", xscale = 365.25, break.time.by = 365.25, 
  ggtheme = theme_bw(), palette = cols, xlab = 'Time(years)', ylab = 'Progression Free Interval(%)')
 
  return(pfi.plot)
 }
}


# disease free interval
DFISurPlot <- function(dfi.data, group.label, risk.table.index = TRUE){

 library("survival") 
 library("survminer")
 
 rownames(dfi.data) <- dfi.data$bcr_patient_barcode
 commonSamples <- intersect(rownames(group.label), rownames(dfi.data))
 
 dfi.data <- dfi.data[commonSamples, , drop = FALSE] 
 group.label <- group.label[commonSamples, , drop = FALSE]
 
 if(length(unique(group.label[, 1]))<=1){
  
  return(NULL)
 }else{
  
  dfi.data$subtype <- group.label[, 1]
  
  cols <- c("#E7B800", "#2E9FDF", "#D1751E","#572875","#D8908F","#E0B665",
  "#ABC47A","#BA2021","#C7B3CD","#ACC7D8","#1E893A","#136295","#A45223", "#F5B041", "#AF7AC5", "#FF00FF", "#800080")
  cols <- cols[seq(unique(dfi.data$subtype))]
  
  fit <- survfit(Surv(DFI.time, DFI) ~ subtype, data = dfi.data)
  
  dfi.plot <- ggsurvplot(data = dfi.data, fit, pval = TRUE, conf.int = FALSE, risk.table = risk.table.index, 
  risk.table.col = "strata", linetype = "strata", surv.median.line = "hv", xscale = 365.25, break.time.by = 365.25, 
  ggtheme = theme_bw(), palette = cols, xlab = 'Time(years)', ylab = 'Disease Free Interval(%)')
  
  return(dfi.plot)
 }
}


SystemSurvivalPlot <- function(OS.data, DSS.data, PFI.data, DFI.data, cancer.subtype, output.file.name, risk.table.index = TRUE){

 os.plot <- OSSurPlot(OS.data, cancer.subtype, risk.table.index)
 dss.plot <- DSSSurPlot(DSS.data, cancer.subtype, risk.table.index)
 pfi.plot <- PFISurPlot(PFI.data, cancer.subtype, risk.table.index)
 dfi.plot <- DFISurPlot(DFI.data, cancer.subtype, risk.table.index)
 
 sur.plots <- list()
 sur.plots[[1]] <- os.plot
 sur.plots[[2]] <- dss.plot
 sur.plots[[3]] <- pfi.plot
 sur.plots[[4]] <- dfi.plot
 
 # Arrange multiple ggsurvplots and print the output
 require("survminer")
 
 res <- arrange_ggsurvplots(sur.plots, print = TRUE, ncol = 2, nrow = 2)
 ggsave(paste0(c(output.file.name, ".pdf"), collapse = ""), res)
}


