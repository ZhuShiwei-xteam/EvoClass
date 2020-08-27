

# Truncation time and events were extracted from clinical characteristic data
ClinicalOutcome <- function(clinical.data, survival.type){
  
  if(survival.type == 'OS'){
    time <- as.numeric(clinical.data$OS.time)
    event <- as.numeric(clinical.data$OS)
  }
  if(survival.type == 'DSS'){
    time <- as.numeric(clinical.data$DSS.time)
    event <- as.numeric(clinical.data$DSS)
  }
  if(survival.type == 'PFI'){
    time <- as.numeric(clinical.data$PFI.time)
    event <- as.numeric(clinical.data$PFI)
  }
  if(survival.type == 'DFI'){
    time <- as.numeric(clinical.data$DFI.time)
    event <- as.numeric(clinical.data$DFI)
  }
  clinical.outcome <- cbind(time,event)
  rownames(clinical.outcome) <- rownames(clinical.data)
  clinical.outcome <- as.data.frame(clinical.outcome)

  return(clinical.outcome)
  
}

library("survival") 
library("survminer")
KMplot <- function(clinical.data, risk.type, colors){
   
  common.samples <- intersect(rownames(risk.type), rownames(clinical.data))
  clinical.data <- clinical.data[common.samples, ,drop = FALSE] 
  risk.type <- risk.type[common.samples, , drop = FALSE]
  
  # If the number of groups is not ≥2 groups, survival judgment cannot be made
  if(length(unique(risk.type[, 1]))<=1){
    return(NULL)
  }else{
    # Add risk types to the data.frame of the clinical data
    clinical.data <- cbind(clinical.data, subtype = risk.type[ , 1])
    # Draw a survival curve
    fit <- survfit(Surv(time, event) ~ subtype, data = clinical.data)
    names(fit$strata) <- unlist(lapply(strsplit(names(fit$strata),"="),function(x) x[[2]]))# Only the name of the subtype is displayed in the graphic presentation
    km.plot <- ggsurvplot(data = clinical.data, 
                          fit, 
                          pval = TRUE, 
                          conf.int = FALSE, 
                          risk.table = TRUE, 
                          risk.table.col = "strata",
                          linetype = "strata", 
                          surv.median.line = "hv", 
                          ggtheme = theme_bw(), 
                          palette = colors)
    return(km.plot)
  }
}




# Univariate Cox regression analysis
CoxUniv <- function(clinical.data, interData=FALSE, clinical.factor=NULL, cancer.subtype,  survival.type){
  
  # Determine if it is a variable within Clinical.data
  if(interData==F){
    # the prognostic factor is not included in Clinical. Data
    commen.sample <- intersect(rownames(clinical.data),rownames(cancer.subtype))
    clinical.data <- cbind(clinical.data[commen.sample, , drop=F],cancer.subtype[commen.sample, , drop=F])
    clinical.factor <- names(cancer.subtype)
  }   
  index <- which(colnames(clinical.data) %in% clinical.factor)
  if(length(index)!=0){

    factor.index <- which(!is.na(clinical.data[,clinical.factor]))
    clinical.data <- clinical.data[factor.index, , drop=F]
    clinical.outcome <- ClinicalOutcome(clinical.data, survival.type)
    # Univariate Cox regression analysis
    surv <- Surv(clinical.outcome$time, clinical.outcome$event == 1)
    colnames(clinical.data)[index] <- "factor.VS."
    coxph.result <- coxph(surv ~ factor.VS., clinical.data)
    
	# check variable names in coxph result
    if(class(clinical.data[,index])=="integer" | class(clinical.data[,index])=="numeric"){ # continuous variable
      names(coxph.result$coefficients) <- clinical.factor
    }else{
      # categorical variable
      VS.names <- names(coxph.result$coefficients)
      all.classes <- unique(clinical.data$factor.VS.)
      names(coxph.result$coefficients) <- AdjustName(VS.names,all.classes) 
    }
    return(coxph.result)
  }else{
    return(NULL)
  }
}



# remove NA before multivariate Cox regression
RmnaBeforeMulCox <- function(clinical.data, clinical.factors, cancer.subtype=NULL){
  
  # Determine if it is a variable within Clinical.data
  if(!is.null(cancer.subtype)){
    commen.sample <- intersect(rownames(clinical.data),rownames(cancer.subtype))
    clinical.data <- cbind(clinical.data[commen.sample, , drop=F],cancer.subtype[commen.sample, , drop=F])
    clinical.factors <- c(clinical.factors, names(cancer.subtype))
  }
  index <- which(colnames(clinical.data) %in% clinical.factors)
  if(length(index)!=0){
    # remove NA
    clinical.data <-  clinical.data[complete.cases(clinical.data[,index]),]
    return(clinical.data)
  }else{
    return(NULL)
  }
}


# Multivariate Cox regression
CoxMulti <- function(clinical.data, clinical.factors, cancer.subtype=NULL, survival.type){
  # remove NA
  clinical.data <- RmnaBeforeMulCox(clinical.data, clinical.factors, cancer.subtype)
  if(!is.null(cancer.subtype)){
    clinical.factors <- c(clinical.factors, names(cancer.subtype))
  }
  
  # surv model
  clinical.outcome <- ClinicalOutcome(clinical.data, survival.type)
  surv <- Surv(clinical.outcome$time, clinical.outcome$event == 1)
  # Change the name of a multi-factor variable
  index <- which(colnames(clinical.data) %in% clinical.factors)
  colnames(clinical.data)[index] <- paste(clinical.factors,".VS.",sep="")
  # Multivariate Cox regression formula
  multiv_formula <- as.formula(paste("surv~", paste(paste(clinical.factors,".VS.",sep=""), collapse="+"), sep=""))
  coxph.result <- coxph(multiv_formula, data = clinical.data)
  
  # check variable names in coxph result
  names.index <- 1
  for(i in 1:length(index)){
    if(class(clinical.data[,index[i]])=="integer" | class(clinical.data[,index[i]])=="numeric"){
      # continuous variable
      names(coxph.result$coefficients)[names.index] <- clinical.factors[i]
      names.index <- names.index+1
    }else{
      # categorical variable
      all.classes <- unique(clinical.data[,index[i]])
      VS.names <- names(coxph.result$coefficients)[names.index:(names.index+length(all.classes)-2)]
      names(coxph.result$coefficients)[names.index:(names.index+length(all.classes)-2)] <- AdjustName(VS.names,all.classes) 
      names.index <- names.index+length(all.classes)-1
    }
  }
  return(coxph.result)
}

# check variable names in coxph result
AdjustName <- function(VS.names, all.classes){
  others <- unlist(lapply(strsplit(VS.names,split=".VS."),function(x)x[2]))

  ref <- setdiff(all.classes,others)

  adjust.names <- paste(ref, others, sep=".VS.")
  return(adjust.names)
}

CoxphResult <- function(cox.model){
  
  # Extract useful information
  tmp <- summary(cox.model);
  
  # Extract the P value
  p.value <- tmp$coefficients[ ,5]
  #p.value <- round(tmp$coefficients[ ,5], digits = 4)
  #p.value[which(p.value < 0.0001)] <- "<0.0001"
  
  # Extract HR
  HR <- round(tmp$coefficients[ ,2], digits = 4)
  
  # extract the upper and lower bounds of the 95% confidence interval
  HR.confint.lower <- round(tmp$conf.int[ ,"lower .95"], digits = 4)
  HR.confint.upper <- round(tmp$conf.int[ ,"upper .95"],digits = 4)
  
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  
  coxph.result <- as.data.frame(cbind(HR, p.value))
  if(dim(coxph.result)[1]==1){
    rownames(coxph.result) <- names(cox.model$coefficients)
  }
  
  return(coxph.result)
}


#  calculate the C-index for a cox model
cIndex <- function(clinical.data, time, event, models, diff = TRUE){
  
  options(stringsAsFactors = FALSE)

  # load survival, Hmisc package
  suppressPackageStartupMessages(require(survival))
  suppressPackageStartupMessages(require(Hmisc))

  ##-------------------------------
  # CoxphObject：create a coxph object
  CoxphObject <- function(data, time, event, variables){
    formula <- as.formula(paste("Surv(time, event)~", paste(variables, collapse=" + ")))
    coxph.object <-coxph(formula, data = data)

    return(coxph.object)
  }

  ##-------------------------------
  # Cindex and confidence interval of the output
  c.ci <- function(rcorrobj){
    CIndex <- round(rcorrobj['C Index'], digits = 4)
    se     <- rcorrobj['S.D.']/2
    Lower <- round(CIndex - 1.96*se, digits = 4)
    Upper <- round(CIndex + 1.96*se, digits = 4)
    result <- c(CIndex, Lower, Upper)
    names(result) <- c("C-Index", "Lower", "Upper")

    return(result)
  }
  
  # Calculate the C-index value and confidence interval for each model
  coxph.list <- lapply(models, function(x){CoxphObject(data = clinical.data, time = time, event = event, variables = x)})
  pred.models.coxph <- lapply(coxph.list, function(x){predict(x, type = "risk")})
  models.result <- lapply(pred.models.coxph, function(x){rcorr.cens(-x, Surv(time = time, event = event))})
  
  models.filter.result <- lapply(models.result, function(x){c.ci(rcorrobj = x)})
  
  # Whether to compare C-Index1
  if (diff == FALSE) {
    result <- do.call(rbind, models.filter.result)
    conf.interval <- paste(result[, 'Lower'], result[, 'Upper'], sep = "-")
    result <- data.frame(result[, 'C-Index'], conf.interval)
    colnames(result) <- c('C-index', '95%CI')

    return(result)
  } else {
    # The p values for the comparison are all compared to the first model
    compare.cindex <- lapply(pred.models.coxph[-1], function(x){rcorrp.cens(pred.models.coxph[[1]], x, Surv(time = time, event = event))})
    p.value <- c("-", unlist(lapply(compare.cindex, function(x)(round(2*(1 - pnorm(x[['Dxy']] / x[['S.D.']])), digits=4))))) 

    result <- do.call(rbind, models.filter.result)
    conf.interval <- paste(result[, 'Lower'], result[, 'Upper'], sep = "-")
    result <- data.frame(result[, 'C-Index'], conf.interval, p.value)
    colnames(result) <- c('C-index', '95%CI', 'p-value')

    return(result)
  }
}  
