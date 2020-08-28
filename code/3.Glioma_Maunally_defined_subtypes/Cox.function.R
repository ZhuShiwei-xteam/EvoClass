

Cox.function <- function(time, event, clinical.data, clinical.variate = NULL)
{
  ###Set up the work environment
  options(stringsAsFactors = FALSE, warn = -1);
  suppressPackageStartupMessages(require(survival));

  ###Determine covariate types: numeric (Num. covariate), non-numeric (Chara. covariate).
  if(is.null(clinical.variate))
  {
      covariates  <- colnames(clinical.data)[-c(1:3)];
  }
  if(is.numeric(clinical.variate))
  {
      covariates  <- colnames(clinical.data)[clinical.variate];
  }
  num.variate <- NULL
  for(i in covariates)
  {
      if(is.numeric(clinical.data[, i]))
      {
          num.variate <- append(num.variate, i)
      }
  }
  chara.variate <- setdiff(covariates, num.variate)

  #### univariate.cox
  univariate.cox <- function(data, num, chara)
  {

    univ_formulas <- sapply(covariates,
                            function(x) as.formula(paste('Surv(time, event)~', x)));
    
  
    univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data)});


    univ_results <- lapply(univ_models, function(x)
                            {                             
                             tmp <-summary(x);
                             
                             #Extract the P value, reserving two significant digits
                             p.value <- round(tmp$coefficients[ ,5], digits = 4);
                             p.value[which(p.value < 0.0001)] <- "<0.0001";
   
                             
                             #Extract HR
                             HR <- round(tmp$coefficients[ ,2], digits = 4);
                             
                             # extract the upper and lower bounds of the 95% confidence interval
                             HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 4);
                             HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 4);    

                             
                             HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")");
                             
                             variate <- rownames(tmp$coefficients);
                             
                             #Combine all values in one matrix
                             all.data <- as.data.frame(cbind(variate, HR, p.value));
                           }
                          );
    
    #Standardized output format
    for(i in num)
    {
        tmp <- univ_results[[i]];
        tmp$type <- " ";
        univ_results[[i]] <- tmp; 
    }

    for(i in chara)
    {
        tmp <- univ_results[[i]]
        tmp.variate <- substr(tmp$variate, start = 0, stop = nchar(i));
        tmp.type <- sub(pattern = i, replacement = "", tmp$variate);
        tmp$variate <- tmp.variate;
        tmp$type <- paste(tmp.type, "VS", setdiff(unique(data[, i]), tmp.type));
        univ_results[[i]] <- tmp; 
    }

    #converte list to data.frame
    univ.result <- do.call(rbind.data.frame, univ_results);
    univ.result <- univ.result[,c(1,4,2,3)];
    colnames(univ.result) <- c('variate', 'type', 'univ HR (95% CI for HR)', 'univ p value')
    rownames(univ.result) <- NULL;
    return(univ.result)  
  };

  ####multivariate.cox
  multivariate.cox <- function(data, num, chara)
  {
    options(stringsAsFactors = FALSE);

    #All the selected covariates were analyzed by multiple factors directly
    multiv_formula <- as.formula(paste("Surv(time, event)~", paste(covariates, collapse="+"), sep=""));
    multiv_model    <- coxph(multiv_formula, data = data);
    
    tmp <- summary(multiv_model);
    
    #Extract the P value, reserving two significant digits
    p.value <- round(tmp$coefficients[ ,5], digits = 4);
    p.value[which(p.value < 0.0001)] <- "<0.0001"; 
    
       
    #Extract HR
    HR <- round(tmp$coefficients[ ,2], digits = 4);
    
    #extract the upper and lower bounds of the 95% confidence interval
    HR.confint.lower <- round(tmp$conf.int[ ,"lower .95"], digits = 4);
    HR.confint.upper <- round(tmp$conf.int[ ,"upper .95"],digits = 4);
    

    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")");
    
    variate <- rownames(tmp$coefficients);
    

    multiv_result <- as.data.frame(cbind(variate, HR, p.value));

    #STEP3:Create a new data frame to store multivariate cox results
    multiv.result <- NULL;

    for(i in num)
    {
        n.row <- grep(pattern = i, multiv_result$variate);
        tmp <- multiv_result[n.row, ];
        tmp$type <- " ";
        multiv.result <- rbind(multiv.result,tmp);
    }

    for(i in chara)
    {
        n.row <- grep(pattern = i, multiv_result$variate);
        tmp <- multiv_result[n.row, ];
        tmp.variate <- substr(tmp$variate, start = 0, stop = nchar(i));
        tmp.type <- sub(pattern = i, replacement = "", tmp$variate);
        tmp$variate <- tmp.variate;
        tmp$type <- paste(tmp.type, "VS", setdiff(unique(data[, i]), tmp.type));
        multiv.result <- rbind(multiv.result,tmp);         
    }
    multiv.result <- multiv.result[,c(1,4,2,3)]
    colnames(multiv.result) <- c("variate", "type","multiv HR (95% CI for HR)", "multiv p value");
    rownames(multiv.result) <- NULL;

    return(multiv.result);
  };

  UniCoxPH <- univariate.cox(data = clinical.data, num = num.variate, chara = chara.variate);
  MultiCoxPH <- multivariate.cox(data = clinical.data, num = num.variate, chara = chara.variate);


  cox.result <- merge(UniCoxPH, MultiCoxPH, by = c("variate", "type"), all = T);
  colnames(cox.result) <- c("variate", " ", "univ HR (95% CI for HR)", "univ p value","multiv HR (95% CI for HR)", "multiv p value");


  for(i in chara.variate)
  {
  tmp.row <- which(cox.result[,1] == i);
  tmp.vec <- c(i, rep(" ", times = 5));
  cox.result <- rbind(cox.result[1 : (tmp.row-1),], tmp.vec, cox.result[tmp.row : nrow(cox.result),]); 
  };
  cox.result[duplicated(cox.result[,1]),1] <- " ";
  rownames(cox.result) <- 1:nrow(cox.result);

  return(cox.result)
}


