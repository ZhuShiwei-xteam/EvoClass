# This wrapper is designed to work with any single sample data

suppressPackageStartupMessages(library(sequenza))

identify.mut.copy.number.ascat <- function(x,sub.mat.mut, sub.mat.copy) {
  
  mut                  <- sub.mat.mut[x,,drop=FALSE]
  ww                   <- which(as.numeric(sub.mat.copy$Chr)==as.numeric(mut$Chr)
                                &as.numeric(sub.mat.copy$Start)<=as.numeric(mut$Start_position)
                                &as.numeric(sub.mat.copy$End)>=as.numeric(mut$Start_position))
  copy                 <- sub.mat.copy[ww,,drop=FALSE]
  
  mutation_id         <- paste(mut$Patient,mut$Chr,mut$Start_position,mut$Reference,sep=":")
  ref_counts          <- mut$Ref_freq
  var_counts          <- mut$Variant_freq
  normal_cn           <- 2

  patient <- mut$Patient
  Reference_Base <- mut$Reference
  Alternate_Base <- mut$Alternate
  Hugo_Symbol <- mut$Hugo_Symbol
  Variant_Classification <- mut$Variant_Classification
  HGVSp_Short <- mut$HGVSp_Short
  
  if (nrow(copy)!=1) {
    minor_cn <- NA
    major_cn <- NA
    output <-  data.frame(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, patient, Reference_Base, Alternate_Base, Variant_Classification, HGVSp_Short, Hugo_Symbol, stringsAsFactors=FALSE)
    return(output)
  }
  
  minor_cn <- min(c(copy$nA,copy$nB))
  major_cn <- max(c(copy$nA,copy$nB))
  output <-  data.frame(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, patient, Reference_Base, Alternate_Base, Variant_Classification, HGVSp_Short, Hugo_Symbol, stringsAsFactors=FALSE)
  return(output)
}

earlyORlate <- function(patient,complete.mutation.table,purity){
  
  # Get the table ready, with only information for specific patient
  mut.table  <- complete.mutation.table[complete.mutation.table$patient==patient,]
  
  # and all the other stuff
  cellularity <- as.numeric(purity)
  minor.cn <- unlist(mut.table$minor_cn)
  major.cn <- unlist(mut.table$major_cn)
  abs.cn   <- unlist(mut.table$minor_cn) + unlist(mut.table$major_cn)
  depth.t  <- unlist(mut.table$ref_counts) + unlist(mut.table$var_counts)
  max.cn   <- max(abs.cn)
  VAF      <- unlist(mut.table$var_counts)/(unlist(mut.table$var_counts)+unlist(mut.table$ref_counts))
  
  # Estimate theoretical VAFs for each type of copy number
  types    <- mufreq.types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 2)
  types.xy <- mufreq.types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 1)
  types    <- rbind(types, types.xy)
  
  types    <- types[types$Mt >= 1, ]
  types$F  <- 0
  for (i in 1:nrow(types)) {
    types$F[i] <- theoretical.mufreq(cellularity = cellularity, CNn = types$CNn[i], CNt = types$CNt[i], Mt = types$Mt[i])
  }
  
  # Let's create some functions that can estimate whether early or late
  
  get.Mt <- function (F, depth.t, types, CNt, CNn, Mt) {
    types <- types[types$CNn == CNn, ]
    l <- sequenza:::mufreq.dpois(mufreq = F, types$F[types$CNt == CNt&types$Mt<=Mt],depth.t = depth.t)
    l <- l/sum(l)
    L <- data.frame(l = l, Mt = types$Mt[types$CNt == CNt&types$Mt<=Mt])
  }
  
  absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number, locus_minor_cn, gene){
    f.function <- function (c,purity,local.copy.number,locus_minor_cn,gene) {
      if(gene %in% "TP53" & local.copy.number == 2 & locus_minor_cn == 0){
        return((2*purity*c) / (2*(1-purity) + purity*local.copy.number))
      }else{
        return((purity*c) / (2*(1-purity) + purity*local.copy.number))
      }
    }
    x <- dbinom(n.alt,depth, prob=sapply(seq(0.01,1,length.out=100),f.function,purity,local.copy.number,locus_minor_cn,gene))
    if (min(x)==0) {
      x[length(x)] <- 1
    }
    
    names(x) <- seq(0.01,1,length.out=100)
    sub.cint <- function(x, prob = 0.95,n.alt,depth) {
      xnorm   <- x/sum(x)
      xsort   <- sort(xnorm, decreasing = TRUE)
      xcumLik <- cumsum(xsort)
      n = sum(xcumLik < prob) + 1
      LikThresh <- xsort[n]
      cint  <- x[xnorm >= LikThresh]
      all   <- as.numeric(names(x))
      cellu <- as.numeric(names(cint))
      l.t   <- cellu[1]
      r.t   <- cellu[length(cellu)]
      m     <- cellu[which.max(cint)]
      
      prob.subclonal <- sum(xnorm[1:90])
      prob.clonal    <- sum(xnorm[91:100])
      
      data.frame(left = l.t, est = m, right = r.t,prob.subclonal=prob.subclonal,prob.clonal=prob.clonal)
    }
    return(sub.cint(x,n.alt=n.alt,depth=depth))
  }
  
  # add an absolute estimate of the cancer cell fraction
  
  get.all.mut.info <- function(i){
    # First estimate the VAF confidence intervals
    obs.VAF <- VAF[i]
    
    if (abs.cn[i]==0) {
      output <- cbind(obs.VAF,absolute.ccf=NA, absolute.ccf.0.05=NA, absoltue.ccf.0.95=NA, prob.subclonal=NA, prob.clonal=NA)
      return(output)
    }
    
    # Next estimate the likelihood relating to which copy number the mutation has
    L <- get.Mt(F = VAF[i], depth.t = depth.t[i], CNt = abs.cn[i], types = types, CNn = unlist(mut.table$normal_cn[i]), Mt=major.cn[i])
    
    if (is.na(L$l)[1]) {
      output <- cbind(obs.VAF, absolute.ccf=NA, absolute.ccf.0.05=NA, absoltue.ccf.0.95=NA, prob.subclonal=NA, prob.clonal=NA)
      return(output)
    }
    
    # Now determine which likelihood should be used
    absolute.calc        <- absolute.cancer.cell.fraction(n.alt=unlist(mut.table$var_counts)[i],depth=depth.t[i],purity=cellularity,local.copy.number=abs.cn[i],locus_minor_cn = minor.cn[i],gene = unlist(mut.table$Hugo_Symbol)[i])
    absolute.ccf.0.05    <- absolute.calc[1]
    absolute.ccf.0.95    <- absolute.calc[3]
    absolute.ccf         <- absolute.calc[2]
    prob.subclonal       <- absolute.calc[4]
    prob.clonal          <- absolute.calc[5]
    
    # Let's put this all together and output it
    output  <- data.frame(obs.VAF, absolute.ccf, absolute.ccf.0.05, absolute.ccf.0.95, prob.subclonal, prob.clonal, stringsAsFactors=FALSE)
    return(output)
  }
  
  output <- t(sapply(1:nrow(mut.table),get.all.mut.info))
  output <- data.frame(output,stringsAsFactors=FALSE)
  
  colnames(output) <- c('obs.VAF', 'absolute.ccf', 'absolute.ccf.0.05', 'absolute.ccf.0.95', 'prob.subclonal', 'prob.clonal')
  out <- cbind(mut.table,output)
  return(out)
}

clonality.estimation <- function(mutation.table, seg.table, data.type, TCGA.barcode, ANALYSIS.DIR, sub.clonal.cut.off = 1) {

  # create the root directory
  data.folder <- paste(ANALYSIS.DIR,data.type, sep="")

  if ( !file.exists(data.folder)) {
    if (!dir.create(data.folder, recursive = TRUE)) {
      stop("Unable to create root directory.\n")
    }
  }

  # create a folder for the specific patient
  patient.folder <- paste(data.folder, "/",TCGA.barcode,sep="")
  
  if (!file.exists(patient.folder)) {
    if (!dir.create(patient.folder, recursive = TRUE)) {
      stop("Unable to create tmp directory.\n")
    }
  }
  
  # load somatic mutations

  req.mut.colnames <- c('Patient', 'Chr', 'Start_position', 'End_position', 'Reference', 'Alternate', 'Variant_freq', 'Ref_freq', 'Hugo_Symbol', 'Variant_Classification', 'HGVSp_Short')
  
  if (length(req.mut.colnames[!which(req.mut.colnames%in%colnames(mutation.table))])!=0) {
    stop(paste('Mutation table not in correct format\nThe following columns are required:\n', PasteVector(req.mut.colnames, sep="\n"),sep=""))
  }
  
  # load somatic copy numbers
  
  seg.mat.copy <- seg.table
  req.seg.colnames <- c('SampleID', 'Chr', 'Start', 'End', 'nProbes', 'cn', 'nA', 'nB', 'Ploidy', 'Aberrant Cell Fraction')
  
  if (length(req.seg.colnames[!which(req.seg.colnames%in%colnames(seg.mat.copy))])!=0) {
    stop(paste('Segment copy number not in correct format\nThe following columns are required:\n', PasteVector(req.seg.colnames, sep="\n"),sep=""))
  }

  # combine copy number and mutation data
  # first of all, select only samples where both available
  barcodes <- intersect(unique(seg.mat.copy$SampleID),unique(mutation.table[,1]))
  cat(paste("\n\nThere are ", length(barcodes), " patients with both copy number and mutation data. \nOnly these will be used "))
  
  if (!TCGA.barcode%in%barcodes) {
    stop("TCGA barcode not found in either mutation table or seg.mat.copy")
  }
  
  # select patient specific data
  sub.mat.copy  <- seg.mat.copy[seg.mat.copy$SampleID==TCGA.barcode,,drop=FALSE]
  sub.mat.mut <- mutation.table[mutation.table[,1]==TCGA.barcode,,drop=FALSE]
  
  # remove sex chromosomes
  sub.mat.copy <- sub.mat.copy[sub.mat.copy$Chr%in%c(1:22),]
  sub.mat.copy$Chr  <- as.numeric(sub.mat.copy$Chr)
  sub.mat.mut <- sub.mat.mut[sub.mat.mut$Chr%in%c(1:22),]
  sub.mat.mut$Chr <- as.numeric(sub.mat.mut$Chr)

  sub.mat.copy  <- unique(sub.mat.copy) 
  
  #combine mutation and copy number
  mut.table <- data.frame(t(sapply(1:nrow(sub.mat.mut), identify.mut.copy.number.ascat, sub.mat.mut, sub.mat.copy)), stringsAsFactors=FALSE)
  mut.table <- mut.table[!is.na(mut.table$minor_cn),]
  mut.table <- mut.table[!is.na(mut.table$ref_counts),]
  mut.table <- mut.table[!duplicated(mut.table$mutation_id),]
  
  TCGA.purity   <- as.character(unique(sub.mat.copy[,grep("Aberrant",colnames(sub.mat.copy))]))
  
  if (TCGA.purity>1) {
    stop('\nYour tumour content is greater than 1! \nYou probably don\'t have an ASCAT estimate. \nThere is an alernative script for this scenario.')
  }
  
  #### EarlyorLate implementation #####
  TCGA.earlyLate <- earlyORlate(patient=TCGA.barcode,complete.mutation.table=mut.table,purity=TCGA.purity)

  # Let's choose the columns of interest
  TCGA.earlyLate.out <- TCGA.earlyLate
  TCGA.earlyLate.out <- TCGA.earlyLate.out[!is.na(TCGA.earlyLate.out$absolute.ccf.0.95), ]

  # predict the clonal and subclonal mut
  TCGA.earlyLate.out$CI95.timing <- 'Subclonal'
  TCGA.earlyLate.out[as.numeric(unlist(TCGA.earlyLate.out$absolute.ccf.0.95)) >= sub.clonal.cut.off,'CI95.timing'] <- 'Clonal'

  TCGA.earlyLate.out$prob.clonal.timing <- 'Subclonal'
  TCGA.earlyLate.out[as.numeric(unlist(TCGA.earlyLate.out$prob.clonal)) > as.numeric(unlist(TCGA.earlyLate.out$prob.subclonal)),'prob.clonal.timing'] <- 'Clonal'
  col.names <- c('patient', 'mutation_id', 'Reference_Base',  'Alternate_Base', 'Hugo_Symbol', 'Variant_Classification', 'HGVSp_Short', 'ref_counts', 'var_counts', 'obs.VAF', 'normal_cn', 'minor_cn',	'major_cn', 'absolute.ccf', 'absolute.ccf.0.05', 'absolute.ccf.0.95', 'CI95.timing', 'prob.clonal', 'prob.subclonal', 'prob.clonal.timing')

  TCGA.earlyLate.out <- TCGA.earlyLate.out[, col.names]
  
	TCGA.earlyLate <- cbind(TCGA.earlyLate, TCGA.purity)
  
  # let's write the file to a location
  
  earlylate.tsv <- paste(patient.folder,"/",TCGA.barcode,".earlylate.tsv",sep="")
  earlylate.out <- apply(TCGA.earlyLate.out,2,as.character)
  write.table(earlylate.out, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, file=earlylate.tsv)
}

