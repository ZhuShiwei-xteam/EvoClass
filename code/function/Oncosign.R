# Main function to run OncoSign and SEA together
oncosign.sea = function(cases=NULL,gam=NULL,max.levels=3,min.alterations=0.01){
	
	# property_file
	# if property file is not provided, build it from arguments
	# if(missing(property_file)){
		
		if(is.null(cases) | is.null(gam)){
			stop("ERROR: Either the property file or cases and gam attributes need to be provided .\n")
		}
		
		props = c("case_file","genomic_alteration_file","max_levels","min_fraction_of_alterations",
				cases,gam,max.levels,min.alterations)		
		P=matrix(props,ncol=2)
		write.table(P,file="oncosign.props",col.names=F,row.names=F,quote=F,sep="=")
		property_file = "oncosign.props"
	# }else{
		# props = read.table(property_file,header=F,sep="=")
		# P = as.matrix(props)
	# }
	
	gam.file = P[P[,1]=="genomic_alteration_file",2]
	cases.file = P[P[,1]=="case_file",2]
	max.levels = P[P[,1]=="max_levels",2]
	min.alterations = P[P[,1]=="min_fraction_of_alterations",2]

	# run OncoSign
	run.oncosign(property_file)

	# load data
	data=read.table(gam.file,header=T,row.names=1,check.names=F)
	classes = read.table("oncosign_results/osc_assignments.txt",sep="\t",row.names=1)
	
	# initialize result object with input data and osc assignments
	os.list = alist()
	os.list[[1]] = data
	os.list[[2]] = classes
	
	actual.levels = ncol(classes)
	
	for(i in 1:actual.levels){
		sea.results = sea(data,classes,level = i)
		os.list[[i+2]] = sea.results$signatures
	}
	
	names(os.list) = c("data","classes",paste("level",1:max.levels,sep=""))
	return(os.list)

}

# Main function to run OncoSign (to determine enriched signatures, SEA needs to be run separately using run.sea)
oncosign = function(property_file,cases=NULL,gam=NULL,max.levels=3,min.alterations=0.01,K=0){

	#if property file is not provided, build it from arguments
	if(missing(property_file)){
		
		if(is.null(cases) | is.null(gam)){
			stop("ERROR: Either the property file or cases and gam attributes need to be provided .\n")
		}
		
		props = c("case_file","genomic_alteration_file","max_levels","min_fraction_of_alterations","scaling_coeff",
				cases,gam,max.levels,min.alterations,K)		
		P=matrix(props,ncol=2)
		write.table(P,file="oncosign.props",col.names=F,row.names=F,quote=F,sep="=")
		property_file = "oncosign.props"
	}
	else{
		props = read.table(property_file,header=F,sep="=")
		P = as.matrix(props)
	}

	gam.file = P[P[,1]=="genomic_alteration_file",2]
	cases.file = P[P[,1]=="case_file",2]
	max.levels = P[P[,1]=="max_levels",2]
	min.alterations = P[P[,1]=="min_fraction_of_alterations",2]

	# run OncoSign
	run.oncosign(property_file)

	# load data
	data=read.table(gam.file,header=T,row.names=1,check.names=F)
	classes = read.table("oncosign_results/osc_assignments.txt",sep="\t",row.names=1)

	# create and return results object	
	os.list = list("data" = data, "classes" = classes)
	
	return(os.list)

}

# R wrapper for Java function
run.oncosign = function(property_file){
	Rpaths = paste(.libPaths(),"oncosign",sep="/")	
	OSpath = Rpaths[which(file.exists(Rpaths))]
	jar.file = paste(OSpath,"/exec/oncosign.jar",sep="")
	command=paste("java -jar",jar.file,property_file)
	system(command)
}

# R function for SEA (Signature Enrichment Analysis)
sea = function(data,classes,level = 1, min.alterations = 0.01,
					max.alterations=1, print.results = TRUE, maxQ = 0.05){

	# by default the first column of classes are the sample names
	class = as.matrix(classes[,level])
	rownames(class) = rownames(classes)
	class = as.matrix(class[class[,1] != "",])
	colnames(class) = c("class")

	MIN_ALTERATIONS = nrow(class)*min.alterations
	MAX_ALTERATIONS = nrow(class)*max.alterations
	
	data = t(data)
	
	merged.data = merge(data,class,by="row.names",all.x=FALSE,all.y=FALSE)

	if(nrow(merged.data) == 0){
		data = t(data)
		merged.data = merge(data,class,by="row.names",all.x=FALSE,all.y=FALSE)	
		if(nrow(merged.data) == 0){
			stop("ERROR: Data ID and Classes ID do not match.\n")
		}
	}
	
	merged.data.rownames = merged.data$Row.names
	merged.data.colnames = colnames(merged.data)

	class.size=table(merged.data$class)
	all=sum(class.size)
	prob = class.size/all
	classes = unique(sort(merged.data$class))
	class.names = sapply(classes, as.character)

	entries = c()
	sea.p.values = c()
	gof.p.values = c()
	max.class = c()

	no_of_events = 0
	event_ids = c()

	sel=c()

	for(i in 2:(ncol(merged.data)-1)){
	
		if(sum(merged.data[,i]) > MIN_ALTERATIONS && sum(merged.data[,i]) < MAX_ALTERATIONS+1){
	
			#cat(c(colnames(merged.data)[i],"\n"))
	
			sel = c(sel,TRUE)	
			mut=c()
			for(j in 1:length(classes))
				mut=c(mut,length(merged.data[,i][merged.data[,i] == 1 & merged.data$class == classes[j]]))
		
			entries=c(entries,mut)
		
			#gof.test = chisq.test(mut,p=prob,simulate.p.value = TRUE)
			gof.test = suppressWarnings(chisq.test(mut,p=prob))
			fraction = mut/class.size
			max.class.index = which.max(fraction)
		
			f.values = c(mut[max.class.index], (class.size[max.class.index] - mut[max.class.index]),
				sum(merged.data[,i]) - mut[max.class.index],
				all - sum(merged.data[,i]) - class.size[max.class.index] + mut[max.class.index]) 
		
			sea.test = fisher.test(matrix(f.values,ncol=2),alternative="greater")
			
			sea.p.values=c(sea.p.values,sea.test$p.value)
			gof.p.values=c(gof.p.values,gof.test$p.value)
			max.class=c(max.class,class.names[max.class.index],class.size[max.class.index])
	
			no_of_events = no_of_events+1
			event_ids = c(event_ids,merged.data.colnames[i])
		}
		else
			sel = c(sel,FALSE)

	}	

	compute.q.values = qvalue(sea.p.values,lambda=0)
	sea.q.values = compute.q.values$qvalues

	E = matrix(entries,nrow=no_of_events,byrow=T)
	MC = matrix(max.class,ncol=2,byrow=T)

	results.colnames = c("Event",paste(classes,sep="\t"),"GOF p.val","SEA p.val","SEA q.val","MaxClass","nsamples")
	results=cbind(event_ids,E,gof.p.values,sea.p.values,sea.q.values,MC)
	
	colnames(E) = classes
	rownames(E) = event_ids
	
	significant = data.frame(sea.p.values,sea.q.values,MC)
	
	colnames(significant) = c("p","q","enriched.osc","nsamples")
	rownames(significant) = event_ids
	
	E = E[order(significant[,"q"]),]
	S = significant[order(significant[,"q"]),]
	
	#sea.results = list("events" = event_ids, "occurrences" = E, "gof.p" = gof.p.values, "significant" = S)
	
	sea.results = merge(E,S,by="row.names",all.x=FALSE,all.y=FALSE)
	sea.results = sea.results[order(sea.results[,"q"]),]
	rownames(sea.results) = sea.results[,1]
	sea.results = sea.results[,-1]

	# here you assemble all the material to provide oncogenic signature table
	# <OSC_ID> <event: q < maxQ in OSC_ID> <% altered in OSC_ID> <% of alteration in event>

	maxcol=c()
	for(i in 1:nrow(sea.results)){
		value = sea.results[i,sea.results[,"enriched.osc"][i]]
		maxcol=c(maxcol,value)
	}
	
	T = cbind(sea.results,maxcol)
	event.size = as.numeric(rowSums(T[,1:4]))
	T = cbind(T,event.size)
	csize = as.numeric(levels(T[,"nsamples"]))[T[,"nsamples"]]
	class.percent = round(100*T$maxcol/csize,digits=2)
	event.percent = round(100*T$maxcol/T$event.size,digits=2)
	T = cbind(T,class.percent,event.percent)

	if(print.results){
		file.name = paste("oncosign_results/sea_report_L",level,".txt",sep="")
		write.table(T,file=file.name,sep="\t",row.names=F,quote=F)
	}
		
	table = T[T$q < maxQ,]
	osc = as.character(table$enriched.osc)
	table = cbind(osc,rownames(table),table$class.percent,table$event.percent)
	colnames(table) = c("OSC_ID","Event","OSC.Percent","Event.Percent")
	table = table[order(table[,"OSC_ID"]),]
	table = as.data.frame(table, stringsAsFactors = FALSE)
	table[] = lapply(table,type.convert)

	results = list("summary" = sea.results, "signatures" = table)

	return(results)
	
}