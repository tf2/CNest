`check_cntype` <- function(cntype) {
	cnty = -1
	if(cntype=="mixed") { cnty = 1; }
	if(cntype=="matched") { cnty = 2; }
	if(cntype=="mismatched") { cnty = 3; }
	if(cnty==-1) { stop("Invalid cntype!") }
return(cnty)
}

`Cbin_extract_item` <- function(samplefile, cbinfile, genome_index_pos, cntype="mixed") {
	samples = read.table(samplefile)
	r = .C("get_cbin_index_pos"
		, "filename"=as.character(cbinfile)
		, "index_pos"=as.integer(genome_index_pos)
		, "number_samples"=as.integer(nrow(samples))
		, "cntype"=as.integer(check_cntype(cntype))
		, "values"=as.double(1:nrow(samples))
		, "PACKAGE" = "Rbin")
	r = r$values; names(r)=samples[,1]
return(r)
}

`Cbin_extract_range` <- function(samplefile, cbinfile, start, stop, cntype="mixed") {
return( lapply(start:stop, function(x) { Cbin_extract_item(samplefile, cbinfile, x, cntype) }) )
}


### cnv gwas
`cnv_gwas` <- function(sample_file, phenotype_file, covariates_file, cbin_file, outputname, test_type="linear", bin_type="mixed", invR=FALSE, chunk_n, chunk_len) {
	my.invnorm = function(x) {
	   res = rank(x); res = qnorm(res/(length(res)+0.5))
	return(res)
	}
	samples = read.table(sample_file)
	phenotypes = read.table(phenotype_file, header=T)
	covars = read.table(covariates_file, header=T)
	phenotypes = phenotypes[phenotypes$sample_id%in%covars$sample_id,]
	covars = covars[covars$sample_id%in%phenotypes$sample_id,]
	phenotypes = phenotypes[order(phenotypes$sample_id),]; covars = covars[order(covars$sample_id),];
	index_file = gsub(".cbin", ".cbii", cbin_file)
	index = read.table(index_file)
	nposes = ceiling(nrow(index)/chunk_len)
	stop_pos = min((nposes*chunk_n), nrow(index))
	start_pos =  (stop_pos- nposes)+1
	genome_pos = index[start_pos:stop_pos,]
	cndata = lapply(start_pos:stop_pos, function(x) { Cbin_extract_item(sample_file, cbin_file, x, bin_type) }) 
	X = as.matrix(covars[,-(1)])	
	for(x in 2:ncol(phenotypes)) {
		print(paste("processing..", colnames(phenotypes)[x]))
		chr = unlist(strsplit(gsub(".cbin", "", basename(cbin_file)), "_"))
		outputname = paste0(outputname, "_chunk", chunk_n, "_", chr[length(chr)], "_", colnames(phenotypes)[x], ".cwas", sep="")
		y1 = phenotypes[,x]
		bs_ps_list = list()
		if(invR==TRUE) { y1 = my.invnorm(y1); }
		for(y in 1:length(cndata)) {
			S = cndata[[y]]
			S = S[names(S)%in%phenotypes$sample_id]
			S = S[order(names(S))]
			if(sum(names(S)==phenotypes$sample_id)!=length(S) | sum(names(S)==covars$sample_id)!=length(S)) { stop("sample id problem") }
			beta_and_p = NULL
			if(test_type=="logistic") {
				beta_and_p = summary(glm(y1~S+X,family = binomial(link="logit"), na.action(na.omit)))$coefficients
			}
			if(test_type=="linear") {
				beta_and_p = summary(glm(y1~S+X))$coefficients
			}
			bs_ps_list[[y]] = c(mad(S), mean(S), beta_and_p[2,c("Estimate", "Pr(>|t|)")])
		}
		cwas=cbind(genome_pos, do.call(rbind, bs_ps_list))
		colnames(cwas) = c("chr", "start", "stop", "n", "mad", "mean", "beta", "p")
		write.table(cwas, file=outputname, sep="\t", row.names=F, quote=F)
		print(paste(colnames(phenotypes)[x], ".. done"))
	}
}


