`check_cntype` <- function(cntype) {
	cnty = -1
	if(cntype=="mixed") { cnty = 1; }
	if(cntype=="matched") { cnty = 2; }
	if(cntype=="mismatched") { cnty = 3; }
	if(cnty==-1) { stop("Invalid cntype!") }
invisible(cnty)
}

`Cbin_extract_item` <- function(samplefile, cbinfile, genome_index_pos, cntype="mixed") {
	samples = read.table(samplefile)
	r = .C("get_cbin_index_pos"
		, "filename"=as.character(cbinfile)
		, "index_pos"=as.integer(genome_index_pos)
		, "number_samples"=as.integer(nrow(samples))
		, "cntype"=as.integer(check_cntype(cntype))
		, "values"=as.double(1:nrow(samples))
		, "PACKAGE" = "CNest")
	r = r$values; names(r)=samples[,1]
invisible(r)
}

`Cbin_extract_range` <- function(samplefile, cbinfile, start, stop, cntype="mixed") {
	samples = read.table(samplefile)
	number_positions = (stop-start)+1
	long_values = as.double(1:(nrow(samples)*number_positions))
	r = .C("get_cbin_index_range"
		, "filename"=as.character(cbinfile)
		, "index_pos"=as.integer(start)
		, "number_samples"=as.integer(nrow(samples))
		, "number_postions"=as.integer(number_positions)
		, "cntype"=as.integer(check_cntype(cntype))
		, "values"=long_values
		, "PACKAGE" = "CNest")$values
invisible(r)
}

### cnv gwas
`cnv_gwas` <- function(sample_file, phenotype_file, covariates_file, cbin_file, outputname, test_type="linear", bin_type="mixed", invR=FALSE, chunk_n, chunk_len) {
	my.invnorm = function(x) {
	   res = rank(x); res = qnorm(res/(length(res)+0.5))
	return(res)
	}
	samples = read.table(sample_file)
	n_samples = nrow(samples)
	phenotypes_all = read.table(phenotype_file, header=T)
	covars_all = read.table(covariates_file, header=T)
	phenotypes_all = phenotypes_all[phenotypes_all$sample_id%in%samples[,1],]
	phenotypes_all = phenotypes_all[phenotypes_all$sample_id%in%covars_all$sample_id,]
	covars_all = covars_all[covars_all$sample_id%in%phenotypes_all$sample_id,]
	phenotypes_all = phenotypes_all[order(phenotypes_all$sample_id),]; covars_all = covars_all[order(covars_all$sample_id),];
	index_file = gsub(".cbin", ".cbii", cbin_file)
	index = read.table(index_file)
	nposes = ceiling(nrow(index)/chunk_len)
	stop_pos = min((nposes*chunk_n), nrow(index))
	start_pos =  (stop_pos- nposes)+1
	genome_pos = index[start_pos:stop_pos,]
	cndata = Cbin_extract_range(sample_file, cbin_file, start_pos, stop_pos, bin_type)	
	for(x in 2:ncol(phenotypes_all)) {
		print(paste("processing..", colnames(phenotypes_all)[x]))
		chr = unlist(strsplit(gsub(".cbin", "", basename(cbin_file)), "_"))
		outname = paste0(outputname, "_chunk", chunk_n, "_", chr[length(chr)], "_", colnames(phenotypes_all)[x], ".cwas", sep="")
		phenotypes = data.frame(phenotypes_all[,1], phenotypes_all[,x])
		colnames(phenotypes) = c(colnames(phenotypes_all)[c(1,x)])
		phenotypes = phenotypes[complete.cases(phenotypes),]
		covars = covars_all[covars_all$sample_id%in%phenotypes[,1],]; phenotypes = phenotypes[phenotypes[,1]%in%covars$sample_id,];
		phenotypes= phenotypes[order(phenotypes[,1]),]; covars = covars[order(covars$sample_id),];
		y1 = phenotypes[,2]
		X = as.matrix(covars[,-(1)])
		bs_ps_list = list()
		if(invR==TRUE) { y1 = my.invnorm(y1); }
		for(y in 1:nrow(genome_pos)) {
			s_pos = ((n_samples*y)-n_samples)+1
			S = cndata[s_pos:(s_pos+n_samples)]
			names(S) = samples[,1]
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
		genome_pos$n = length(y1)
		cwas=cbind(genome_pos, do.call(rbind, bs_ps_list))
		colnames(cwas) = c("chr", "start", "stop", "n", "mad", "mean", "beta", "p")
		write.table(cwas, file=outname, sep="\t", row.names=F, quote=F)
		print(paste(colnames(phenotypes_all)[x], ".. done"))
	}
}


