`RbinConvert_aCGH` <- function(inputname, binname) {
	.C("makeaCGHbin", "inputname"=as.character(inputname), "biname"=as.character(binname),"PACKAGE" = "Rbin")
}

`RbinConvert_exome` <- function(inputname, binname) {
	.C("makeexomebin", "inputname"=as.character(inputname), "biname"=as.character(binname),"PACKAGE" = "Rbin")
}

`RbinConvert_exome_ratio` <- function(inputname, binname) {
    .C("make_exome_ratio_reference_bin", "inputname"=as.character(inputname), "biname"=as.character(binname),"PACKAGE" = "Rbin")
}

`RbinRead_aCGH` <- function(chr, start, stop, indexfile = NULL, filenames=NULL, verbose=TRUE) {
	
	offsets <- .C( "getOffsets",
			"filename" = as.character(indexfile),
			"chr" = as.integer(chr), 
	  		"start" =as.integer(start),
	   		"stop" = as.integer(stop),
	    	"rows" = as.integer(0),
	    	"offset" = as.double(0),
		"PACKAGE" = "Rbin"
			)
	if(verbose) {
		print("Calculated Offsets and Length!!!")
	}
	positions <- .C( "getPositions",
			"filename" = as.character(indexfile),
			"chr" = as.integer(chr), 
	  		"start" =as.integer(start),
	   		"stop" = as.integer(stop),
	    	"chr_values" = as.integer(1:offsets$rows),
	    	"start_values" = as.integer(1:offsets$rows),
	    	"stop_values" = as.integer(1:offsets$rows),
		"PACKAGE" = "Rbin"
			)
	if(verbose) {
		print("Retrived Genomic Positions!!!")
	}
	matr = matrix(nrow=offsets$rows, ncol=length(filenames))
	if(verbose) {
		pb <- txtProgressBar(min = 0, max = length(filenames), style = 3);
	}
	for(x in 1:length(filenames)) {
		rr= .C( "getValues",
			"filename" = as.character(filenames[x]),
	 		"start_position" = as.double(offsets$offset), 
	  		"number_row_to_read" =as.integer(offsets$rows),
	   		"values" = as.double(1:offsets$rows),
			"PACKAGE" = "Rbin"
			)
			matr[,x] = rr$values; 
			if(verbose) {
				setTxtProgressBar(pb, x); rr=NULL;
			}
	   } 
	   data = cbind(positions$chr_values, positions$start_values, positions$stop_values, matr)
	   colnames(data)=c("chr", "start", "stop", filenames)
		if(verbose) {   
			close(pb);
		}
	if(verbose) {
		print(paste("Retirved ", dim(data)[1]*dim(data)[2], " datapoints with ", length(filenames), " sequential file reads!!!",sep=""))
	}
	
	positions=NULL; offsets=NULL; rr=NULL;
	invisible(data)
}

`RbinRead_exome` <- function(chr, start, stop, typ="adm3score", indexfile = NULL, filenames=NULL) {
	
	offsets <- .C( "geteOffsets",
				  "filename" = as.character(indexfile),
				  "chr" = as.integer(chr), 
				  "start" =as.integer(start),
				  "stop" = as.integer(stop),
				  "rows" = as.integer(0),
				  "offset" = as.double(0),
				"PACKAGE" = "Rbin"
				  )
	print("Calculated Offsets and Length!!!")
	
	positions <- .C( "getePositions",
					"filename" = as.character(indexfile),
					"chr" = as.integer(chr), 
					"start" =as.integer(start),
					"stop" = as.integer(stop),
					"chr_values" = as.integer(1:offsets$rows),
					"start_values" = as.integer(1:offsets$rows),
					"stop_values" = as.integer(1:offsets$rows),
					"PACKAGE" = "Rbin"
					)
	print("Retrived Genomic Positions!!!")
	
	matr = matrix(nrow=offsets$rows, ncol=length(filenames))
	pb <- txtProgressBar(min = 0, max = length(filenames), style = 3);
	
	for(x in 1:length(filenames)) {
		rr = NULL
		if(typ=="adm3score") {
			rr= .C( "geteValues2",
			   "filename" = as.character(filenames[x]),
			   "start_position" = as.double(offsets$offset), 
			   "number_row_to_read" =as.integer(offsets$rows),
			   "values" = as.double(1:offsets$rows),
				"PACKAGE" = "Rbin"
			   )
			
		} else {
			rr= .C( "geteValues1",
				   "filename" = as.character(filenames[x]),
				   "start_position" = as.double(offsets$offset), 
				   "number_row_to_read" =as.integer(offsets$rows),
				   "values" = as.double(1:offsets$rows),
				"PACKAGE" = "Rbin"
				   )
			
		}
		matr[,x] = rr$values; setTxtProgressBar(pb, x); rr=NULL;
	} 
	data = cbind(positions$chr_values, positions$start_values, positions$stop_values, matr)
	colnames(data)=c("chr", "start", "stop", filenames)
	close(pb);
	print(paste("Retirved ", dim(data)[1]*dim(data)[2], " datapoints with ", length(filenames), " sequential file reads!!!",sep=""))
	positions=NULL; offsets=NULL; rr=NULL;
	invisible(data)
}


`RbinRead_exome_ratio` <- function(chr, start, stop, typ="mixed", indexfile = NULL, filenames=NULL) {
	
	offsets <- .C( "geteOffsets_ratio",
				  "filename" = as.character(indexfile),
				  "chr" = as.integer(chr), 
				  "start" =as.integer(start),
				  "stop" = as.integer(stop),
				  "rows" = as.integer(0),
				  "offset" = as.double(0),
				"PACKAGE" = "Rbin"
				  )
	print("Calculated Offsets and Length!!!")
	
	positions <- .C( "getePositions_ratio",
					"filename" = as.character(indexfile),
					"chr" = as.integer(chr), 
					"start" =as.integer(start),
					"stop" = as.integer(stop),
					"chr_values" = as.integer(1:offsets$rows),
					"start_values" = as.integer(1:offsets$rows),
					"stop_values" = as.integer(1:offsets$rows),
					"PACKAGE" = "Rbin"
					)
	print("Retrived Genomic Positions!!!")
	
	matr = matrix(nrow=offsets$rows, ncol=length(filenames))
	pb <- txtProgressBar(min = 0, max = length(filenames), style = 3);
	
	for(x in 1:length(filenames)) {
		rr = NULL
		if(typ=="mixed") {
			rr= .C( "geteValues1_ratio",
			   "filename" = as.character(filenames[x]),
			   "start_position" = as.double(offsets$offset), 
			   "number_row_to_read" =as.integer(offsets$rows),
			   "values" = as.double(1:offsets$rows),
				"PACKAGE" = "Rbin"
			   )
		}
		if(typ=="matched") {
			rr= .C( "geteValues2_ratio",
				   "filename" = as.character(filenames[x]),
				   "start_position" = as.double(offsets$offset), 
				   "number_row_to_read" =as.integer(offsets$rows),
				   "values" = as.double(1:offsets$rows),
				"PACKAGE" = "Rbin"
				   )
			
		}
		if(typ=="mismatched") {
			rr= .C( "geteValues3_ratio",
				   "filename" = as.character(filenames[x]),
				   "start_position" = as.double(offsets$offset), 
				   "number_row_to_read" =as.integer(offsets$rows),
				   "values" = as.double(1:offsets$rows),
				"PACKAGE" = "Rbin"
			)
			
		}


		matr[,x] = rr$values; setTxtProgressBar(pb, x); rr=NULL;
	} 
	data = cbind(positions$chr_values, positions$start_values, positions$stop_values, matr)
	colnames(data)=c("chr", "start", "stop", filenames)
	close(pb);
	print(paste("Retirved ", dim(data)[1]*dim(data)[2], " datapoints with ", length(filenames), " sequential file reads!!!",sep=""))
	positions=NULL; offsets=NULL; rr=NULL;
	invisible(data)
}

processtobin <- function(inputfile, tempfile, binfile, indexfile) {
		data = read.table(inputfile, sep="\t", skip=1)
		# split index
		index = read.table(indexfile)
		parts = unlist(strsplit(as.character(index[,1]), ":"))
		chr = parts[seq(1,length(parts), by=2)]
        chr = gsub("chr", "", chr)
		pos = unlist(strsplit(parts[seq(2,length(parts), by=2)], "-"))
		start = as.numeric(pos[seq(1, length(pos), by=2)])
		end = as.numeric(pos[seq(2, length(pos), by=2)])
		# output data
		data = data.frame(chr, start, end, data)
		data = data.frame(data, 1,1)
		chr = as.character(data[,1])
		chr[chr=="X"] = "23"
		chr[chr=="Y"] = "24"
		data[,1] = as.numeric(chr)
		data = data[order(data[,1], data[,2], data[,3]),]
		write.table(data,  file=tempfile, sep="\t", row.names=F, col.names=F, quote=F)
		rname = RbinConvert_aCGH(tempfile, binfile)
		system(paste("rm ", tempfile, sep=""))
}

processtobin_fast <- function(tempfile, binfile) {
	# For using with fast mode
	rname = RbinConvert_aCGH(tempfile, binfile)
}

### Gender classification
classify_gender <- function(bin_dir, index_file, qc_file) {
	ids = dir(bin_dir)
	filenames = paste(bin_dir, "/", ids, sep="")
	index = read.table(index_file)

	x_med = vector()
	x_mad = vector()
	mads = vector()
	mens = vector()
	sds = vector()
	co = vector()
	for(x in 1:length(filenames)) {
		rr2= .C( "getValues",
				"filename" = as.character(filenames[x]),
		 		"start_position" = as.double(0), 
		  		"number_row_to_read" =as.integer(nrow(index)),
		   		"values" = as.double(1:nrow(index)),
		   		"PACKAGE" = "Rbin"
				)$values

				mens[x] = mean(rr2)
				sds[x] = sd(rr2)
				mads[x] = mad(rr2)
				x_med[x] = mean(rr2[index[,1]==23])
				x_mad[x] = mad(rr2[index[,1]==23])
				print(x)
	}

	dat = data.frame(ids, mens, sds, mads, x_med, x_mad)
	colnames(dat) = c("sample", "mean", "sd", "mad", "Xmean", "Xmad")
	write.table(dat, file=qc_file, sep="\t", row.names=F, quote=F)
}

assign_gender <- function(qc_file, gender_file) {
	qc = read.table(qc_file, header=T)
	kk = qc$mean/qc$Xmean
	k = kmeans(kk, 2)
	gender = as.character(k$cluster)
	wmale = names(which.min(k$centers[,1]))
	gender[gender==wmale]  = "female"
	gender[gender!="female"]  = "male"
	gen_dat = data.frame(qc[,1], gender)
	colnames(gen_dat) = c("sample",	"gender")
	write.table(gen_dat, file=gender_file, sep="\t", row.names=F, quote=F)
}

generate_correlation <- function(bin_dir, cor_dir, sample_name, index_file) {
	index = read.table(index_file)
	filenames = dir(bin_dir, full.names=TRUE)
	cor_file = paste0(cor_dir, "/", sample_name)
	rr1= .C( "getValues",
				"filename" = paste0(bin_dir, '/', sample_name),
		 		"start_position" = as.double(0), 
		  		"number_row_to_read" =as.integer(nrow(index)),
		   		"values" = as.double(1:nrow(index)),
		   		"PACKAGE" = "Rbin"
				)$values
	co = vector()
	for(x in 1:length(filenames)) {
		rr2= .C( "getValues",
				"filename" = filenames[x],
		 		"start_position" = as.double(0), 
		  		"number_row_to_read" =as.integer(nrow(index)),
		   		"values" = as.double(1:nrow(index)),
		   		"PACKAGE" = "Rbin"
				)$values
		co[x] = cor(rr1, rr2)
	}
	dataset = data.frame(basename(filenames), co)
	write.table(dataset, file=cor_file, sep="\t", row.names=F, col.names=F, quote=F)
}

generate_correlation_fast <- function(bin_dir, cor_dir, index_file) {
	# read everything only once but require big mem
	index = read.table(index_file)
	filenames = dir(bin_dir, full.names=TRUE)
	rr_mat = matrix(0, nrow=nrow(index), ncol=length(filenames))
	for(x in 1:length(filenames)) {
		rr_mat[, x] = .C( "getValues",
				"filename" = filenames[x],
		 		"start_position" = as.double(0), 
		  		"number_row_to_read" =as.integer(nrow(index)),
		   		"values" = as.double(1:nrow(index)),
		   		"PACKAGE" = "Rbin"
				)$values
	}
	colnames(rr_mat) = basename(filenames)
	cor_mat = cor(rr_mat)
	# TODO: output to each file
	dataset = data.frame(filenames, co)
	write.table(dataset, file=cor_file, sep="\t", row.names=F, col.names=F, quote=F)
}

get_references <- function(sample_name, index_file, gender_file, outputdir,
							correlation_path, bin_dir, batch_size = 1000) {
	get_sample <- function(filename, index) {
	return(rr1= .C( "getValues",
				"filename" = as.character(filename),
		 		"start_position" = as.double(0), 
		  		"number_row_to_read" =as.integer(nrow(index)),
		   		"values" = as.double(1:nrow(index)),
		   		"PACKAGE" = "Rbin"
				)$values)
	}
	mean_norm <- function(ref_filenames, index) {
		sum_cov = rep(0, nrow(index))
		for(x in 1:length(ref_filenames)) {
			cov = get_sample(ref_filenames[x], index)
			sum_cov = sum_cov + cov
		}
	return(sum_cov)
	}
	med_norm <- function(ref_filenames, index, index_file) {
		meds = NULL
		chrs = unique(index[,1])
		for(x in 1:length(chrs)) {
			r = RbinRead_aCGH(chrs[x], min(index[index[,1]==chrs[x],2]), max(index[index[,1]==chrs[x],3]), index_file, ref_filenames)
			m = apply(r[,-(1:3)], 1, median)
			meds = c(meds, m)
		}
	return(meds)
	}
	get_reference <- function(sample_name, sample_file, index, index_file, correlation_path, bin_dir, gender_file, batch_size = 5000, cor_cut = 0.9) {
		# sample correlation file processing
		# V1 is sample name; V2 is cor coeff
		d_corl = read.table(paste0(correlation_path, "/", sample_name))
		d_corl['full_path'] = paste0(bin_dir, '/', d_corl$V1)
		# excluding self
		d_corl = subset(d_corl, V1!=sample_name)
		d = d_corl[complete.cases(d_corl),] # na remove
		dd = d[d[,2]>cor_cut,] # keep samples above cutoff

		mixmdl = normalmixEM(dd[,2], k=3, maxrestarts=1000)
		comps = sapply(1:nrow(mixmdl$posterior), function(x) which.max(mixmdl$posterior[x,]))
		best_comp =  which.max(mixmdl$mu)
		ref_files = dd[comps==best_comp,]

		# for highest correlated samples - excluding self
		ref_files = ref_files[order(ref_files[,2], decreasing=T),]
		ref_filenames  = as.character(ref_files[1:(batch_size), 'full_path'])
		ref_filenames = ref_filenames[complete.cases(ref_filenames)]

		# for gender matched highest crrleted reference - excluding self
		gens = read.table(gender_file, header=T)
		gens = gens[complete.cases(gens),] # Note that this is going to remove the low covergae data which did not get gender classification
		this_gender = gens[gens[,1]==sample_name,]
		matched_ref_files = dd[dd$V1 %in% gens[gens$gender==this_gender$gender,1],]
		matched_ref_files = matched_ref_files[order(matched_ref_files[,2], decreasing=T),]
		matched_ref_file_names  = as.character(matched_ref_files[1:(batch_size), 'full_path'])
		matched_ref_file_names = matched_ref_file_names[complete.cases(matched_ref_file_names)]

		# for gender mismatched highest crrleted reference - excluding self
		mismatched_ref_files = dd[dd$V1 %in% gens[gens$gender!=this_gender$gender,1],]
		mismatched_ref_files = mismatched_ref_files[order(mismatched_ref_files[,2], decreasing=T),]
		mismatched_ref_files_names  = as.character(mismatched_ref_files[1:(batch_size), 'full_path'])
		mismatched_ref_files_names = mismatched_ref_files_names[complete.cases(mismatched_ref_files_names)]

		mean_cov = mean_norm(ref_filenames, index)
		median_cov = med_norm(ref_filenames, index, index_file)
		matched_mean_cov = mean_norm(matched_ref_file_names, index)
		matched_med_cov = med_norm(matched_ref_file_names, index, index_file)
		mismatched_mean_cov = mean_norm(mismatched_ref_files_names, index)
		mismatched_med_cov = med_norm(mismatched_ref_files_names, index, index_file)

	return(list("ref_mean_cov"=mean_cov, "ref_med_cov"=median_cov, "ref_n"=length(ref_filenames)
							, "matched_mean_cov"=matched_mean_cov, "matched_med_cov"=matched_med_cov, "matched_n"=length(matched_ref_file_names)
							, "mismatched_mean_cov"=mismatched_mean_cov, "mismatched_med_cov"=mismatched_med_cov, "mismatched_n"=length(mismatched_ref_files_names)))
	}

	# bin_dir was called log2_path before
	output_file = paste0(outputdir, "/", sample_name)
	sample_files = dir(bin_dir, full.names=TRUE)
	sample_file = paste0(bin_dir, '/', sample_name)

	index = read.table(index_file)
	data = get_sample(sample_file, index)
	ref = get_reference(sample_name, sample_file, index, index_file, correlation_path, bin_dir, gender_file, batch_size)

	dat = data.frame(index, data, ref$ref_mean_cov, ref$ref_med_cov, ref$matched_mean_cov, ref$matched_med_cov, ref$mismatched_mean_cov, ref$mismatched_med_cov)
	write.table(dat, file=output_file, sep="\t", row.names=F, quote=F, col.names=F)

}

processLogRtoBin <- function(logr_dir, rbin_dir, sample_name) {
	## input
	# V1-V4: chrom, start, end, count
	# V5, V6: mix_mean, mix_median
	# V7, V8: gender_match_mean, gender_match_median
	# V9, V10: gender_mismatch_mean, gender_mismatch_median
	data = read.table(paste0(logr_dir, '/', sample_name), sep="\t")
	## tmp & output
	tempfile = paste0(rbin_dir, "/", sample_name, ".tmp")
	binfile = paste0(rbin_dir, "/", sample_name)
	## TODO: These should be redunant
	# chr = as.character(data[,1]) 
	# chr[chr=="X"] = "23"
	# chr[chr=="Y"] = "24"
	# data[,1] = as.numeric(chr)
	# data = data[order(data[,1], data[,2], data[,3]),]
	l = log2((data[,4]+1)/(data[,6]+1)) # log2(count/mix_median)
	l1= l-median(l[data[,1]<23])
	l = log2((data[,4]+1)/(data[,8]+1)) # log2(count/gender_match_median)
	l2= l-median(l[data[,1]<23])
	l = log2((data[,4]+1)/(data[,10]+1)) # log2(count/gender_mismatch_median)
	l3 = l-median(l[data[,1]<23])
	data = data.frame(data[,1:3], l1, l2, l3)
	write.table(data, file=tempfile, sep="\t", row.names=F, col.names=F, quote=F)
	rname = RbinConvert_exome_ratio(tempfile, binfile)
}

run_hmm_rbin <- function(rbin_path, sample_name, indexfile, cov_file, cor_path, gender_file, outpath,	batch_size=1000, cov_cut = 20, cor_cut = 0.9) {
	hmm_call <- function(f, active=FALSE) {
			process = data.frame(as.vector(sapply(1:ncol(f), function(x) rep(x, nrow(f))))
														,	as.vector(sapply(1:ncol(f), function(x) 1:nrow(f)))
														, as.vector(f)
														)
			v = ViteRbi(process, active=FALSE)
			colnames(v) = c("sample", "index", "log2", "state")
			v$state = as.factor(v$state)
			sams = unique(v$sample)
			if(active==TRUE) {
				for(x in 1:length(sams)) {
					vv = v[v[,1]==sams[x],]
					g1 = ggplot(data=vv, aes(x=state, y= log2, fill= state)) + geom_violin()+theme_bw()
					g2 = ggplot(data=vv, aes(x=index, y= log2, color= state)) + geom_point(size=0.1)+ylim(-2,2)+theme_bw()
					multiplot(g1, g2, cols=1)
				}
			}
		return(v)
	}
	hmm_call_all_together <- function(f, active=FALSE) {
			sample_index = as.vector(sapply(1:ncol(f), function(x) rep(x, nrow(f))))
			process = data.frame(1, 1, as.vector(f))
			v = ViteRbi(process, active=FALSE)
			colnames(v) = c("sample", "index", "log2", "state")
			v$sample = sample_index
			v$state = as.factor(v$state)
			sams = unique(v$sample)
			if(active==TRUE) {
				for(x in 1:length(sams)) {
					vv = v[v[,1]==sams[x],]
					g1 = ggplot(data=vv, aes(x=state, y= log2, fill= state)) + geom_violin()+theme_bw()
					g2 = ggplot(data=vv, aes(x=index, y= log2, color= state)) + geom_point(size=0.1)+ylim(-2,2)+theme_bw()
					multiplot(g1, g2, cols=1)
				}
			}
		return(v)
	}
	get_reference_filenames_by_type <- function(sample_name, rbin_path, batch_size = 1000, cor_cut = 0.9, type="mixed", cor_path, gender_file) {
		d_corl = read.table(paste0(cor_path, "/", sample_name))
		# d_corl[,1] = basename(as.character(d_corl[,1]))
		d = d_corl[complete.cases(d_corl),]
		dd = d[d[,2]>cor_cut,]
		mixmdl = normalmixEM(dd[,2], k=3, maxrestarts=1000)
		comps = sapply(1:nrow(mixmdl$posterior), function(x) which.max(mixmdl$posterior[x,]))
		best_comp =  which.max(mixmdl$mu)
		ref_files = dd[comps==best_comp,]
		ref_filenames = NULL
		# for random index
		if(type == "mixed") {
			file_index = unique(round(runif(batch_size, 1, nrow(ref_files))))
			ref_filenames = as.character(ref_files[file_index,1])
			# for highest correlated samples - excluding self
			ref_files = ref_files[order(ref_files[,2], decreasing=T),]
			ref_filenames  = as.character(ref_files[2:(batch_size+1),1])
			ref_filenames = ref_filenames[ref_filenames!=sample_name]
		} else {
			gens = read.table(gender_file, header=T)
			gens = gens[complete.cases(gens),] # Note that this is going to remove the low covergae data which did not get gender classification
			this_gender = gens[gens[,1]==sample_name,]
			# for gender matched highest crrleted reference - excluding self
			if(type == "matched") {
				matched_ref_files = dd[dd[,1]%in%gens[gens$gender==this_gender$gender,1],]
				matched_ref_files = matched_ref_files[order(matched_ref_files[,2], decreasing=T),]
				matched_ref_file_names  = as.character(matched_ref_files[2:(batch_size+1),1])
				matched_ref_file_names = matched_ref_file_names[matched_ref_file_names!=sample_name]
				ref_filenames = matched_ref_file_names
			}
			# for gender mismatched highest crrleted reference - excluding self
			if(type == "mismatched") {
				mismatched_ref_files = dd[dd[,1]%in%gens[gens$gender!=this_gender$gender,1],]
				mismatched_ref_files = mismatched_ref_files[order(mismatched_ref_files[,2], decreasing=T),]
				mismatched_ref_files_names  = as.character(mismatched_ref_files[1:(batch_size),1])
				mismatched_ref_files_names = mismatched_ref_files_names[mismatched_ref_files_names!=sample_name]
				ref_filenames = mismatched_ref_files_names
			}
		}
	return(paste(rbin_path, "/", ref_filenames, sep=""))
	}
	get_all_data_by_type <- function(indexfile, binfiles, type="mixed") {
		index = read.table(indexfile)
		chrs = unique(index[,1])
		return(do.call(rbind, sapply(chrs, function(x) RbinRead_exome_ratio(x, min(index[index[,1]==x,2]), max(index[index[,1]==x,3]), type, indexfile, binfiles))))
	}
	process_and_collect <- function(m) {
			pos = m[,1:3]
			m = as.matrix(m[,-(1:3)])
			### NOTE: this is the HMM call 
			m_v = hmm_call_all_together(m, active=FALSE)
			state_calls = m_v[1:nrow(pos),]
			#state_calls = m_v[m_v[,1]==1,]
			state_calls[,1:2] = pos[,1:2]
			calls = extract_calls(state_calls)
		return(list("states"=state_calls, "calls"=calls))
	}
	hmm_by_reference_types <- function(sample_name, rbin_path, cov_file, indexfile, batch_size=1000, cor_cut = 0.9, cov_cut = 20, cor_path, gender_file, states_outname, call_outname) {
		# the mean coverage of each bait location across all samples
		cov = read.table(cov_file)

		# get mismatched filenames for training
		binfiles = get_reference_filenames_by_type(sample_name, rbin_path, batch_size, cor_cut, "mismatched", cor_path, gender_file)
		all_data = get_all_data_by_type(indexfile, binfiles, "mismatched")

		sample_data = get_all_data_by_type(indexfile, paste0(rbin_path, "/", sample_name), "mixed")
		all_data = data.frame(sample_data, all_data[,-(1:3)])

		# run hmm trained on all data that made up each reference separatly
		state_calls = process_and_collect(all_data[cov[,4]>cov_cut,])
		write.table(state_calls$states, file=states_outname, sep="\t", row.names=F, quote=F)
		write.table(state_calls$calls, file=call_outname, sep="\t", row.names=F, quote=F)
	}

	require(ViteRbi)
	print (sample_name)
	system(paste("mkdir ", outpath, "/", sample_name, sep=""))
	states_output_file = paste0(outpath, "/", sample_name, "/", sample_name, "_mixed_states.txt", sep="")
	calls_output_file = paste0(outpath, "/", sample_name, "/", sample_name, "_mixed_calls.txt", sep="")
	hmm_by_reference_types(sample_name, rbin_path, cov_file, indexfile, batch_size, cor_cut, cov_cut, cor_path, gender_file, states_output_file, calls_output_file)
}


`Rbin_example` <- function(chr=1, start=1, stop=20000000, n=1000) {
	
	# Make the example binary file format
	path = system.file("extdata", package="Rbin");
	binname = paste(path, "/", "bin.dat", sep=""); inname = paste(path, "/", "Example_Input.txt", sep="");
	RbinConvert_aCGH(inname,binname);
	
	# Make the example index file and ensure (just for example) that it is sorted correctly.
	index = read.table(inname); index=index[order(index[,1], index[,2], index[,3]),];
	indexfile = paste(path, "/", "Example_Index.txt", sep="");
	write.table(index[,1:3], indexfile, sep="\t", row.names=F, col.names=F, quote=F);
	
	# Replicate the binary filename by n and execute n sequetial file reads returning a matrix for the given position.
	filenames = rep(binname, n);
	r = RbinRead_aCGH(chr, start, stop, indexfile, filenames);
	
	invisible(r)
}
