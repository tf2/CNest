library(Rbin)

args = commandArgs(trailingOnly=TRUE)

print(args[1])

if(args[1] == "processtobin") {
	project_root = args[2]
	sample_id = args[3]
	inputfile= paste(project_root, "/txt/", sample_id, sep="")
	tempfile= paste(project_root, "/tmp/", sample_id, sep="")
	binfile= paste(project_root, "/bin/", sample_id, sep="")
	indexfile = paste(project_root, "/index.txt", sep="")
	processtobin(inputfile, tempfile, binfile, indexfile)
	print ("processtobin")
}

if(args[1] == "processtobin_fast") {
	project_root = args[2]
	sample_id = args[3]
	tempfile= paste(project_root, "/tmp/", sample_id, sep="")
	binfile= paste(project_root, "/bin/", sample_id, sep="")
	processtobin_fast(tempfile, binfile)
	print ("processtobin (fast mode)")
}

if(args[1] == "classify_gender") {
	project_root = args[2]
	index_file = args[3]
	qc_file = args[4]
	gender_file = args[5]
	bin_dir= paste(project_root, "/bin", sep="")
	classify_gender(bin_dir, index_file, qc_file)
	assign_gender(qc_file, gender_file)
}

if(args[1] == "generate_correlation") {
	project_root = args[2]
	sample_index = as.numeric(as.character(args[3]))
	indexfile = args[4]
	bin_dir = paste(project_root, "/bin", sep="")
	cor_dir = paste(project_root, "/cor", sep="")
	if(file.exists(paste(project_root,"/txt", sep=""))) {
		system(paste("rm -r ", paste(project_root,"/txt", sep="")))
	}
	if(file.exists(paste(project_root,"/tmp", sep=""))) {
		system(paste("rm -r ", paste(project_root,"/tmp", sep="")))
	}
	generate_correlation(bin_dir, cor_dir, sample_index, indexfile)
	print ("generate_correlation")
}

if(args[1] == "get_references") {
	project_root = args[2]
	sample_index = as.numeric(as.character(args[3]))
	indexfile = args[4]
	gender_file = args[5]
	batchsize = as.numeric(as.character(args[6]))
	bin_dir = paste(project_root, "/bin", sep="")
	cor_dir = paste(project_root, "/cor", sep="")
	outputdir = paste(project_root, "/logr", sep="")

	get_references(sample_index, indexfile, gender_file, outputdir, cor_dir, bin_dir, batchsize)
	print ("get_references")
}

if(args[1] == "processLogRtoBin") {
	project_root = args[2]
	sample_index = as.numeric(as.character(args[3]))
	logR_dir = paste(project_root, "/logr", sep="")
	rbin = paste(project_root, "/rbin", sep="")
	processLogRtoBin(logR_dir, rbin, sample_index)
	print ("processLogRtoBin")
}

if(args[1] == "run_hmm_rbin") {
	project_root = args[2]
	sample_index = as.numeric(as.character(args[3]))
	indexfile = args[4]
	cov_file = args[5]
	gender_file = args[6]
	outpath = args[7]
	batch_size = args[8]
	rbin = paste(project_root, "/rbin", sep="")
	cor_dir = paste(project_root, "/cor", sep="")
	if(file.exists(paste(project_root,"/logr", sep=""))) {
		system(paste("rm -r ", paste(project_root,"/logr", sep="")))
	}
	print(rbin)
	print(sample_index)
	run_hmm_rbin(rbin, sample_index, indexfile, cov_file, cor_dir, gender_file, outpath, batch_size)
	print ("run_hmm_rbin")
}




