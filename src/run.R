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
	bin_dir = args[2]
	index_tab = args[3]
	qc_file = args[4]
	gender_file = args[5]
	classify_gender(bin_dir, index_tab, qc_file)
	assign_gender(qc_file, gender_file)
}

if(args[1] == "generate_coverage") {
	bin_dir = args[2]
	index_tab = args[3]
	cov_file = args[4]
	generate_coverage(bin_dir, index_tab, cov_file)
}

if(args[1] == "generate_correlation") {
	bin_dir = args[2]  # path to all bin files
	cor_dir = args[3] # output cor dir
	sample_name = args[4] # sample name 
	index_file = args[5] # index_tab.txt path
	generate_correlation(bin_dir, cor_dir, sample_name, index_file)
	print ("generate_correlation")
}

if(args[1] == "get_references") {
	# ['Rscript', '/resources/run.R', 'get_references', bin_dir, cor_dir, logr_dir, sample_name, index_tab, gender_file, batch_size]
	bin_dir = args[2]
	cor_dir = args[3]
	logr_dir = args[4]
	sample_name = args[5]
	index_file = args[6]
	gender_file = args[7]
	batchsize = as.numeric(as.character(args[8]))
	get_references(sample_name, index_file, gender_file, logr_dir, cor_dir, bin_dir, batchsize)
	print ("get_references")
}

if(args[1] == "processLogRtoBin") {
	# cmd6 = ['Rscript', '/resources/run.R', 'processLogRtoBin', logr_dir, rbin_dir, sample_name]
	logr_dir = args[2]
	rbin_dir = args[3]
	sample_name = args[4]
	processLogRtoBin(logr_dir, rbin_dir, sample_name)
	print ("processLogRtoBin")
}

if(args[1] == "run_hmm_rbin") {
	#   Rscript /resources/run.R run_hmm_rbin ${bin_dir} ${cor_dir} ${cnv_dir} ${sample_name} \
    #        ${index_tab} ${cov_file} ${gender_file} ${batch_size}
	rbin_dir = args[2]
	cor_dir = args[3]
	cnv_dir = args[4]
	sample_name = args[5]
	index_file = args[6]
	cov_file = args[7] # cov_file=/resources/txt_covered_all_coverage.txt 
	gender_file = args[8]
	batch_size = args[9]
	print(sample_name)
	run_hmm_rbin(rbin_dir, sample_name, index_file, cov_file, cor_dir, gender_file, cnv_dir, batch_size)
	print ("run_hmm_rbin")
}