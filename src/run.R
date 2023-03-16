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

if(args[1] == "generate_correlation_batch") {
	bin_dir = args[2]  # path to all bin files
	cor_dir = args[3] # output cor dir
	batch_size = as.numeric(as.character(args[4])) # index_tab.txt path
	start_pos = as.numeric(as.character(args[5])) # index_tab.txt path
	index_file = args[6] # index_tab.txt path
	generate_correlation_batch(bin_dir, cor_dir, batch_size, start_pos, index_file)
	print ("generate_correlation batch")
}

if(args[1] == "generate_correlation_chunk") {
	bin_dir = args[2]  # path to all bin files
	cor_dir = args[3] # output cor dir
	target_size = as.numeric(as.character(args[4])) # size of target samples
	start_pos = as.numeric(as.character(args[5])) # starting position in file list - should be chunked rel to batch_size
	index_file = args[6] # index_tab.txt path
	generate_correlation_chunk(bin_dir, cor_dir, target_size, start_pos, index_file)
	print ("generate_correlation batch")
}

if(args[1] == "generate_correlation_chunk_batch") {
	bin_dir = args[2]  # path to all bin files
	cor_dir = args[3] # output cor dir
	target_size = as.numeric(as.character(args[4])) # size of target samples
	batch_size = as.numeric(as.character(args[5])) # size of target samples
	start_pos = as.numeric(as.character(args[6])) # starting position in file list - should be chunked rel to batch_size
	index_file = args[7] # index_tab.txt path
	generate_correlation_chunk_batch(bin_dir, cor_dir, target_size, start_pos, index_file)
	print ("generate_correlation_chunk_batch")
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
    cor_cut = as.numeric(as.character(args[9]))
    skip_em = as.character(args[10])
    if(skip_em=="True") {
        skip_em = TRUE
    } else {
        skip_em = FALSE
    }
	get_references(sample_name, index_file, gender_file, logr_dir, cor_dir, bin_dir, batchsize, cor_cut, skip_em)
	print ("get_references")
}


if(args[1] == "get_references_to_rbin") {
	#index_file, gender_file, rbin_dir,cor_dir, bin_dir, batch_size = 1000, target_size, start_pos, cor_cut = 0.9, skip_em=FALSE
	index_file = args[2]
	gender_file = args[3]
	rbin_dir = args[4]
	cor_dir = args[5]
	bin_dir = args[6]
	batchsize = as.numeric(as.character(args[7]))
	target_size = as.numeric(as.character(args[8]))
	start_pos = as.numeric(as.character(args[9]))
    cor_cut = as.numeric(as.character(args[10]))
    skip_em = as.character(args[11])
    if(skip_em=="True") {
        skip_em = TRUE
    } else {
        skip_em = FALSE
    }
	batch_get_references_to_rbin(index_file, gender_file, rbin_dir,
                                    	cor_dir, bin_dir, batchsize, 
                                        target_size, start_pos, cor_cut, skip_em)
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
    cov_cut = as.numeric(as.character(args[10]))
    cor_cut = as.numeric(as.character(args[11]))
    skip_em = as.character(args[12])
    if(skip_em=="True") {
        skip_em = TRUE
    } else {
        skip_em = FALSE
    }
	print(sample_name)
	run_hmm_rbin(rbin_dir, sample_name, index_file, cov_file, cor_dir, gender_file, cnv_dir, batch_size, cov_cut, cor_cut, skip_em)
	print ("run_hmm_rbin")
}
