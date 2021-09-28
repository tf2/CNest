version 1.0

import "./cnest_part0.wdl" as CnestPart0
import "./cnest_part1.wdl" as CnestPart1
import "./cnest_part2.wdl" as CnestPart2
import "./cnest_part3.wdl" as CnestPart3
import "./cnest_part4.wdl" as CnestPart4

workflow CnestWorkflow {
  input {

    String project
    Array[String] samples
    
    Int preemptible_tries = 3
    
    #part0
    File bedgz

    #part1
    Array[File] bam_file
    Array[File] bai_file
    File ref
  
    Boolean test = false
    
    #part3 and part 4
    Int batch
    Int wgs = 1
    Float cor_cut = 0.9
    Int cov_cut = 20
    Boolean skipem = false
    
    Int part3_addtional_mem_gb = 0
    Int part4_addtional_disk_gb = 0
    Int part4_addtional_mem_gb = 0
  }
  
    call CnestPart0.step0 {
    input:
      test = test,
      bedgz = bedgz
  }

  call CnestPart0.step1 {
    input: 
      bed = step0.ch_bed,
      project = project
  }
  
  scatter (i in range(length(samples))) {
      call CnestPart1.step2 {
        input:
          test = test,
          project = project,
          name = samples[i],
          file_path = bam_file[i],
          file_path_index = bai_file[i],
          ch_ref = ref,
          ch_index_bed = step1.ch_index_bed,
          preemptible_tries = preemptible_tries
      }
  }
  
  Array[File] list_of_bins = select_all(step2.step2_out)
  
  call CnestPart1.tarzip_bins {
    input:
      list_of_bins = list_of_bins,
      project = project
  }

  call CnestPart2.gender_qc {
    input:
      bin_dir = tarzip_bins.bin_zip,
      index_tab = step1.ch_index_tab,
      preemptible_tries = preemptible_tries
  }
   
  scatter (sample_name in samples) {
     call CnestPart3.logR_ratio {
       input:
         project = project,
         bin_dir = tarzip_bins.bin_zip,
         index_tab = step1.ch_index_tab,
         gender = gender_qc.gender_class,
         sample_name = sample_name,
         batch = batch,
         wgs = wgs,
         cor_cut = cor_cut,
         skipem = skipem, 
         preemptible_tries = preemptible_tries,
         addtional_mem = part3_addtional_mem_gb
     }
  }
    
  scatter (sample_name in samples) {
     call CnestPart4.hmm_calls {
       input:
         project = project,
         rbinfiles = logR_ratio.rbin,
         corfiles = logR_ratio.cor,
         cov_file = gender_qc.mean_coverage,
         index = step1.ch_index_tab,
         gender = gender_qc.gender_class,
         batch = batch,
         wgs = wgs,
         sample_name = sample_name,
         cor_cut = cor_cut,
         cov_cut = cov_cut,
         skipem = skipem,
         addtional_disk_gb = part4_addtional_disk_gb,
         preemptible_tries = preemptible_tries,
         addtional_mem = part4_addtional_mem_gb
     }
  }

  output{
    #part0
    File out_ch_index_tab = step1.ch_index_tab
    File out_ch_index = step1.ch_index
    File out_ch_index_bed = step1.ch_index_bed
    #part1
    File out_bin_zip = tarzip_bins.bin_zip
    #part2
    File out_gender_qc = gender_qc.gender_qc
    File out_gender_class = gender_qc.gender_class
    File out_mean_coverage = gender_qc.mean_coverage
    #part3
    Array[File] out_cor = logR_ratio.cor
    Array[File] out_logr = logR_ratio.logr
    Array[File] out_rbin = logR_ratio.rbin
    #part4
    Array[File] mixed_calls = hmm_calls.mixed_calls
    Array[File] mixed_stats = hmm_calls.mixed_stats
  }
}

