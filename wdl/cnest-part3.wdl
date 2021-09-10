version 1.0

workflow CnestWorkflow {
  input {
    
    String project
    File index_tab
    File bin_dir
    File gender
    Int batch
    Array[String] samples
    Int wgs = 1
    Float cor_cut = 0.9
    Boolean skipem = false
    
    Int preemptible_tries = 3
    Int addtional_mem_gb = 0

  }
    
   scatter (sample_name in samples) {
    call logR_ratio {
      input:
        project = project,
        bin_dir = bin_dir,
        index_tab = index_tab,
        gender = gender,
        sample_name = sample_name,
        batch = batch,
        wgs = wgs,
        cor_cut = cor_cut,
        skipem = skipem, 
        preemptible_tries = preemptible_tries,
        addtional_mem = addtional_mem_gb
    }
  }

  output {
    Array[File] out_cor = logR_ratio.cor
    Array[File] out_logr = logR_ratio.logr
    Array[File] out_rbin = logR_ratio.rbin
  }


}

task logR_ratio {
  input {
    String project
    File bin_dir
    File index_tab 
    File gender 
    String sample_name
    Int batch
    Int wgs
    Float cor_cut
    Boolean skipem
    Int addtional_mem
    
    Int preemptible_tries
  }
  String bin_dir_name = basename(bin_dir, ".tar.gz")
  
  Int cal_memory = ceil(5.6 * batch * wgs / 100) + addtional_mem
  Int memory_gb = if cal_memory > 2 then cal_memory else 2
  
  Float bin_dir_size = size(bin_dir, "GiB")
  # Additional 20G for text output
  Int disk_size = ceil(bin_dir_size + 20)

  command {
    tar -xvf ~{bin_dir}
    
    mkdir -p ./~{project}/cor/ ./~{project}/logr/ ./~{project}/rbin/
        
    cnest.py step4 \
      --bindir ./~{bin_dir_name}/ \
      --indextab ~{index_tab} \
      --gender ~{gender} \
      --sample ~{sample_name} \
      --batch ~{batch} \
      --cordir ./~{project}/cor/ \
      --logrdir ./~{project}/logr/ \
      --rbindir ./~{project}/rbin/ \
      ~{true="--skipem" false="" skipem} \
      --cor ~{cor_cut} 
      
    }
    
  output {
    File cor = "./~{project}/cor/~{sample_name}"
    File logr = "./~{project}/logr/~{sample_name}"
    File rbin = "./~{project}/rbin/~{sample_name}"
  }
    
  runtime {
    docker: "tomas81/cnest:dev"
    memory: memory_gb + " GB"
    preemptible: preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
  }
}

