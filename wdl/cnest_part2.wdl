version 1.0

workflow CnestWorkflow {
  input {
    
    String project
    File bin_dir
    File index_tab
    
    Int preemptible_tries = 3

  }
    

  call gender_qc {
    input:
      bin_dir = bin_dir,
      index_tab = index_tab,
      preemptible_tries = preemptible_tries
  }


  output {
    File out_gender_qc = gender_qc.gender_qc
    File out_gender_class = gender_qc.gender_class
    File out_mean_coverage = gender_qc.mean_coverage
  }


}

task gender_qc {
  input {
    File bin_dir  
    File index_tab
    Int preemptible_tries
  }
  String bin_dir_name = basename(bin_dir, ".tar.gz")
  
  Float bin_dir_size = size(bin_dir, "GiB")
  # Additional 20G for text output
  Int disk_size = ceil(bin_dir_size + 20)

  command {
    set -euo pipefail
    
    tar -xvf ~{bin_dir}
    
    cnest.py step3 \
      --indextab ~{index_tab} \
      --bindir ./~{bin_dir_name}/ \
      --qc gender_qc.txt \
      --gender gender_classification.txt \
      --cov mean_coverage.txt
    }
    
  output {
    File gender_qc = "gender_qc.txt"
    File gender_class = "gender_classification.txt"
    File mean_coverage = "mean_coverage.txt"
  }
    
  runtime {
    docker: "quay.io/smshuai/cnest:dev2"
    preemptible: preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
  }
}

