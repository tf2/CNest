version 1.0

workflow CnestWorkflow {
  input {
    
    String project
    File index
    Array[File] rbinfiles
    Array[File] corfiles
    File cov
    File gender
    Int batch
    Array[String] samples
    Float cor_cut = 0.9
    Int cov_cut = 20
    Boolean skipem = false

    Int wgs = 1
    Int addtional_disk_gb = 0
    Int addtional_mem_gb = 0
    Int preemptible_tries = 3

  }
    
  scatter (sample_name in samples) {
    call hmm_calls {
      input:
        project = project,
        rbinfiles = rbinfiles,
        corfiles = corfiles,
        cov_file = cov,
        index = index,
        gender = gender,
        batch = batch,
        wgs = wgs,
        sample_name = sample_name,
        cor_cut = cor_cut,
        cov_cut = cov_cut,
        skipem = skipem,
        addtional_disk_gb = addtional_disk_gb,
        preemptible_tries = preemptible_tries,
        addtional_mem = addtional_mem_gb
    }
  }

  output {
    Array[File] mixed_calls = hmm_calls.mixed_calls
    Array[File] mixed_stats = hmm_calls.mixed_stats
  }


}

task hmm_calls {
  input {
    String project
    Array[File] rbinfiles
    Array[File] corfiles
    File index 
    File gender
    File cov_file
    Int batch
    String sample_name 
    Int wgs
    Float cor_cut
    Int cov_cut
    Boolean skipem
    
    Int addtional_disk_gb
    Int preemptible_tries
    Int addtional_mem
  }
  Int cal_memory = ceil(5.6 * batch * wgs / 100) + addtional_mem
  Int memory_gb = if cal_memory > 2 then cal_memory else 2

  Int disk_size = 50 + addtional_disk_gb

  command {
    echo "Processing sample $sample_name"
    mkdir -p ~{project}/cnv/ rbindir/ cordir/
    
    echo "~{sep='\n' rbinfiles}" > rbinfiles.txt
    echo "~{sep='\n' corfiles}" > corfiles.txt
    
    for file in $(cat rbinfiles.txt); do mv "$file" rbindir; done
    for file in $(cat corfiles.txt); do mv "$file" cordir; done
    
    cnest.py step5 \
      --indextab ~{index} \
      --rbindir rbindir/ \
      --cordir cordir/ \
      --cnvdir ~{project}/cnv/ \
      --cov    ~{cov_file} \
      --sample ~{sample_name} \
      --gender ~{gender} \
      --batch ~{batch} \
      ~{true="--skipem" false="" skipem} \
      --cor ~{cor_cut} \
      --covc ~{cov_cut} 
    }
    
  output {
    File mixed_calls = "~{project}/cnv/${sample_name}/${sample_name}_mixed_calls.txt"
    File mixed_stats = "~{project}/cnv/${sample_name}/${sample_name}_mixed_states.txt"
  }
    
  runtime {
    docker: "tomas81/cnest:dev"
    memory: memory_gb + " GB"
    preemptible: preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
  }
}

