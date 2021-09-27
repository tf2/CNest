version 1.0

workflow CnestWorkflow {
  input {
    String project
    Array[String] sample_name
    Array[File] bam_file
    Array[File] bai_file
    File ref
    File indexb
  
    Boolean test = false
    
    Int preemptible_tries = 3
  }
    
# Calls 

  # The original nf workflow had a lines that would take in a CSV file of samples 
  # to be used in step 2 task, but since we can use Terra's Data tables as
  # the CSV and approriately assign each shard with the bam and sample name
  # this task has been commented out. It might be useful later if the
  # workflow needs to be run outside of Tera and needs to use CSVs again. 
  
  #call csv2tsv {
  #  input:
  #    csv_file=design
  #}

  #Array[Array[String]] ch_files_sets = read_tsv(csv2tsv.design_tsv)
  
  #scatter (sample in ch_files_sets) {
  #  if (sample[0] != "name") { #Used to skip header
  #    call step2 {
  #      input:
  #        test = test,
  #        project = project,
  #        name = sample[0],
  #        file_path = sample[1],
  #        file_path_index = sample[2],
  #        ch_ref = ref,
  #        ch_index_bed = indexb
  #    }
  #  }
  #}
  
  scatter (i in range(length(sample_name))) {
      call step2 {
        input:
          test = test,
          project = project,
          name = sample_name[i],
          file_path = bam_file[i],
          file_path_index = bai_file[i],
          ch_ref = ref,
          ch_index_bed = indexb,
          preemptible_tries = preemptible_tries
      }
  }
  
  Array[File] list_of_bins = select_all(step2.step2_out)
  
  call tarzip_bins {
    input:
      list_of_bins = list_of_bins,
      project = project
  }

  output {
    #Array[File?] out_step2out = step2.step2_out
    File out_bin_zip = tarzip_bins.bin_zip
  }

}
# Tasks

task step2 {
  input {
    String project
    String name
    Boolean test
    File file_path
    File file_path_index
    File ch_ref
    File ch_index_bed
    
    Int preemptible_tries
  }

  Float file_path_size = size(file_path, "GB")
  Float ch_ref_size = size(ch_ref, "GB")
  # Sometimes the output is larger than the input, or a task can spill to disk.
  # In these cases we need to account for the input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
  Float disk_multiplier = 2.5
  Int disk_size = ceil(file_path_size + ch_ref_size + (disk_multiplier * file_path_size) + 20)

  command {
    set -euo pipefail
    
    mkdir -p ~{project}/tmp/ ~{project}/bin/
    cp  ~{ch_index_bed} ./~{project}/index.bed
    
    export INDEX_DIR=$(readlink -f ~{file_path_index} | xargs dirname)
    mv ~{file_path_index} $INDEX_DIR/~{basename(file_path)}.bai
    
    if [[ "~{test}" = "true" ]]
    then
      cnest.py --debug step2 --project ~{project} --sample ~{name} --input ~{file_path} --fasta ~{ch_ref} --fast
    else
      cnest.py step2 --project ~{project} --sample ~{name} --input ~{file_path} --fasta ~{ch_ref} --fast
    fi
  }
    
  output {
    File step2_out = "~{project}/bin/~{name}"
  }
    
  runtime {
    docker: "tomas81/cnest:dev"
    preemptible: preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
  }
}

task tarzip_bins {
  input {
    Array[File] list_of_bins
    String project
    
  }
  command {
    mkdir ~{project}_bindir
    echo "~{sep='\n' list_of_bins}" > bin.txt
    
    for file in $(cat bin.txt); do mv "$file" ./~{project}_bindir; done
    
    tar -czvf ~{project}_bindir.tar.gz ~{project}_bindir
    ls
  }
  output {
    File bin_zip = "~{project}_bindir.tar.gz"
  }
    
  runtime {
    docker: "tomas81/cnest:dev"
  }
}

#task csv2tsv {
#  input {
#    File csv_file
#  }
#  command {
#  
#    set -euo pipefail
#    
#    python <<CODE
#import csv
#csv.writer(open('design.tsv', 'w+'), delimiter='\t').writerows(csv.reader(open("~{csv_file}"))) 
#CODE
#  }
#    output {
#    File design_tsv = "design.tsv"
#  }
#    
#  runtime {
#    docker: "python:3.8"
#  }
#}

