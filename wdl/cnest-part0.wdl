version 1.0

workflow CnestWorkflow {
  input {
    
    String project
    File bedgz
  
    Boolean test = false
  }
    

  call step0 {
    input:
      test = test,
      bedgz = bedgz
  }

  call step1 {
    input: 
      bed = step0.ch_bed,
      project = project
  }

  output{
    File out_ch_index_tab = step1.ch_index_tab
    File out_ch_index = step1.ch_index
    File out_ch_index_bed = step1.ch_index_bed
  }


}

task step0 {
  input {
    Boolean test
    File bedgz
    
  }
  String bed_name = basename(bedgz, ".gz")


  command {
    if [[ "~{test}" = "true" ]]
    then
      gzip -cd ~{bedgz} | head -1000 > ~{bed_name}
    else
      gzip -cd ~{bedgz} > ~{bed_name}
    fi
  }
    
  output{
    File ch_bed = "~{bed_name}"
  }
    
  runtime{
    docker: "quay.io/smshuai/cnest:dev2"
  }
}


task step1 {
  input {
    File bed
    String project
  }

  command {

    cnest.py step1 --project ~{project} --bed ~{bed}

  }
    
  output{
    File ch_index_tab = "~{project}/index_tab.txt"
    File ch_index = "~{project}/index.txt" 
    File ch_index_bed = "~{project}/index.bed"
  }
    
  runtime{
    docker: "quay.io/smshuai/cnest:dev2"
  }
}


