#!/usr/bin/env nextflow

// Re-usable componext for adding a helpful help message in our Nextflow script
def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --project ukb --bed index.bed
    Mandatory arguments:
      --fastq_list                  [file] A comma seperated file with all the fastq files locations
                                    A header is expected, and 3 columns that define the following:
                                    - accession number, identifier for fastq file pair
                                    - link (ftp, https) to the fastq 1
                                    - link (ftp, https) to the fastq 2
                                    
                                    A file could look like this:
                                    accession,fastq_1,fastq_2
                                    ERR908503,ftp://ERR908503_1.fastq.gz,ftp://ERR908503_2.fastq.gz

    """.stripIndent()
}

if (params.index) ch_index = Channel.value(file(params.index))

// Re-usable process skeleton that performs a simple operation, listing files
// nextflow run main.nf --project ukb
process step1 {
  tag "${index}"
  echo true
  publishDir "results", mode: 'copy'


  input: 
  file(index) from ch_index

  output: 
  set file("${params.project}/index.txt"), file("${params.project}/index_tab.txt") into ch_indexs
  file "${params.project}/bin/" into ch_bin
  file "${params.project}/txt/" into ch_txt
  file "${params.project}/tmp/" into ch_tmp
  file ("${params.project}") into ch_project

  script:
  """
  mv ./${index} /input_location/
  cnest.py step1 --project ${params.project} --bed ${index}
  mv /output_location/${params.project} .
  """
}

process step2 {
  tag "${index}"
  echo true
  publishDir "results", mode: 'copy'


  input: 
  file(project) from ch_project

  output: 



  script:
  """
    ls -alL ukb
  """
}

// Step1
// ```bash
// docker run -v "${index_path}:/input_location" -v "${output_path}:/output_location" -it --rm smshuai/cnest:dev step1 --project ukbb_wes --index index.txt
// ```

// Step2
// ```bash
// # BAM
// docker run -v "${input_path}:/input_location" -v "${output_path}/output_location" -it --rm smshuai/cnest:dev step1 -project ukbb_wes --sample 'A' --input 'a.bam'

// # CRAM (Need to mount ref path)
// docker run -v "${input_path}:/input_location" -v "${output_path}/output_location" -v "${ref_path}:/ref" -it --rm smshuai/cnest:dev step1 -project ukbb_wes --sample 'A' --input 'a.bam'
// ```