#!/usr/bin/env nextflow

/*
Part3 is used to run step4 - from bin to logR to logRbin (sample by sample). 
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run part3.nf --project test_proj --binDir ./results/test_proj/bin --index ./results/test_proj/index_tab.txt --gender ./results/gender_classification.txt
    
    Mandatory arguments:
      --project      [string] Name of the project
      --bindir       [path] Path to the bin file folder
      --index        [file] Path to index_tab.txt
      --gender       [file] Path to gender_classification.txt
      --batch        [int] Batch size used in reference search
    
    Optional arguments:
      --test         [flag] test mode (only 3 samples are used)

    """.stripIndent()
}

if (params.bindir) {
  ch_bin = Channel.value(file(params.bindir))
  all_samples = file(params.bindir).list()
  ch_sample_names = Channel.from(all_samples)
}
if (params.test) ch_sample_names = ch_sample_names.take(3)
if (params.index) ch_index = Channel.value(file(params.index))
if (params.gender) ch_gender = Channel.value(file(params.gender))

process step4 {
  tag "${sample_name}"
  echo true
  publishDir "results/", mode: "copy"
  memory { 2.GB * params.batch / 100 }
  time { 20.m * params.batch / 100  }

  input:
  path bin_dir from ch_bin
  path index from ch_index
  path gender from ch_gender
  val sample_name from ch_sample_names

  output:
  path "${params.project}/cor/$sample_name" into ch_cor
  path "${params.project}/logr/$sample_name" into ch_logr
  path "${params.project}/rbin/$sample_name" into ch_rbin

  script:
  """
    echo "Processing sample $sample_name"
    mkdir -p ${params.project}/cor ${params.project}/logr ${params.project}/rbin
    cnest.py step4 \
      --bindir $bin_dir \
      --indextab $index \
      --gender $gender \
      --sample $sample_name \
      --batch ${params.batch} \
      --cordir ${params.project}/cor \
      --logrdir ${params.project}/logr \
      --rbindir ${params.project}/rbin
  """
}
