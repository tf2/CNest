#!/usr/bin/env nextflow

/*
Part2 is used to run step3 - gender classification (all samples together). QC is required after part2.
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run part2.nf --binDir ./results/test_proj/bin --index ./results/test_proj/index_tab.txt
    
    Mandatory arguments:
      --binDir       [path] Path to the bin file folder
      --index        [file] Path to index_tab.txt
    """.stripIndent()
}

if (params.binDir) ch_bin = Channel.value(file(params.binDir))
if (params.index) ch_index = Channel.value(file(params.index))

process step3 {
  echo true
  publishDir "results/", mode: "copy"

  input:
  path bin_dir from ch_bin
  path index from ch_index

  output:
  path "gender_qc.txt" into ch_gender_qc
  path "gender_classification.txt" into ch_gender_file
  path "mean_coverage.txt" into ch_cov_file

  script:
  """
    cnest.py step3 \
    --indextab $index \
    --bindir $bin_dir \
    --qc gender_qc.txt \
    --gender gender_classification.txt \
    --cov mean_coverage.txt
  """
}