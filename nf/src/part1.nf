#!/usr/bin/env nextflow

/*
* Part1 is used to run step1 & 2 - index generation & bam/cram2bin.
*/


// Re-usable componext for adding a helpful help message in our Nextflow script
def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run part1.nf --project test --design design_file.csv --ref ref.fa [--bed bait.bed]
    
    Mandatory arguments:
      --project       [string] Name of the project
      --design        [file] A csv file with sample name, CRAM path and CRAI path
                      A file could look like this:
                      name,cram,crai
                      test,test.cram,test.cram.crai
      --ref           [file] Path for the genome FASTA. Used for CRAM decoding.
    
    Optional arguments:
      --bed           [file] Path for the bait file (in BED3 format; hg38 WGS baits will be used by default)
      --test          [flag] test mode (only 10 samples are used)
    """.stripIndent()
}

if (params.bed) ch_bed = Channel.value(file(params.bed))
if (params.ref) ch_ref = Channel.value(file(params.ref))
if (params.design) {
  Channel.fromPath(params.design)
    .splitCsv(sep: ',', skip: 1)
    .map { name, file_path, index_path -> [ name, file(file_path), file(index_path) ] }
    .set { ch_files_sets }
}
// In test mode, a max of 10 samples are used
if (params.test) {
  ch_files_sets = ch_files_sets.take(10)
}

ch_bedgz = Channel.value(file("$baseDir/data/hg38.1kb.baits.bed.gz"))
process step0 {
  tag "${params.project}"
  echo true

  input:
  file(bedgz) from ch_bedgz

  output:
  file("hg38.1kb.baits.bed") into ch_bed

  when:
  !params.bed

  script:
  if (params.test)
    """
    gzip -cd ${bedgz} | head -1000 > "hg38.1kb.baits.bed"
    """
  else
    """
    gzip -cd ${bedgz} > "hg38.1kb.baits.bed"
    """
}

// Step1 create work directory
process step1 {
  tag "${params.project}"
  publishDir "results/", mode: "copy", pattern: "${params.project}/index*"
  echo true

  input: 
  file(bed) from ch_bed

  output: 
  file ("${params.project}") into ch_proj
  path "${params.project}/index_tab.txt" into ch_index_tab
  path "${params.project}/index.txt" into ch_index
  path "${params.project}/index.bed" into ch_index_bed

  script:
  """
  cnest.py step1 --project ${params.project} --bed ${bed}
  """
}

process step2 {
  tag "id:${name}-file:${file_path}-index:${index_path}"
  publishDir "results/", mode: "copy"
  errorStrategy 'ignore'
  echo true

  input:
  set val(name), file(file_path), file(index_path) from ch_files_sets
  file("genome.fa") from ch_ref
  file(project) from ch_proj

  output:
  path "${params.project}/bin/$name" into ch_bin

  script:
  if (params.test)
    """
    cnest.py --debug step2 --project ${params.project} --sample ${name} --input ${file_path} --fasta genome.fa --fast
    """
  else
    """
    cnest.py step2 --project ${params.project} --sample ${name} --input ${file_path} --fasta genome.fa --fast
    """
}