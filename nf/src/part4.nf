#!/usr/bin/env nextflow

/*
Part4 is used to run step5 - CNV HMM calls (sample by sample). 
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --project test --design design_file.csv --ref ./ref/
    
    Mandatory arguments:
      --project       [string] Name of the project
      --design        [file] A csv file with sample name, CRAM path and CRAI path
                      A file could look like this:
                      name,cram,crai
                      test,test.cram,test.cram.crai
      --ref           [file] Path for the genome FASTA. Used for CRAM decoding.
    
    Optional arguments:
      --test          [flag] test mode

    """.stripIndent()
}