#!/usr/bin/env bash
ref_path=$1 # Change this to your genome reference path
ref_path=`realpath $ref_path`
echo "REF_PATH=$ref_path"
executor=lsf # set this to slurm for SLURM HPCs

if [[ -d results ]]
then
  rm -r ./results/
fi

if [[ -d work ]]
then
  rm -r ./work/
fi

nextflow run ../main.nf --part 0 --project test_proj --bed ./test.bed -profile $executor 

nextflow run ../main.nf --part 1 --project test_proj --design design.csv --ref $ref_path --indexb ./results/test_proj/index.bed

# nextflow run ../main.nf --part 2 --bindir ./results/test_proj/bin --index ./results/test_proj/index_tab.txt

# nextflow run ../main.nf --part 3 --project test_proj --bindir ./results/test_proj/bin --index ./results/test_proj/index_tab.txt --gender ./results/gender_classification.txt
