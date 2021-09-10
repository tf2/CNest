# Quick start with `CNest` 

To run `CNest`, either singularity or docker can be used.
Below you can find examples to run each step with singularity or docker. 

Alternatively, if you have many samples, you can simply use our [Nextflow pipeline](https://github.com/smshuai/CNest-nf) to run across all samples on a HPC environment with job scheduler like SLURM (`sbatch`) or LSF (`bsub`) etc., or a cloud environment supporting Nextflow.

## Using Singularity
### Getting help
```bash
# Tested on singularity version 3.5.0 and 3.7.0-1.el7
singularity pull docker://tomas81/cnest:dev
# Show help of available commands
singularity run docker://tomas81/cnest:dev -h
# Show help for a given step
singularity run docker://tomas81/cnest:dev step2 -h
```

### Run Step1
```bash
singularity run -B "${index_path}:/input,${output_path}:/output" --pwd "/output/" \
    docker://tomas81/cnest:dev step1 \
    --project $project_name \
    --bed "/input/$index_bed"
```

### Run Step2
```bash
# BAM
singularity run -B "${input_path}:/input,${output_path}:/output" --pwd /output/ \
    docker://tomas81/cnest:dev step2 \
    --project 'test_proj' \
    --sample 'test_bam' \
    --input '/input/test.bam'

# CRAM (Need to mount ref path)
singularity run -B "${input_path}:/input,${output_path}:/output,${ref_path}:/ref" --pwd /output/ \
    docker://tomas81/cnest:dev step2 \
    --project $project_name \
    --sample $sample_name \
    --input "/input/$cram_name"

# CRAM (fast mode)
singularity run -B "${input_path}:/input,${output_path}:/output,${ref_path}:/ref" --pwd /output/ \
    docker://tomas81/cnest:dev step2 \
    --project $project_name \
    --sample $sample_name \
    --input "/input/$cram_name" \
    --fasta "/ref/$fasta" \
    --fast
```

### Run Step3
```bash
singularity run -B "${output_path}:/output" --pwd "/output/" \
    docker://tomas81/cnest:dev step3 \
    --project $project_name
```

### Run Step4
```bash
index_tab=$project_name/index_tab.txt
bin_dir=$project_name/bin
cor_dir=$project_name/cor
logr_dir=$project_name/logr
rbin_dir=$project_name/rbin
gender_file=$project_name/gender_classification.txt
singularity run -B "${output_path}:/output" --pwd "/output/" \
    docker://tomas81/cnest:dev step4 \
    --indextab $index_tab \
    --bindir $bin_dir \
    --cordir $cor_dir \
    --logrdir $logr_dir \
    --rbindir $rbin_dir \
    --sample $sample_name \
    --gender $gender_file \
    --batch 1000
```

### Run Step5
```bash
index_tab=$project_name/index_tab.txt
bin_dir=$project_name/bin
cor_dir=$project_name/cor
cnv_dir=$project_name/cnv
rbin_dir=$project_name/rbin
gender_file=$project_name/gender_classification.txt
cov_file=$project_name/mean_coverage.txt
singularity run -B "${output_path}:/output" --pwd "/output/" \
    docker://tomas81/cnest:dev step5 \
    --indextab $index_tab \
    --rbindir $rbin_dir \
    --cordir $cor_dir \
    --cnvdir $cnv_dir \
    --cov    $cov_file \
    --sample $sample_name \
    --gender $gender_file \
    --batch 100
```

## Using Docker
### Getting help
```bash
# Show help of available commands
docker run -it --rm tomas81/cnest:dev -h
# Show help for a given step
docker run -it --rm tomas81/cnest:dev step2 -h
```

### Run Step1
Note: `index.bed` must have the same chromosome names as the BAM/CRAM file.
```bash
docker run -v "${index_path}:/input" -v "${output_path}:/output" -w "/output" -it --rm \
    tomas81/cnest:dev step1 \
    --project test_proj \
    --bed /wkdir/index.bed
```

### Run Step2
```bash
# BAM
docker run -v "${input_path}:/input" -v "${output_path}:/output" -w "/output" -it --rm \
    tomas81/cnest:dev step2 \
    --project test_proj \
    --sample 'test_bam' \
    --input '/input/test.bam'

# CRAM (Need to mount ref path)
docker run -v "${input_path}:/input" -v "${output_path}:/output" -v "${ref_path}:/ref" -w "/output" -it --rm \
    tomas81/cnest:dev step2 \
    --project test_proj \
    --sample 'test_cram' \
    --input '/input/test.cram'

# CRAM fast mode
docker run -v "${input_path}:/input" -v "${output_path}:/output" -v "${ref_path}:/ref" -w "/output" -it --rm \
    tomas81/cnest:dev step2 \
    --project test_proj \
    --sample 'test_cram' \
    --input '/input/test.cram' \
    --fast \
    --fasta '/ref/GCA_000001405.15_GRCh38_full_analysis_set.fna'
```

### Run Step3
```bash
index_tab=$project_name/index_tab.txt
bin_dir=$project_name/bin
qc_file=$project_name/gender_qc.txt
gender_file=$project_name/gender_classification.txt
cov_file=$project_name/mean_coverage.txt
docker run -v "${output_path}:/output" -w "/output" -it --rm \
    tomas81/cnest:dev step3 \
    --indextab $index_tab \
    --bindir $bin_dir \
    --qc $qc_file \
    --gender $gender_file \
    --cov $cov_file
```

### Run Step4
```bash
index_tab=$project_name/index_tab.txt
bin_dir=$project_name/bin
cor_dir=$project_name/cor
logr_dir=$project_name/logr
rbin_dir=$project_name/rbin
gender_file=$project_name/gender_classification.txt
docker run -v "${output_path}:/output" -w "/output" -it --rm \
    tomas81/cnest:dev step4 \
    --indextab $index_tab \
    --bindir $bin_dir \
    --cordir $cor_dir \
    --logrdir $logr_dir \
    --rbindir $rbin_dir \
    --sample $sample_name \
    --gender $gender_file \
    --batch 1000
```

### Run Step5
```bash
index_tab=$project_name/index_tab.txt
bin_dir=$project_name/bin
cor_dir=$project_name/cor
cnv_dir=$project_name/cnv
rbin_dir=$project_name/rbin
gender_file=$project_name/gender_classification.txt
cov_file=$project_name/mean_coverage.txt
docker run -v "$PWD:$PWD" -w "$PWD" -it --rm \
    tomas81/cnest:dev step5 \
    --indextab $index_tab \
    --rbindir $rbin_dir \
    --cordir $cor_dir \
    --cnvdir $cnv_dir \
    --cov    $cov_file \
    --sample $sample_name \
    --gender $gender_file \
    --batch 300
```
