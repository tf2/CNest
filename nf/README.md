# A Nextflow pipeline for copy number estimation and analysis using CNest

## Run with LSF/SLURM and Singularity
Test with singularity version `3.5.0` and `3.7.0-1.el7`.

```bash
#####################
## Required arguments
#####################

module load singularity/3.5.0 
ref_path=/hps/research1/birney/users/shimin/reference/UKBB_FE/genome.fa
bed_path=/hps/research1/birney/users/shimin/CNest/index_files/ukbb_wes_index.bed
project_name=test_proj
design_file=design.csv
batch_size=1000
executor=lsf # set this to slurm for SLURM HPCs

##########

# Part 0 : Make index
nextflow run -with-report p1_report.html -with-trace p1_trace  -profile $executor \
    -r main -latest smshuai/CNest-nf \
    --part 0 \
    --project $project_name \
    --bed $bed_path

# Part 1 : Bait read count
nextflow run -with-report p1_report.html -with-trace p1_trace  -profile $executor \
    -r main -latest smshuai/CNest-nf \
    --part 1 \
    --project $project_name \
    --design $design_file \
    --ref $ref_path \
    --indexb ./results/$project_name/index.bed

# Part 2 : Gender QC
nextflow run -with-report p2_report.html -with-trace p2_trace -profile $executor \
    -r main -latest smshuai/CNest-nf \
    --part 2 \
    --project $project_name \
    --bindir ./results/$project_name/bin/ \
    --index ./results/$project_name/index_tab.txt

# ! Do QC here before continue

# Part 3 : logR-ratio calculation
nextflow run -with-report p3_report.html -with-trace p3_trace -profile $executor \
    -r main -latest smshuai/CNest-nf \
    --part 3 \
    --project $project_name \
    --bindir ./results/$project_name/bin/ \
    --index ./results/$project_name/index_tab.txt \
    --gender ./results/gender_classification.txt \
    --batch $batch_size

# Part 4 : HMM call
nextflow run -with-report p4_report.html -with-trace p4_trace -profile $executor \
    -r main -latest smshuai/CNest-nf \
    --part 4 \
    --project $project_name \
    --rbindir ./results/$project_name/rbin/ \
    --cordir ./results/$project_name/cor/ \
    --index ./results/$project_name/index_tab.txt \
    --gender ./results/gender_classification.txt \
    --cov ./results/mean_coverage.txt \
    --batch $batch_size
```

## Run with AWS
```
# Part 2 : Gender QC
# ! $binfiles is a .txt with one bin file path per line
nextflow run smshuai/CNest-nf \
    --part 2 \
    --project $project_name \
    --binlist $binfiles \
    --index ./results/$project_name/index_tab.txt


# Part 3 : logR-ratio calculation
nextflow run smshuai/CNest-nf \
    --part 3 \
    --project $project_name \
    --bindir $bin_dir \
    --index index_tab.txt \
    --gender gender_classification.txt \
    --batch $batch_size

```
