# Using Singularity
## Getting help
```bash
# Tested on singularity version 3.7.0-1.el7
# Show help of available commands
singularity run docker://smshuai/cnest:dev2 -h
# Show help for a given step
singularity run docker://smshuai/cnest:dev2 step2 -h
```

## Run Step1
```bash
singularity run -B "${index_path}:/input,${output_path}:/output" --pwd /output/ docker://smshuai/cnest:dev2 step1 --project 'test_proj' --bed '/input/test.bed'
```

## Run Step2
```bash
# BAM
singularity run -B "${input_path}:/input,${output_path}:/output" --pwd /output/ docker://smshuai/cnest:dev2 step2 --project 'test_proj' --sample 'test_bam' --input '/input/test.bam'

# CRAM (Need to mount ref path)
singularity run -B "${input_path}:/input,${output_path}:/output,${ref_path}:/ref" --pwd /output/ docker://smshuai/cnest:dev2 step2 --project 'test_proj' --sample 'test_cram' --input '/input/test.cram'
```

# Using Docker
## Getting help
```bash
# Show help of available commands
docker run -it --rm smshuai/cnest:dev -h
# Show help for a given step
docker run -it --rm smshuai/cnest:dev step2 -h
```

## Run Step1
Note: `index.bed` must have the same chromosome names as the BAM/CRAM file.
```bash
docker run -v "${index_path}:/input" -v "${output_path}:/output" -w "/output" -it --rm smshuai/cnest:dev2 step1 --project test_proj --bed /wkdir/index.bed
```

## Run Step2
```bash
# BAM
docker run -v "${input_path}:/input" -v "${output_path}:/output" -w "/output" -it --rm smshuai/cnest:dev2 step2 --project test_proj --sample 'test_bam' --input '/input/test.bam'

# CRAM (Need to mount ref path)
docker run -v "${input_path}:/input" -v "${output_path}:/output" -v "${ref_path}:/ref" -w "/output" -it --rm smshuai/cnest:dev2 step2 --project test_proj --sample 'test_cram' --input '/input/test.cram'
```