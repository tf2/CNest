# Using Singularity
Getting help
```bash
# Pull from docker hub (dev version)
singularity pull -F docker://smshuai/cnest:dev
# Overall
singularity run /hps/nobackup2/singularity/shimin/cnest-dev.simg -h
# For a given step
singularity run /hps/nobackup2/singularity/shimin/cnest-dev.simg step2 -h
```

Step1
```bash
singularity run -B "/hps/research1/birney/users/shimin/CNest/index_files:/input_location,/hps/research1/birney/users/shimin/CNest/benchmark:/output_location" /hps/nobackup2/singularity/shimin/cnest-dev.simg step1 --project ukbb_wes --bed ukbb_wes_index.bed
```

Step2
```bash
# BAM
singularity run -B "${input_path}:/input_location,${output_path}/output_location" /hps/nobackup2/singularity/shimin/cnest-dev.simg python3.8 /resources/cnest.py step2 --project ukbb_wes --sample 'A' --input 'a.bam'

# CRAM (Need to mount ref path)
singularity run -B "${input_path}:/input_location,${output_path}/output_location,${ref_path}:/ref" /hps/nobackup2/singularity/shimin/cnest-dev.simg python3.8 /resources/cnest.py step2 --project ukbb_wes --sample 'A' --input 'a.cram'
```


# Using Docker
Getting help
```bash
# Overall
docker run -it --rm smshuai/cnest:dev -h
# For a given step
docker run -it --rm smshuai/cnest:dev step2 -h
```

Step1
```bash
docker run -v "${index_path}:/input_location" -v "${output_path}:/output_location" -it --rm smshuai/cnest:dev step1 --project ukbb_wes --bed index.bed
```

Step2
```bash
# BAM
docker run -v "${input_path}:/input_location" -v "${output_path}:/output_location" -it --rm smshuai/cnest:dev step2 --project ukbb_wes --sample 'A' --input 'a.bam'

# CRAM (Need to mount ref path)
docker run -v "${input_path}:/input_location" -v "${output_path}:/output_location" -v "${ref_path}:/ref" -it --rm smshuai/cnest:dev step2 --project ukbb_wes --sample 'A' --input 'a.bam'
```