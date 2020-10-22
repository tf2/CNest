########## SETUP #############
## Note: this simply defines the file systems root for this work and copy some needed resources into the rigt place 
docker run -v /Users/tomas/projects/test_bams/:/input_location -v /Users/tomas/projects/cnv-docker/output_location:/output_location -it cnv:1.0 ./resources/run step1 cnv-project

########## CONVERT #############
## Note: here its best to provide input bam file pathes - needs to be relative to <INPUT DIR> - and the linked UKBB ID - this is important to link at this stage and will be carried through the later stages
docker run -v /Users/tomas/projects/test_bams/:/input_location -v /Users/tomas/projects/cnv-docker/output_location:/output_location -it cnv:1.0 ./resources/run step2 cnv-project 1000704_23163_0_0.bam 1000704

########## RUN - each step to run completion across all input samples #############
## Note: i.e. we need to run all samples at each step and we cant move on it the nest step until all have completed - here we just use an index 1:number_of samples - indenpendant job for every sample and every step
### EVERY STEP (3-6) MUST COMPLETE ACROSS ALL SAMPLES BEFORE THE NEXT STEP IS RUN
# i.e.
sample_index=1
docker run -v /Users/tomas/projects/test_bams/:/input_location -v /Users/tomas/projects/cnv-docker/output_location:/output_location -it cnv:1.0 ./resources/run step3 cnv-project ${sample_index}
## WAIT FOR ALL SAMPLES
docker run -v /Users/tomas/projects/test_bams/:/input_location -v /Users/tomas/projects/cnv-docker/output_location:/output_location -it cnv:1.0 ./resources/run step4 cnv-project ${sample_index}
## WAIT FOR ALL SAMPLES
docker run -v /Users/tomas/projects/test_bams/:/input_location -v /Users/tomas/projects/cnv-docker/output_location:/output_location -it cnv:1.0 ./resources/run step5 cnv-project ${sample_index}
## WAIT FOR ALL SAMPLES
docker run -v /Users/tomas/projects/test_bams/:/input_location -v /Users/tomas/projects/cnv-docker/output_location:/output_location -it cnv:1.0 ./resources/run step5 cnv-project ${sample_index}
## WAIT FOR ALL SAMPLES
docker run -v /Users/tomas/projects/test_bams/:/input_location -v /Users/tomas/projects/cnv-docker/output_location:/output_location -it cnv:1.0 ./resources/run step6 cnv-project ${sample_index}
## WAIT FOR ALL SAMPLES

# COMPLETE


### Example running with singularity image

### docker and singularity image creation - note creating the old style .img file 
docker save cnv:1.0 | gzip > cnv_1.0.tar.gz
docker run -v /var/run/docker.sock:/var/run/docker.sock  -v /Users/tomas/projects/cnv-docker2singularity/simg:/output --privileged -t --rm singularityware/docker2singularity:1.11 cnv:1.0

### Running using singularity on EBI cluster
project_name="cnv-project"
run_script="/resources/run"
singularity_image="cnv_1.0-2020-07-19-3a6e6c3510fe.img"

# step 1 - create output directory structure and resources - set your path to bam files here (root path) and the desired location of the output files
bam_path=/hps/nobackup/birney/users/tomas/ukbb/exome_cram_cnv/process/
output_path=/hps/nobackup/birney/users/tomas/ukbb/exome_cram_cnv/example/
bind_dirs="${bam_path}:/input_location,${output_path}:/output_location"
singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} step1 ${project_name}

# step 2 - loop all bam files - submitting a single job for each - important here is that we link that file with the sample id - in this case the ukbb eid - make sure that link is correct
mem="-M10000"
job_out="/dev/null"
for f in $bam_path/*.bam
do
  echo "Processing $f file..."
  file=$(basename "$f")
  sample=`sed s/_23163_0_0.bam//g <<< "$file"`
  command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} step2 ${project_name} ${file} ${sample}"
  command="bsub ${mem} -o ${job_out} ${command}"
  echo $command
  eval $command
done

# step 3 - process all files now available in output directory structure
number_of_files=`ls -1q ${output_path}/${project_name}/bin | wc -l`
for (( sample_index=1; sample_index<=${number_of_files}; sample_index++ ))
do  
   command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} step3 ${project_name} ${sample_index}"
   command="bsub ${mem} -o ${job_out} ${command}"
   echo ${command}
   eval $command
done

#### At this step we need to do some QC - so stop and we discuss
# step 4 - process all files now available in output directory structure
number_of_files=`ls -1q ${output_path}/${project_name}/bin | wc -l`
search_size=1000

##########################################################################################
### Set gender file haf to male and female for example - DO NOT DO THIS!
gender_file="${output_path}/${project_name}/gender_qc.txt"
printf '%s\t%s\n' "sample gender" > ${gender_file}
half_n=$(( $number_of_files / 2 ))
search_size=$((${half_n}-1))
i=0
for f in ${output_path}/${project_name}/bin/*
do
  file=$(basename "$f")
  if [[ $i -lt $half_n ]]
  then
    printf '%s\t%s\n' "${file} male" >> ${gender_file}
  else
    printf '%s\t%s\n' "${file} female" >> ${gender_file}
  fi
  i=$((i+1))
done
##########################################################################################

# run step4 with gender_file (after qc checking!) and search_size
for (( sample_index=1; sample_index<=${number_of_files}; sample_index++ ))
do  
   command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} step4 ${project_name} ${sample_index} ${gender_file} ${search_size}"
   command="bsub ${mem} -o ${job_out} ${command}"
   echo ${command}
   eval $command
done


# run step5 conversion
number_of_files=`ls -1q ${output_path}/${project_name}/logr | wc -l`
for (( sample_index=1; sample_index<=${number_of_files}; sample_index++ ))
do  
   command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} step5 ${project_name} ${sample_index}"
   command="bsub ${mem} -o ${job_out} ${command}"
   echo ${command}
   eval $command
done


## HMM - step6 - cnv calling - last step until we do merging and assocition testing genome wide
number_of_files=`ls -1q ${output_path}/${project_name}/rbin | wc -l`
for (( sample_index=1; sample_index<=${number_of_files}; sample_index++ ))
do  
   command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} step6 ${project_name} ${sample_index} ${gender_file} ${search_size}"
   command="bsub ${mem} -o ${job_out} ${command}"
   echo ${command}
   eval $command
done










