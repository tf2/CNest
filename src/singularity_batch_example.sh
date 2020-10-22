########################################################################
# Titel: Example running ngscnv with a singularity image
# Author: Tomas W Fitzgerald
# Date: 21/07/2020
# Email: tomas@ebi.ac.uk
########################################################################

### Running using singularity on EBI cluster
singularity_image="/nfs/software/birney/CNest/cnv-docker-cnv1.0.simg"
project_name="cnv-project"
run_script="/resources/run"

## Note: this is an important parameter for a 50K sample run keep this set to 1000
search_size=1000

# Job schulder parameters
mem="-M10000"
job_out="/dev/null"

# BIND directories - set your bind dirs here - set bam_path to bam file location (root path) and set output_path to the desired output location (root path)
# Note: do not modify the "bind_dirs" variable - only change "bam_path" and "output_path"
bam_path=/hps/nobackup/birney/users/tomas/ukbb/exome_cram_cnv/process/
output_path=/hps/nobackup/birney/users/tomas/ukbb/exome_cram_cnv/example/
bind_dirs="${bam_path}:/input_location,${output_path}:/output_location"

# STEP 1 - create output directory structure and resources
step="step1"
singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} ${step} ${project_name}

# STEP 2 - linking all bam files - submitting a single job for each - important here is that we link each bam file with its sample id - in this case the ukbb eid - make sure that the link is correct!
# Note: this example assumes that all bam files are stored within the same directory - that is unlikely - so here you need to write something that finds the individual bam files
# IMPORTANT: bam file paths need to be relative to the bind location of the "bam_path" (input_location) and the sample ID must be set here - in this cases (its the ukb eid)
# Note: in this example we just remove a trailing string from the file names to get the correct ids - you need to make sure that the link to bam file and eid is correct at this step!
mem = "-M2000" # 2GB
step="step2"
for f in $bam_path/*.bam
do
  echo "Processing $f file..."
  bam_file=$(basename "$f")
  sample_eid=${bam_file/_23163_0_0.bam/}
  command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} ${step} ${project_name} ${bam_file} ${sample_eid}"
  command="bsub ${mem} -o ${job_out} ${command}"
  echo $command
  eval $command
done

#### STEP 3 - classify gender QC step
#### CAUTION!!!! After this step we need to do some QC - so after this step (step3) completes - stop and check 
step="step3"
command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} ${step} ${project_name}"
command="bsub ${mem} -o ${job_out} ${command}"
echo ${command}
eval $command

# STEP 4 - start processing all files in the directory structure - jobs are now based on an internal index.
# Note: all next steps could be submitted as jobs arrays - we need to know the index variable prior to submission and set it - its not possible to reference at submission (sadly!)
# TODO: understand which job schduler is in use and then we simply set the index based on that - this is probably the best way to do it and would make all next steps a single command!
step="step4"
number_of_files=`ls -1q ${output_path}/${project_name}/bin | wc -l`
for (( sample_index=1; sample_index<=${number_of_files}; sample_index++ ))
do  
   command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} ${step} ${project_name} ${sample_index}"
   command="bsub ${mem} -o ${job_out} ${command}"
   echo ${command}
   eval $command
done

# RUN STEP 5 - with gender classification (after QC checking!) and search_size 
# - this is a command that could be used once the correct gender_file is in place and we set the "search_size" parameter
step="step5"
mem = "-M10000" # 10GB
gender_file="${output_path}/${project_name}/gender_classification.txt"
number_of_files=`ls -1q ${output_path}/${project_name}/bin | wc -l`

for (( sample_index=1; sample_index<=${number_of_files}; sample_index++ ))
do  
   command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} ${step} ${project_name} ${sample_index} ${gender_file} ${search_size}"
   command="bsub ${mem} -o ${job_out} ${command}"
   echo ${command}
   eval $command
done


# STEP 6 - conversion to rbin
# Note: ultimately this could be combined with step7
step="step6"
mem = "-M2000" # 2GB
number_of_files=`ls -1q ${output_path}/${project_name}/logr | wc -l`
for (( sample_index=1; sample_index<=${number_of_files}; sample_index++ ))
do  
   command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} ${step} ${project_name} ${sample_index}"
   command="bsub ${mem} -o ${job_out} ${command}"
   echo ${command}
   eval $command
done


## STEP 7 - HMM - CNV calling - this is the last step until we do merging and association testing genome wide...
## NB - here is where we potentially need large amounts of memory - e.g. for search_size 1000 we probably need approx 50GB for each job
step="step7"
mem = "-M50000" # 50GB
number_of_files=`ls -1q ${output_path}/${project_name}/rbin | wc -l`
for (( sample_index=1; sample_index<=${number_of_files}; sample_index++ ))
do  
   command="singularity exec -B ${bind_dirs} ${singularity_image} ${run_script} ${step} ${project_name} ${sample_index} ${gender_file} ${search_size}"
   command="bsub ${mem} -o ${job_out} ${command}"
   echo ${command}
   eval $command
done

