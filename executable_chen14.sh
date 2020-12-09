#!/bin/bash

job_num=$1 #job number currently being processed
input_file="chen14_kin_input_v2.csv" #input file to process
main_dir="chen14" #name of dir to store files in
main_name="chen14_unc" #name-tag to assign to output files
num_per=2 #stars to process per job
base_add=0 #num jobs prev submitted (if multiple submit files needed)

start=$(echo "$num_per*($job_num+$base_add)" | bc) #first star to run this job
stop=$(echo "$start+$num_per" | bc) #last star to run this job

condor_base=/afs/crc.nd.edu/user/s/sdietz/condor_stuff/

#make new dir for every 1000 files
dir_num=$(echo "($job_num+$base_add)/1000./1" | bc)
dir_name="${main_dir}/runs_$dir_num/"

file_num=$(echo "$job_num+$base_add" | bc)
file_name="${main_name}_run_$file_num" #where to store kin output

#save stuff in a temporary dir (transfer at end)
tmp_dir=$(/usr/bin/mktemp -d)
tmp_save_dir="$tmp_dir/$dir_name"
mkdir -p $tmp_save_dir

module load python/2.7.14

# Copy necessary files over for local processing (/tmp)
cp ${condor_base}/dietz_mini_pipeline_v1.py ${tmp_dir}
cp ${condor_base}/${input_file} ${tmp_dir}
cp ${condor_base}/staeckel_orbit.py ${tmp_dir}

# Begin local processing
cd ${tmp_dir}

python dietz_mini_pipeline_v1.py $input_file $file_name $tmp_save_dir $start \
    $stop

# Copy over tmp_dir results into AFS
cp -r ${tmp_save_dir}/ ${condor_base}/${main_dir}/

# cleanup tmp_dir
cd -
rm -rf $tmp_dir

