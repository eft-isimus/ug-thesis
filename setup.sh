#!/bin/bash
source /etc/profile.d/modules.sh

module load compiler/gcc-11.2.0
module load compiler/python3.8 

python3 create_input_files.py

read -p "Enter path to param_dir: " path
read -p "Enter the compute node: " node

# python3 create_input_files.py
echo "$path" | ./move_files.sh
./make_jobs.sh $path $node

echo 'done'

