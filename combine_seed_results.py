#----------------------
#NOTE: This script combines the production results for all seeds for that parameter value
#NOTE: It will be placed in the RUN directory (the parent directory for all par dirs)
#----------------------
import os
import numpy as np
import ast
import sys

#----------------------
# 1.1 Find if brush is mixed and set up header (gives format)
#----------------------
run_dir      = os.getcwd()
run_dir_name = os.path.basename(run_dir)

par_dirs = [dir for dir in os.listdir(run_dir) if os.path.isdir(os.path.join(run_dir, dir))]
mixed = any(f.endswith('.seq') for f in os.listdir(os.path.join(run_dir, par_dirs[0])))

# set up header and output file
output_file = os.path.join(run_dir, f'{run_dir_name}_results.txt')
if mixed:
    header = 'par_value,all_mean,all_std,long_mean,long_std,short_mean,short_std'
else:
    header = 'par_value,all_mean,all_std'

# write header once before the loop
with open(output_file, 'w') as f:
    f.write(header + '\n')

#----------------------
# 1.2 Main loop which goes through each par_dir
#----------------------
for par_dir_i in range(len(par_dirs)):
    par_dir = par_dirs[par_dir_i]
    par_dir_path = os.path.join(run_dir, par_dir)
    seed_dirs = [dir for dir in os.listdir(par_dir_path) if os.path.isdir(os.path.join(par_dir_path, dir))]
    N_seeds = len(seed_dirs)

    if mixed:
        par_array = np.zeros((N_seeds, 6))
    else:
        par_array = np.zeros((N_seeds, 2))
    
    for seed_dir_i in range(len(seed_dirs)):
        seed_dir = seed_dirs[seed_dir_i]
        seed_dir_path = os.path.join(par_dir_path, seed_dir)
        prod_file_path = os.path.join(seed_dir_path, 'prod_results.txt')

        par_array[seed_dir_i] = np.genfromtxt(prod_file_path, delimiter=',', skip_header=2)

    all_mean = np.mean(par_array[:, 0])
    all_std  = np.std(par_array[:, 0])

    if not mixed:
        with open(output_file, 'a') as f:
            f.write(f'{par_dir},{all_mean},{all_std}\n')
    if mixed:
        long_mean  = np.mean(par_array[:, 2])
        long_std   = np.std(par_array[:, 2])
        short_mean = np.mean(par_array[:, 4])
        short_std  = np.std(par_array[:, 4])
        with open(output_file, 'a') as f:
            f.write(f'{par_dir},{all_mean},{all_std},{long_mean},{long_std},{short_mean},{short_std}\n')