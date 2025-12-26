import os
import numpy as np

#--------------------------------
#--------------------------------
# Defining paths
#--------------------------------
#--------------------------------
cwd         = os.getcwd()              # e.g., /home/user/grandparent/parent/child
cwd_name = cwd.split('/')[-1]

#--------------------------------
#--------------------------------
# Defining required functions, and main analysis function
#--------------------------------
#--------------------------------
def removesuffix(s, suffix):
    return s[:-len(suffix)] if s.endswith(suffix) else s

def results_funct(NM, NP, filename):

    with open(filename, 'r') as f:
        lines = f.readlines()
    count = 0
    for i in range(len(lines)):
        if lines[i].startswith('ITEM: TIMESTEP'):
            count += 1

    r_e = np.zeros((count, NP)) # each row is one dumpstep, each col is one polymer

    count = 0 # so first dumpstep gets index 0
    for i in range(len(lines)): # reads through the whole file
        if lines[i].startswith('ITEM: TIMESTEP'):
            step = lines[i + 1].strip() # keeps track of what dumpstep is at
            count += 1
            if lines[i + 8].startswith('ITEM: ATOMS id mol type x y z'): #
                for j in range(NP): # loops over number of polymers
                    max_z = 0
                    min_z = 0
                    for k in range(NM):
                        parts = lines[i + 8 + (j*NM) + 1 + k].split()
                        z_pos = float(parts[5]) # first monomer z
                        if z_pos > max_z:
                            max_z = z_pos
                        if z_pos < min_z:
                            min_z = z_pos
                    r_e[count-1][j] = max_z - min_z

    mean_vals = np.mean(r_e, axis=1)

    return mean_vals[-1]

#--------------------------------
#--------------------------------
# Looping over every parameter directory, and every seed directory inside each one
#--------------------------------
#--------------------------------
for param_dir in os.listdir():
    if not os.path.isdir(param_dir):     # <-- skip anything that isn’t a folder
        continue

    results_file = os.path.join(param_dir, f"{param_dir}_results.txt")
    if os.path.exists(results_file):   # skip if already analysed
        continue

    with open(results_file, "w") as f:
        f.write(f"{param_dir}\n")

        seed_dirs = os.listdir(param_dir)

        for seed in seed_dirs:
            seed_path = os.path.join(param_dir, seed)
            if not os.path.isdir(seed_path):     # <-- skip anything that isn’t a folder
                continue

            output_files = [file for file in os.listdir(seed_path) if file.endswith(".poly")]
            if len(output_files) == 0:
                continue
            output_file = output_files[0]
            filename = os.path.join(seed_path, output_file)
            parts = output_file.split("_")

            NM = int(removesuffix(parts[0], 'Nm'))
            NP = int(removesuffix(parts[1], 'Np'))

            r_e_mean = results_funct(NM, NP, filename)

            f.write(f"{seed}, {r_e_mean}\n")


#--------------------------------
#--------------------------------
# Making long combined results file
#--------------------------------
#--------------------------------
combined_file = os.path.join(cwd, f"{cwd_name}_results_full.txt")
with open(combined_file, "w") as f_out:
    f_out.write("seed, r_e_mean\n") # first line for format
    for param_dir in os.listdir():
        if not os.path.isdir(param_dir): # skip if not a folder
            continue
        results_file = os.path.join(param_dir, f"{param_dir}_results.txt")
        if not os.path.exists(results_file): # if no results file, skip
            continue
        with open(results_file, "r") as f_in:
            f_out.writelines(f_in)


#--------------------------------
#--------------------------------
# Making shorter results file with min, max, mean, std
#--------------------------------
#--------------------------------
with open(combined_file, 'r') as f:
    lines = f.readlines()
lines = lines[1:] # skip first line which gives the format

# getting the number of parameter-value folders
params = [line.strip('\n') for line in lines if ',' not in line] # using [] list with generator
n_params = len(params)
param_indices = [index for index in range(len(lines)) if ',' not in lines[index]]

# calculating n_seeds from how many lines there are between 2 parameter lines
n_seeds = (param_indices[1] - param_indices[0]) - 1

# dictionary for storing values cleanly
dict = {}

# initialize loop with zeros
for i in param_indices:
    dict[lines[i].strip('\n')] = 0

# add values for all seeds in a list, add this to parameter key
for i in param_indices:
    vals_list = []
    for j in range(n_seeds):
        line = lines[i + j + 1]
        line = line.split()[1]
        vals_list.append(float(line.strip('\n')))
    dict[lines[i].strip('\n')] = vals_list

short_file = os.path.join(cwd, f"{cwd_name}_results.txt")
with open(short_file, 'w') as f:
    for i in param_indices:
        key = lines[i].strip()

        f.write(f"{key}\n")
        f.write(f"max = {np.max(dict[key])}\n")
        f.write(f"min = {np.min(dict[key])}\n")
        f.write(f"mean = {np.mean(dict[key])}\n")
        f.write(f"std = {np.std(dict[key])}\n")
                                                           
