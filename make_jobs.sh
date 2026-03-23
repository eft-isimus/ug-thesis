#!/bin/bash

main_dir="$1" 
comp_node="$2"

cd "$main_dir" || exit

for dir in */; do
    cd "$dir" || continue
    simfile=$(ls *.sim 2>/dev/null)
    if [ -n "$simfile" ]; then
        cat > single_seed_analysis.py <<'PYEOF'
#----------------------
#NOTE: This script analyses a single seed's output files
#NOTE: It will be placed in the seed file
#----------------------
import os
import numpy as np
import matplotlib.pyplot as plt
import ast
import sys

seed_dir        = os.getcwd()
seed_dir_name   = seed_dir.split('/')[-1]
par_dir         = os.path.dirname(seed_dir)
par_dir_name    = par_dir.split('/')[-1]

seq_files = [file for file in os.listdir(par_dir) if file.endswith('.seq')]
mixed = bool(len(seq_files))
seq_file = seq_files[0]
if mixed:
    with open(seq_file, 'r') as f:
        seq = f.readlines()
    seq = np.array(ast.literal_eval(seq[0]))

eq_file   = [file for file in os.listdir(seed_dir) if file.endswith('eq.poly')][0]
prod_file = [file for file in os.listdir(seed_dir) if file.endswith('.poly') and not file.endswith('eq.poly')][0]

def process_dump(filename, seq=None, mixed=False, run_type='eq'):
    all_heights   = []
    long_heights  = []
    short_heights = []

    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith('ITEM: TIMESTEP'):
            timestep = int(lines[i + 1].strip())
            while not lines[i].startswith('ITEM: ATOMS'):
                i += 1
            i += 1
            atom_lines = []
            while i < len(lines) and not lines[i].startswith('ITEM:'):
                atom_lines.append(lines[i].split())
                i += 1
            mol_z = {}
            for atom in atom_lines:
                mol_id = int(atom[1])
                z      = float(atom[5])
                if mol_id not in mol_z:
                    mol_z[mol_id] = []
                mol_z[mol_id].append(z)
            heights = [max(mol_z[mid]) for mid in sorted(mol_z.keys())]
            all_heights.append([timestep, np.mean(heights), np.std(heights)])
            if mixed:
                long_h  = [h for h, s in zip(heights, seq) if s == 1]
                short_h = [h for h, s in zip(heights, seq) if s == 0]
                long_heights.append( [timestep, np.mean(long_h),  np.std(long_h)])
                short_heights.append([timestep, np.mean(short_h), np.std(short_h)])
            continue
        i += 1

    all_heights = np.array(all_heights)
    if run_type == 'eq':
        if mixed:
            return all_heights, np.array(long_heights), np.array(short_heights)
        else:
            return all_heights
    elif run_type == 'prod':
        all_result = (np.mean(all_heights[:, 1]), np.std(all_heights[:, 1]))
        if mixed:
            long_heights  = np.array(long_heights)
            short_heights = np.array(short_heights)
            return all_result, (np.mean(long_heights[:, 1]), np.std(long_heights[:, 1])), (np.mean(short_heights[:, 1]), np.std(short_heights[:, 1]))
        else:
            return all_result

if mixed:
    eq_all, eq_long, eq_short       = process_dump(eq_file,   seq=seq, mixed=True, run_type='eq')
    prod_all, prod_long, prod_short = process_dump(prod_file, seq=seq, mixed=True, run_type='prod')
    np.savetxt('eq_results.txt',   np.hstack([eq_all, eq_long[:,1:], eq_short[:,1:]]), delimiter=',', header='mixed\ntimestep,all_mean,all_std,long_mean,long_std,short_mean,short_std', comments='')
    np.savetxt('prod_results.txt', np.array([list(prod_all) + list(prod_long) + list(prod_short)]), delimiter=',', header='mixed\nall_mean,all_std,long_mean,long_std,short_mean,short_std', comments='')
else:
    eq_all   = process_dump(eq_file,   mixed=False, run_type='eq')
    prod_all = process_dump(prod_file, mixed=False, run_type='prod')
    np.savetxt('eq_results.txt',   eq_all,                     delimiter=',', header='non-mixed\ntimestep,all_mean,all_std', comments='')
    np.savetxt('prod_results.txt', np.array([list(prod_all)]), delimiter=',', header='non-mixed\nall_mean,all_std', comments='')
PYEOF

        cat > job.pbs <<EOF
#!/bin/bash
#PBS -N ${main_dir//\//_}_${dir%/}
#PBS -l nodes=1:ppn=8
#PBS -l walltime=15:00:00
#PBS -l mem=4gb
#PBS -l host=compute${comp_node}
#PBS -q cpu
#PBS -j oe
cd \$PBS_O_WORKDIR
module load compiler/gcc-12.0
module load compiler/openmpi-3.1.6
module load apps/lammps-stable_2Aug2023
mpirun -np 4 lmp_mpi < $simfile
python3 single_seed_analysis.py
EOF
    fi
    cd ..
done