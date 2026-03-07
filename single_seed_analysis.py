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
par_dir         = os.path.dirname(seed_dir) # parent / parameter directory containing .seq file
par_dir_name    = par_dir.split('/')[-1]

#----------------------
# 1. determine if seed is for a mixed brush, open seq file if needed
#----------------------
seq_files = [file for file in os.listdir(par_dir) if file.endswith('.seq')]
mixed = bool(len(seq_files)) # true if there is one seq file
seq_file = os.path.join(par_dir, seq_files[0])
if mixed:
        with open(seq_file, 'r') as f:
            seq = f.readlines()
        seq = np.array(ast.literal_eval(seq[0])) # sequence of long and short polymers (0 = short)

#----------------------
# 2. Opening both eq and prod dumps
#----------------------
eq_file = [file for file in os.listdir(seed_dir) if file.endswith('eq.poly')][0]
prod_file = [file for file in os.listdir(seed_dir) if file.endswith('.poly') and not file.endswith('eq.poly')][0]

#----------------------
# 2.1 Defining analyse function
#----------------------
def process_dump(filename, seq=None, mixed=False, run_type='eq'):
    """
    Parameters:
        filename   : path to LAMMPS dump file
        n_molecules: total number of molecules in the brush
        seq        : array of length n_molecules, 1=long, 0=short (only needed if mixed=True)
        mixed      : boolean, whether the brush is mixed
        run_type   : 'eq'   -> return per-timestep arrays of shape (n_timesteps, 3): [timestep, mean_h, std_h]
                     'prod' -> return scalar (mean_h, std_h) averaged over all timesteps

    Returns (mixed=True):
        eq:   all_heights, long_heights, short_heights  -- each (n_timesteps, 3)
        prod: (all_mean, all_std), (long_mean, long_std), (short_mean, short_std)

    Returns (mixed=False):
        eq:   all_heights -- (n_timesteps, 3)
        prod: (all_mean, all_std)
    """
    all_heights   = []
    long_heights  = []
    short_heights = []

    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):

        if lines[i].startswith('ITEM: TIMESTEP'):
            timestep = int(lines[i + 1].strip())

            # scan forward to ITEM: ATOMS
            while not lines[i].startswith('ITEM: ATOMS'):
                i += 1
            i += 1  # step past ITEM: ATOMS header to first atom line

            # read all atom lines for this timestep
            atom_lines = []
            while i < len(lines) and not lines[i].startswith('ITEM:'):
                atom_lines.append(lines[i].split())
                i += 1
            # i now points to the next ITEM: line, do NOT increment at bottom of loop

            # group z-coords by mol-id
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

            continue  # i already at next ITEM: line, skip i += 1

        i += 1

    all_heights = np.array(all_heights)

    if run_type == 'eq':
        if mixed:
            long_heights  = np.array(long_heights)
            short_heights = np.array(short_heights)
            return all_heights, long_heights, short_heights
        else:
            return all_heights

    elif run_type == 'prod':
        all_result = (np.mean(all_heights[:, 1]), np.std(all_heights[:, 1]))
        if mixed:
            long_heights  = np.array(long_heights)
            short_heights = np.array(short_heights)
            long_result  = (np.mean(long_heights[:,  1]), np.std(long_heights[:,  1]))
            short_result = (np.mean(short_heights[:, 1]), np.std(short_heights[:, 1]))
            return all_result, long_result, short_result
        else:
            return all_result

#----------------------
# 2.2 Using the function on both files and saving outputs
#----------------------
if mixed:
    eq_all, eq_long, eq_short     = process_dump(eq_file,   seq=seq, mixed=True, run_type='eq')
    prod_all, prod_long, prod_short = process_dump(prod_file, seq=seq, mixed=True, run_type='prod')
    np.savetxt('eq_results.txt',   np.hstack([eq_all, eq_long[:,1:], eq_short[:,1:]]), delimiter=',', header='mixed\ntimestep,all_mean,all_std,long_mean,long_std,short_mean,short_std', comments='')
    np.savetxt('prod_results.txt', np.array([list(prod_all) + list(prod_long) + list(prod_short)]), delimiter=',', header='mixed\nall_mean,all_std,long_mean,long_std,short_mean,short_std', comments='')
else:
    eq_all   = process_dump(eq_file,   mixed=False, run_type='eq')
    prod_all = process_dump(prod_file, mixed=False, run_type='prod')
    np.savetxt('eq_results.txt',   eq_all,                        delimiter=',', header='non-mixed\ntimestep,all_mean,all_std', comments='')
    np.savetxt('prod_results.txt', np.array([list(prod_all)]),    delimiter=',', header='non-mixed\nall_mean,all_std', comments='')

#NOTE: this gives an eq results file which can be used to plot, and a prod file which has only one line in it