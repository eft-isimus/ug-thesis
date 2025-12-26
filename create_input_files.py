#!/usr/bin/env python3
#=================
# Input file creation script
#=================

#-----------------
# Importing packages
#-----------------
import sys
import numpy as np
import ast
from collections import defaultdict
import sys
import json

inputs = json.loads(sys.argv[1])
path   = sys.argv[2]

# Example input:
# inputs = system_inputs   = {"ad_strength":0, "ad_cutoff": 2.5, "m": 10000, "polymer_seperation": 6, "N_p": 400, "N_m": 20, "bond_length": 1.12246, "k": 0, "T":0.1, "dt":0.001, "t_f":10000, "m":10000, "num_runs":1}

def create_input_file(run_index, **kwargs):
    ad_strength        = kwargs.get('ad_strength', 1)          # NOTE: NOT NEEDED, only there due to laziness
    ad_cutoff          = kwargs.get('ad_cutoff', 2.5)          # NOTE: NOT NEEDED, only there due to laziness
    sigma              = kwargs.get('sigma', 1)                # size of monomers
    polymer_seperation = kwargs.get('polymer_seperation', 10)  # distance between a polymer and its 4 nearest neighbors
    N_p                = kwargs.get('N_p', 2)                  # total number of polymers
    N_m                = kwargs.get('N_m', 10)                 # number of monomers in each polymer
    bond_length        = kwargs.get('bond_length', 1)          # SET TO 1 ALWAYS
    box_len            = (np.sqrt(N_p)) * polymer_seperation   # calculating box length to make box
    m                  = kwargs.get('m', 10_000)
    # the box size needs to be such that the edge and corner polymers also see the same brush
    box = [[-box_len/2, box_len/2], [-box_len/2, box_len/2], [-5, 2*N_m*bond_length]]  # box size
    k   = kwargs.get('k', 0)                                                                     # stiffness of polymers
    rho = N_p/(box_len**2)                                                         # density of polymer brush

    # first set of commands define the units, styles, sim box, computes etc.
    commands1 = f"""
#---------------
# Basic attributes of the simulation
#---------------
units lj
atom_style molecular
boundary p p f    
dimension 3

#---------------
# Bond, angle and pair styles
#---------------  
bond_style fene
angle_style cosine/squared

# WCA potential (LJ 12-6 with diff. cutoff and shift)
pair_style lj/cut {1.12246*sigma}                      # cutoff = 2^(1/6) sigma
pair_modify shift yes                                       # shift the potential to remove discontinuity

read_data ../{N_m}Nm_{N_p}Np_{np.round(rho, decimals=4)}rho.data

mass 1 1.0
bond_coeff 1 30.0 1.5 1.0 1.0
angle_coeff 1 {k} 180.0
pair_coeff 1 1 1.0 {sigma} {1.12246*sigma}        # set cutoff as 2^(1/6) sigma again

special_bonds fene

#---------------
# Defining required computes [to be used later]
#---------------  
# compute pe_all all pe
# compute pe_pair all pe pair
# compute pe_angle all pe angle
# compute pe_bond all pe bond
# compute comp_ang polymer_atoms angle/local theta

"""

    # ----------
    # NOTE:*NEW* Code for creating required polymers
    # ----------

    data_file = f"""
8000 atoms 
7600 bonds
7200 angles

1 atom types
1 bond types
1 angle types

{box[0][0]} {box[0][1]} xlo xhi
{box[1][0]} {box[1][1]} ylo yhi
{box[2][0]} {box[2][1]} zlo zhi

Atoms

"""
    for xi in range(int(np.sqrt(N_p))):
        for yi in range(int(np.sqrt(N_p))):
            x_coord = (xi * polymer_seperation) - ((int(np.sqrt(N_p)) - 1) * polymer_seperation)/2 # centers at 0
            y_coord = (yi * polymer_seperation) - ((int(np.sqrt(N_p)) - 1) * polymer_seperation)/2 # centers at 0
            mol_id = (xi * N_m + yi) + 1 # to ensure mol_id starts from 1 not 0

            for atm in range(N_m):
                z_coord = atm * bond_length
                # atom-id mol-id type x y z
                data_file += f"""{(mol_id - 1) * 20 + atm + 1} {mol_id} 1 {x_coord:.6f} {y_coord:.6f} {z_coord:.6f}\n"""

    bond_count = N_p * (N_m - 1)
    data_file += f"""\nBonds\n\n"""

    bond_id = 1
    for pol in range(N_p):
        first_atom = pol * N_m
        for bond in range(N_m - 1):
            a1 = first_atom + bond + 1
            a2 = first_atom + bond + 2
            data_file += f"""{bond_id} 1 {a1} {a2}\n"""
            bond_id += 1

    data_file += f"""\nAngles\n\n"""
    angle_id = 1
    for pol in range(N_p):
        first_atom = pol * N_m
        for angle in range(N_m - 2):
            a1 = first_atom + angle + 1
            a2 = first_atom + angle + 2
            a3 = first_atom + angle + 3
            data_file += f"""{angle_id} 1 {a1} {a2} {a3}\n"""
            angle_id += 1
    
    with open(path + f"/{N_m}Nm_{N_p}Np_{np.round(rho, decimals=4)}rho.data", "w") as f:
        f.write(data_file)
    
    commands2 = f"""         
    # creating adsorption groups, plane and walls
    group base_atoms id 1:{N_p*N_m}:{N_m}                        # fixing all of the base atoms
    group mobile_atoms subtract all base_atoms # creating group of atoms which can move

    # no adsorption plane for now                                                                
    # fix wall all wall/lj93 zlo 0.0 {ad_strength} 1.0 {ad_cutoff}    # adsorption wall
    # fix_modify wall energy yes

    # fix zwall all wall/reflect zlo 0 zhi {2*N_m*bond_length}                                  # impenetrable wall
    fix zwall all wall/reflect zlo 0 zhi {2*N_m*1.5}                                  # 1.5 is the max length of a bond in FENE
    
    """

    T = kwargs.get('T', 0)
    dt = kwargs.get('dt', 0.0001)
    t_f = kwargs.get('t_f', 0.0001)
    m = kwargs.get('m', 1000)

    commands3 = f"""
    # for post-processing in python
    dump 1 all custom {m} {N_m}Nm_{N_p}Np_{np.round(rho, decimals=4)}rho_{T}T_{t_f}tf_{run_index}ri.poly id mol type x y z    
    # for dumping atom coords in a sensible order (id 1, 2, 3, ...)
    dump_modify 1 sort id                           
    velocity all create {T} {np.random.randint(100_000, 999_999)} mom yes rot yes          # random seed = statistically ind. runs
    neighbor 1.0 bin
    neigh_modify delay 100 every 5 check yes                                                      # hard requirement for minimize
    comm_modify cutoff 3.0
    timestep {dt}

    fix 1 all langevin {T} {T} {100*dt} {np.random.randint(100_000, 999_999)}    # random seed = statistically ind. runs
    fix 2 all nve

    # fix just the base atoms to always be attached
    fix freeze_force base_atoms setforce 0.0 0.0 0.0
    fix_modify freeze_force energy no
    velocity base_atoms set 0.0 0.0 0.0
    restart {int(t_f/(dt*10))} {N_m}Nm_{N_p}Np_{np.round(rho, decimals=4)}rho_{T}T_{t_f}tf_{run_index}ri_restart.bin # creates restart files every tenth of the way
    run {int(t_f/dt)}
    """

    string = commands1 + commands2 + commands3
    
    with open(path + f"/{N_m}Nm_{N_p}Np_{np.round(rho, decimals=4)}rho_{T}T_{t_f}tf_{run_index}ri_input.sim", "w") as f:
        f.write(string)
    
    # return string

for i in range(int(inputs['num_runs'])):
    create_input_file(i + 1, **inputs)

