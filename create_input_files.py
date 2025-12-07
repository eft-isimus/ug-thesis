#!/usr/bin/env python3
#=================
# Input file creation script
#=================

#-----------------
# Importing packages
#-----------------
# import sys
import numpy as np
import ast
from collections import defaultdict
import sys
import json

inputs = json.loads(sys.argv[1])
path   = sys.argv[2]

# inputs = ast.literal_eval(input("Enter input parameters: "))
# path = input("Enter path to directory in which to make inputs: ")

#-----------------
# Defining required functions
#-----------------
def make_grid(N, sep):
    n = int(np.sqrt(N))
    coords = (np.arange(n) - (n-1)/2) * sep
    X, Y = np.meshgrid(coords, coords)
    return np.column_stack([X.ravel(), Y.ravel()])

# Example input:
# inputs = system_inputs   = {'ad_strength':0, 'ad_cutoff': 2.5, 'm': 10000, 'polymer_seperation': 6, 'N_p': 64, 'N_m': 10, 'bond_length': 1.12246, 'k': 0, 'T':0.1, 'dt':0.001, 't_f':10000, 'm':10000, 'num_runs':10}

def create_input_file(run_index, **kwargs):
    ad_strength        = kwargs.get('ad_strength', 1)          # NOTE: NOT NEEDED, only there due to laziness
    ad_cutoff          = kwargs.get('ad_cutoff', 2.5)          # NOTE: NOT NEEDED, only there due to laziness
    sigma              = kwargs.get('sigma', 1)                # size of monomers
    polymer_seperation = kwargs.get('polymer_seperation', 10)  # distance between a polymer and its 4 nearest neighbors
    N_p                = kwargs.get('N_p', 2)                  # total number of polymers
    N_m                = kwargs.get('N_m', 10)                 # number of monomers in each polymer
    bond_length        = kwargs.get('bond_length', 2)          # bonds between successive monomers
    box_len            = (np.sqrt(N_p)) * polymer_seperation   # calculating box length to make box
    # box_len = (np.sqrt(N_p) + 1) * polymer_seperation  # calculating box length to make box
    m                  = kwargs.get('m', 10_000)
    # the box size needs to be such that the edge and corner polymers also see the same brush
    box = [[-box_len/2, box_len/2], [-box_len/2, box_len/2], [-5, 2*N_m*bond_length]]  # box size
    # box = [[-box_len/2 - bond_length*N_m, box_len/2 + bond_length*N_m], 
                # [-box_len/2 - bond_length*N_m, box_len/2 + bond_length*N_m], 
                # [-5, 2*N_m*bond_length]]  # box size
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
    # boundary f f f
    
    dimension 3
    region box block {box[0][0]} {box[0][1]} {box[1][0]} {box[1][1]} {box[2][0]} {box[2][1]}
    create_box 1 box bond/types 1 angle/types 1 extra/bond/per/atom 2 extra/angle/per/atom 3 extra/special/per/atom 10

    #---------------
    # Bond, angle and pair styles
    #---------------  
    mass 1 1.0
    
    bond_style harmonic
    bond_coeff 1 100.0 {bond_length}

    angle_style cosine/squared
    angle_coeff 1 {k} 180.0

    # WCA potential (LJ 12-6 with diff. cutoff and shift)
    pair_style lj/cut {1.12246*sigma}                      # cutoff = 2^(1/6) sigma
    pair_modify shift yes                                       # shift the potential to remove discontinuity
    pair_coeff 1 1 1.0 {sigma} {1.12246*sigma}        # set cutoff as 2^(1/6) sigma again

    #---------------
    # Defining required computes
    #---------------  
    compute pe_all all pe
    compute pe_pair all pe pair
    compute pe_angle all pe angle
    compute pe_bond all pe bond
    # compute comp_ang polymer_atoms angle/local theta

    """

    # ----------
    # NOTE:*NEW* Code for creating required polymers
    # ----------

    mol_file = f"""
    # header section:
{N_m} atoms
{N_m - 1} bonds
{N_m - 2} angles

# body section:
Coords

""" # empty string
    
    for i in range(N_m):
        mol_file += f"""{i + 1} 0.00000 0.00000 {np.round(i * bond_length, decimals=5)}
"""
    mol_file += f"""
Bonds

"""
    for i in range(N_m - 1):
        # id    type    atom_1  atom_2
        mol_file += f"""{i + 1} 1 {i + 1} {i + 2}
"""
    mol_file += f"""
Types

"""
    for i in range(N_m):
        # id type
        mol_file += f"""{i + 1} 1 
"""    
    mol_file += f"""
Angles

"""
    for i in range(N_m - 2):
        # ID type atom1 atom2 atom3        
        mol_file += f"""{i + 1} 1 {i + 1} {i + 2} {i + 3}
"""
    with open(path + f"/{N_m}Nm_{N_p}Np_{np.round(rho, decimals=4)}rho.mol", "w") as f:
        f.write(mol_file)
    
    commands2 = f"""molecule polymer ../{N_m}Nm_{N_p}Np_{np.round(rho, decimals=4)}rho.mol

""" 

    for i in range(int(np.sqrt(N_p))):
        for j in range(int(np.sqrt(N_p))):
            commands2 += f"""create_atoms 0 single {(i * polymer_seperation) - ((np.sqrt(N_p) - 1) * polymer_seperation)/2} {(j * polymer_seperation)  - ((np.sqrt(N_p) - 1) * polymer_seperation)/2} {((N_m - 1)*bond_length)/2} mol polymer 123456 rotate 0 0 0 1
"""
    
    commands3 = f"""         
    # creating adsorption groups, plane and walls
    group base_atoms id 1:{N_p*N_m}:{N_m}                        # fixing all of the base atoms
    group mobile_atoms subtract all base_atoms # creating group of atoms which can move

    # no adsorption plane for now                                                                
    # fix wall all wall/lj93 zlo 0.0 {ad_strength} 1.0 {ad_cutoff}    # adsorption wall
    # fix_modify wall energy yes

    fix zwall all wall/reflect zlo 0 zhi {2*N_m*bond_length}                                  # impenetrable wall
    
    """

    T = kwargs.get('T', 0)
    dt = kwargs.get('dt', 0.0001)
    t_f = kwargs.get('t_f', 0.0001)
    m = kwargs.get('m', 1000)

    commands4 = f"""
    # for post-processing in python
    dump 1 all custom {m} {N_m}Nm_{N_p}Np_{np.round(rho, decimals=4)}rho_{T}T_{t_f}tf_{run_index}ri.poly id mol type x y z    
    # for dumping atom coords in a sensible order (id 1, 2, 3, ...)
    dump_modify 1 sort id                           
    # dump 2 all custom {m} {N_m}Nm_{N_p}Np_{np.round(rho, decimals=4)}rho_{T}T_{t_f}tf_{run_index}ri.lammpstrj id type x y z
    # dump_modify 2 sort id
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

    string = commands1 + commands2 + commands3 + commands4
    
    with open(path + f"/{N_m}Nm_{N_p}Np_{np.round(rho, decimals=4)}rho_{T}T_{t_f}tf_{run_index}ri_input.sim", "w") as f:
        f.write(string)
    
    # return string

for i in range(int(inputs['num_runs'])):
    create_input_file(i + 1, **inputs)
