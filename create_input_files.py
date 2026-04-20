#!/usr/bin/env python3
#-----------------
# Importing packages
#-----------------
import sys
import numpy as np
import json

inputs = json.loads(sys.argv[1])
path   = sys.argv[2]

# inputs = system_inputs   = {"ad_strength":0, "ad_cutoff": 2.5, "m": 10000, "polymer_seperation": 6, 
# "N_p": [12, 4], "N_m": [20, 40], "bond_length": 1.12246, "k": 0, "T":0.1, "dt":0.001, "t_f":100, 
# "m":10000, "num_runs":1, "mixed":1, "t_eq":50}

def calculate_virus_positions(total_N_p, polymer_seperation, n_viruses):
    grid_size = int(np.sqrt(total_N_p))
    box_len = grid_size * polymer_seperation
    
    # Define 4 equi-spaced positions on the grid
    # These are at 1/4 and 3/4 positions in both x and y directions
    qs = box_len / 4 # the quarter-spacing
    positions_4 = [
        (-qs, -qs),  # bottom-left quadrant
        (qs, -qs),   # bottom-right quadrant
        (-qs, qs),   # top-left quadrant
        (qs, qs)  ]  # top-right quadrant
    
    # Return only the requested number of positions
    return positions_4[:n_viruses]

def create_input_file(run_index, **kwargs):
    ad_strength        = kwargs.get('ad_strength', 10)          # for virus adsorption onto grafting plane
    ad_cutoff          = kwargs.get('ad_cutoff', 2.5)          # for virus adsorption onto grafting plane
    sigma              = kwargs.get('sigma', 1)                # size of monomers
    polymer_seperation = kwargs.get('polymer_seperation', 10)  # distance between a polymer and its 4 nearest neighbors
    N_p                = kwargs.get('N_p', 4)                  # list of number of polymers of each type [small to large]
    N_m                = np.array(kwargs.get('N_m', [10]))     # list of number of monomers in each polymer [small to large]
    bond_length        = kwargs.get('bond_length', 1)          # SET TO 1 ALWAYS
    m                  = kwargs.get('m', 10_000)
    mixed              = kwargs.get('mixed', False)
    T                  = kwargs.get('T', 0)
    dt                 = kwargs.get('dt', 0.0001)
    t_f                = kwargs.get('t_f', 0.0001)
    t_eq               = kwargs.get('t_eq', t_f * 0.2)         # default 20% of total time for eq
    m_eq               = kwargs.get('m_eq', m * 10)            # default 10x lower frequency
    dz                 = kwargs.get('dz', 0.2)                 # bin size (delta z) for native density profile

    virus_height       = kwargs.get('virus_height', 50)        # height at which the viruses are added
    virus_sigma        = kwargs.get('virus_sigma', 50)         # size of the virus (for virus-virus repulsion)
    virus_number       = kwargs.get('virus_number', 4)         # number of viruses

    morse_D0           = kwargs.get('morse_D0', 5) # potential well depth
    morse_alpha        = kwargs.get('morse_alpha', 5) # controls potential well width (higher alpha = narrower well, shorter range)
    morse_r0           = kwargs.get('morse_r0', 1) # controls size of the virus (r0 is the minima of the potential)
    morse_cutoff       = kwargs.get('morse_cutoff', morse_r0 + 2.0) # potential cutoff (r0 + 2.0 is good enough for narrow well)

    # 1 = mixed, 0 = not mixed
    mixed = bool(mixed)
    total_N_p          = np.sum(N_p)
    box_len            = (np.sqrt(total_N_p)) * polymer_seperation   # calculating box length to make box
    rho                = total_N_p/(box_len**2)   # density of polymer brush
    
    input_data         = kwargs.get('input_data', None)       # eqbm. snapshot of brush to be used
    eq_snapshot        = kwargs.get('eq_snapshot', None)       # eqbm. snapshot of brush to be used
    
    if not mixed:
        n_m = N_m[0] # for convenience when the brush is not mixed
        
    if mixed:
        common_prefix = f'{N_m[0]}_{N_m[1]}_Nm_{N_p[0]}_{N_p[1]}_Np_{np.round(rho, decimals=4)}rho'
    else:
        common_prefix = f'{N_m[0]}_Nm_{N_p[0]}_Np_{np.round(rho, decimals=4)}rho'
    data_filename = f'{common_prefix}.data'
    seq_filename = f'{common_prefix}.seq'
    input_filename = f'{common_prefix}_{T}T_{t_f}tf_{run_index}ri_input.sim'
    dump_filename_eq = f'{common_prefix}_{T}T_{t_f}tf_{run_index}ri_eq.poly'
    dump_filename = f'{common_prefix}_{T}T_{t_f}tf_{run_index}ri.poly'
    dump_filename_v = f'{common_prefix}_{T}T_{t_f}tf_{run_index}ri_v.poly'
    early_density_filename = f'{common_prefix}_{T}T_{t_f}tf_{run_index}ri_early.profile'
    prod_density_filename = f'{common_prefix}_{T}T_{t_f}tf_{run_index}ri_prod.profile'
    
    # the box size needs to be such that the edge and corner polymers also see the same brush
    box = [[-box_len/2, box_len/2], [-box_len/2, box_len/2], [-5, 1.2*np.max(N_m)*bond_length]]  # box size
    k   = kwargs.get('k', [0])                                                                     # stiffness of polymers

    if virus_number > 0:
        virus_positions = calculate_virus_positions(total_N_p, polymer_seperation, virus_number)

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
"""
    if virus_number > 0:
        commands1 += f"""
# WCA potential (LJ 12-6 with diff. cutoff and shift)
# l-j cutoff = 1.122 * sigma, morse cutoff = 2.5
pair_style hybrid/overlay lj/cut {1.12246*sigma} morse {morse_cutoff}                      # cutoff = 2^(1/6) sigma
pair_modify shift yes
"""
    else:
        commands1 += f"""
# WCA potential (LJ 12-6 with diff. cutoff and shift)
# l-j cutoff = 1.122 * sigma, morse cutoff = 2.5
pair_style lj/cut {1.12246*sigma}                           # cutoff = 2^(1/6) sigma
pair_modify shift yes

"""
    if input_data != None:
        commands1 += f"""
read_data ../{input_data} extra/atom/types 1
read_dump ./{eq_snapshot} 0 x y z
""" #NOTE: temporary, only for testing rn, CHANGE LATER
        if mixed:
            commands1 += f"""
create_atoms 5 single {virus_positions[0][0]} {virus_positions[0][1]} {virus_height}
create_atoms 5 single {virus_positions[1][0]} {virus_positions[1][1]} {virus_height}
create_atoms 5 single {virus_positions[2][0]} {virus_positions[2][1]} {virus_height}
create_atoms 5 single {virus_positions[3][0]} {virus_positions[3][1]} {virus_height}

mass 1 1.0
mass 2 1.0
mass 3 1.0
mass 4 1.0
mass 5 1.0

"""
        
        elif not mixed:
            commands1 += f"""
create_atoms 3 single {virus_positions[0][0]} {virus_positions[0][1]} {virus_height}
create_atoms 3 single {virus_positions[1][0]} {virus_positions[1][1]} {virus_height}
create_atoms 3 single {virus_positions[2][0]} {virus_positions[2][1]} {virus_height}
create_atoms 3 single {virus_positions[3][0]} {virus_positions[3][1]} {virus_height}

mass 1 1.0
mass 2 1.0
mass 3 1.0

"""

    else:
        commands1 += f"""
read_data ../{data_filename}"""

        if mixed and virus_number > 0:
            commands1 += f"""
    mass 1 1.0
    mass 2 1.0
    mass 3 1.0
    mass 4 1.0
    mass 5 1.0

    """
        elif not mixed and virus_number > 0:
            commands1 += f"""
    mass 1 1.0
    mass 2 1.0
    mass 3 1.0

    """


        elif mixed and virus_number == 0:
            commands1 += f"""
    mass 1 1.0
    mass 2 1.0
    mass 3 1.0
    mass 4 1.0

    """
        elif not mixed and virus_number == 0:
            commands1 += f"""
    mass 1 1.0
    mass 2 1.0

    """

    if mixed and virus_number > 0:
        commands1 += f"""
bond_coeff 1 30.0 1.5 1.0 1.0
angle_coeff 1 {k[0]} 180.0
angle_coeff 2 {k[1]} 180.0

# Pair coefficients for WCA (polymer-polymer interactions)
pair_coeff 1 1 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 1 2 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 2 2 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 1 3 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 1 4 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 2 3 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 2 4 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 3 3 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 3 4 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 4 4 lj/cut 1.0 {sigma} {1.12246*sigma}

# virus-chain morse interaction
#format: pair_coeff 1 5 morse morse_D0 morse_alpha morse_r0 morse_cutoff
pair_coeff 1 5 morse {morse_D0[0]} {morse_alpha} {morse_r0} {morse_cutoff}
pair_coeff 2 5 morse {morse_D0[0]} {morse_alpha} {morse_r0} {morse_cutoff}
pair_coeff 3 5 morse {morse_D0[1]} {morse_alpha} {morse_r0} {morse_cutoff}
pair_coeff 4 5 morse {morse_D0[1]} {morse_alpha} {morse_r0} {morse_cutoff}

# # virus-chain WCA interaction
# pair_coeff 1 5 lj/cut 1.0 {virus_sigma} {1.12246*virus_sigma}
# pair_coeff 2 5 lj/cut 1.0 {virus_sigma} {1.12246*virus_sigma}
# pair_coeff 3 5 lj/cut 1.0 {virus_sigma} {1.12246*virus_sigma}
# pair_coeff 4 5 lj/cut 1.0 {virus_sigma} {1.12246*virus_sigma}

# virus-virus WCA interaction
pair_coeff 5 5 lj/cut 1.0 {virus_sigma} {1.12246*virus_sigma}
"""

    elif not mixed and virus_number > 0:
                commands1 += f"""
bond_coeff 1 30.0 1.5 1.0 1.0
angle_coeff 1 {k[0]} 180.0

# Pair coefficients for WCA (polymer-polymer interactions)
pair_coeff 1 1 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 1 2 lj/cut 1.0 {sigma} {1.12246*sigma}
pair_coeff 2 2 lj/cut 1.0 {sigma} {1.12246*sigma}

# virus-chain morse interaction
pair_coeff 1 3 morse {morse_D0[0]} {morse_alpha} {morse_r0} {morse_cutoff}
pair_coeff 2 3 morse {morse_D0[0]} {morse_alpha} {morse_r0} {morse_cutoff}

# virus-virus WCA interaction
pair_coeff 3 3 lj/cut 1.0 {virus_sigma} {1.12246*virus_sigma}
"""
    elif mixed and virus_number == 0:
        commands1 += f"""
bond_coeff 1 30.0 1.5 1.0 1.0
angle_coeff 1 {k[0]} 180.0
angle_coeff 2 {k[1]} 180.0

# Pair coefficients for WCA (polymer-polymer interactions)
pair_coeff 1 1 1.0 {sigma} {1.12246*sigma}
pair_coeff 1 2 1.0 {sigma} {1.12246*sigma}
pair_coeff 2 2 1.0 {sigma} {1.12246*sigma}
pair_coeff 1 3 1.0 {sigma} {1.12246*sigma}
pair_coeff 1 4 1.0 {sigma} {1.12246*sigma}
pair_coeff 2 3 1.0 {sigma} {1.12246*sigma}
pair_coeff 2 4 1.0 {sigma} {1.12246*sigma}
pair_coeff 3 3 1.0 {sigma} {1.12246*sigma}
pair_coeff 3 4 1.0 {sigma} {1.12246*sigma}
pair_coeff 4 4 1.0 {sigma} {1.12246*sigma}
"""        

    elif not mixed and virus_number == 0:
        commands1 += f"""
bond_coeff 1 30.0 1.5 1.0 1.0
angle_coeff 1 {k[0]} 180.0

# Pair coefficients for WCA (polymer-polymer interactions)
pair_coeff 1 1 1.0 {sigma} {1.12246*sigma}
pair_coeff 1 2 1.0 {sigma} {1.12246*sigma}
pair_coeff 2 2 1.0 {sigma} {1.12246*sigma}
"""        

    commands1 += """

special_bonds fene

"""
    if input_data == None:
        data_file = f"""
{np.sum(N_m * N_p) + virus_number} atoms 
{np.sum((N_m - 1) * N_p)} bonds
{np.sum((N_m - 2) * N_p)} angles
    """
    # need 2 atom types, one for base atoms and one for non-base atoms
        if not mixed and virus_number > 0:
            data_file += f""" 
3 atom types
1 bond types
1 angle types
    """

    # need 4 atom types for mixed brushes, two for base atoms and two for non-base atoms
        elif mixed and virus_number > 0:
            data_file += f""" 
5 atom types
1 bond types
2 angle types
    """

        elif not mixed and virus_number == 0:
            data_file += f""" 
2 atom types
1 bond types
1 angle types
    """

    # need 4 atom types for mixed brushes, two for base atoms and two for non-base atoms
        elif mixed and virus_number == 0:
            data_file += f""" 
4 atom types
1 bond types
2 angle types
    """

        data_file += f"""
{box[0][0]} {box[0][1]} xlo xhi
{box[1][0]} {box[1][1]} ylo yhi
{box[2][0]} {box[2][1]} zlo zhi

Atoms

    """
        if mixed:
            choice_array = [0]*N_p[0] + [1]*N_p[1]
            np.random.shuffle(choice_array)
            # create a sequence file to store which polymer is long
            with open(path + f"/{seq_filename}", "w") as f:
                    f.write(f'{choice_array}')

        mol_id = 1
        atom_id = 1
        for xi in range(int(np.sqrt(total_N_p))):
            for yi in range(int(np.sqrt(total_N_p))):
                # x-y coords remain same regardless of mixed brushes
                x_coord = (xi * polymer_seperation) - ((int(np.sqrt(total_N_p)) - 1) * polymer_seperation)/2 # centers at 0
                y_coord = (yi * polymer_seperation) - ((int(np.sqrt(total_N_p)) - 1) * polymer_seperation)/2 # centers at 0
                if not mixed:
                    for atm in range(n_m):
                        z_coord = atm * bond_length
                        # atom-id mol-id type x y z
                        # base atoms = type 2, non-base = type 1
                        if atm == 0:
                            data_file += f"""{atom_id} {mol_id} 2 {x_coord:.6f} {y_coord:.6f} {z_coord:.6f}\n"""
                        else:
                            data_file += f"""{atom_id} {mol_id} 1 {x_coord:.6f} {y_coord:.6f} {z_coord:.6f}\n"""
                        atom_id += 1
                    mol_id += 1
                
                # pol-1 = 1, base-atom-1 = 2, pol-2 = 3, base-2 = 4
                else:
                    num_monomers = N_m[choice_array[mol_id - 1]]
                    if num_monomers == N_m[0]:
                        for atm in range(num_monomers):
                            z_coord = atm * bond_length
                            # atom-id mol-id type x y z
                            if atm == 0:
                                data_file += f"""{atom_id} {mol_id} 2 {x_coord:.6f} {y_coord:.6f} {z_coord:.6f}\n"""
                            else:
                                data_file += f"""{atom_id} {mol_id} 1 {x_coord:.6f} {y_coord:.6f} {z_coord:.6f}\n"""
                            atom_id += 1
                    elif num_monomers == N_m[1]:
                        for atm in range(num_monomers):
                            z_coord = atm * bond_length
                            # atom-id mol-id type x y z
                            if atm == 0:
                                data_file += f"""{atom_id} {mol_id} 4 {x_coord:.6f} {y_coord:.6f} {z_coord:.6f}\n"""
                            else:
                                data_file += f"""{atom_id} {mol_id} 3 {x_coord:.6f} {y_coord:.6f} {z_coord:.6f}\n"""
                            atom_id += 1
                    mol_id += 1
        
        if mixed and virus_number > 0:
            for i in range(virus_number):
                data_file += f"""{atom_id} {mol_id} 5 {virus_positions[i][0]:.6f} {virus_positions[i][1]:.6f} {virus_height:.6f}\n"""
                atom_id += 1
                mol_id += 1
        
        elif not mixed and virus_number > 0:
            for i in range(virus_number):
                data_file += f"""{atom_id} {mol_id} 3 {virus_positions[i][0]:.6f} {virus_positions[i][1]:.6f} {virus_height:.6f}\n"""
                atom_id += 1
                mol_id += 1
        
        bond_count = {np.sum((N_m - 1) * N_p)}
        data_file += f"""\nBonds\n\n"""

        mol_id = 1 
        bond_id = 1
        first_atom = 1
        for pol in range(total_N_p):
            if not mixed:
                first_atom = pol * n_m
                for bond in range(n_m - 1):
                    a1 = first_atom + bond + 1
                    a2 = first_atom + bond + 2
                    data_file += f"""{bond_id} 1 {a1} {a2}\n"""
                    bond_id += 1
                mol_id += 1
            else:
                num_monomers = N_m[choice_array[mol_id - 1]]
                for bond in range(num_monomers - 1):
                    a1 = first_atom + bond 
                    a2 = first_atom + bond + 1
                    data_file += f"""{bond_id} 1 {a1} {a2}\n"""
                    bond_id += 1
                mol_id += 1
                first_atom += num_monomers

        data_file += f"""\nAngles\n\n"""
        mol_id = 1
        angle_id = 1
        first_atom = 1
        for pol in range(total_N_p):
            if not mixed:
                first_atom = pol * n_m
                for angle in range(n_m - 2):
                    a1 = first_atom + angle + 1
                    a2 = first_atom + angle + 2
                    a3 = first_atom + angle + 3
                    data_file += f"""{angle_id} 1 {a1} {a2} {a3}\n"""
                    angle_id += 1
                mol_id += 1
            else:
                angle_type = choice_array[mol_id - 1] + 1 # short polymers get angle 1, long get 2
                num_monomers = N_m[choice_array[mol_id - 1]]
                for angle in range(num_monomers - 2):
                    a1 = first_atom + angle 
                    a2 = first_atom + angle + 1
                    a3 = first_atom + angle + 2
                    data_file += f"""{angle_id} {angle_type} {a1} {a2} {a3}\n"""
                    angle_id += 1
                mol_id += 1
                first_atom += num_monomers
        with open(path + f"/{data_filename}", "w") as f:
            f.write(data_file)
    
    if mixed:
        commands2 = f"""         
# creating adsorption groups, plane and walls
group base_atoms type 2 4                        # fixing all of the base atoms 
"""
    else:
        commands2 = f"""         
# creating adsorption groups, plane and walls
group base_atoms type 2                        # fixing all of the base atoms 
"""     
    commands2 += f"""
group mobile_atoms subtract all base_atoms # creating group of atoms which can move

fix zwall all wall/reflect zlo 0 zhi {np.max(N_m) * bond_length + 2*virus_sigma}                                  # impenetrable wall to keep viruses from diffusing
# fix zwall all wall/reflect zlo 0 zhi {1.2*np.max(N_m)*1.5}                                  # 1.5 is the max length of a bond in FENE

    """
    if mixed and virus_number > 0:
        commands2 += f"""
group virus_atoms type 5 # creating group of virus atoms
"""

    elif not mixed and virus_number > 0:
        commands2 += f"""
group virus_atoms type 3 # creating group of virus atoms
"""

    if virus_number > 0:
        commands2 += f"""
fix wall virus_atoms wall/lj93 zlo 0.0 {ad_strength} {virus_sigma} {ad_cutoff}    # adsorption wall
fix_modify wall energy yes
"""

    # Calculate step counts ensuring they are perfect multiples of 10 for accurately rolling averaging dumps
    steps_eq = (int(t_eq/dt) // 10) * 10
    steps_prod = (int((t_f - t_eq)/dt) // 10) * 10

    commands3 = f"""
velocity all create {T} {np.random.randint(100_000, 999_999)} mom yes rot yes          # random seed = statistically ind. runs
neighbor 1.0 bin
neigh_modify delay 100 every 5 check yes                                               # hard requirement for minimize
comm_modify cutoff 3.0
timestep {dt}

fix 1 all langevin {T} {T} {100*dt} {np.random.randint(100_000, 999_999)}              # random seed = statistically ind. runs
fix 2 all nve

fix freeze_force base_atoms setforce 0.0 0.0 0.0                                       # fix just the base atoms to always be attached
fix_modify freeze_force energy no
velocity base_atoms set 0.0 0.0 0.0
restart {int(t_f/(dt*10))} {common_prefix}_{T}T_{t_f}tf_{run_index}ri_restart.bin # creates restart files every tenth of the way

# Compute 1D chunk grid mapping to compute density along Z
compute zchunks all chunk/atom bin/1d z lower {dz} units box

# Equilibration phase (low resolution dump)
dump 1 all custom {m_eq} {dump_filename_eq} id mol type x y z
dump_modify 1 sort id

# Early density profile (running average for equilibration phase)
fix early_density all ave/chunk 10 1 {steps_eq if steps_eq > 0 else 10} zchunks density/number ave running file {early_density_filename}

run {steps_eq}
unfix early_density

# Production phase (high resolution dump)
undump 1
dump 1 all custom {m} {dump_filename} id mol type x y z
dump_modify 1 sort id

# Production density profile (running average for entire production phase)
fix prod_density all ave/chunk 10 1 {steps_prod if steps_prod > 0 else 10} zchunks density/number ave running file {prod_density_filename}

"""

    if virus_number > 0:
        commands3+= f"""dump 2 virus_atoms custom {m} {dump_filename_v} id mol type x y z
dump_modify 2 sort id
"""

    commands3 += f"""    
run {steps_prod}
    """
        
    string = commands1 + commands2 + commands3
    
    with open(path + f"/{input_filename}", "w") as f:
        f.write(string)