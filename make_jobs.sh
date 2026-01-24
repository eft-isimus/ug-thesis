#!/bin/bash

main_dir="$1" 
comp_node="$2"

cd "$main_dir" || exit

for dir in */; do
    cd "$dir" || continue
    simfile=$(ls *.sim 2>/dev/null)
    if [ -n "$simfile" ]; then
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
EOF
   # qsub job.pbs
    fi
    cd ..
done
