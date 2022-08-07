#! /bin/bash

for node in 1 2 4 8
    do
        printf -v slurm_pure 'slurm_pure_mpi_%d_node' $node
        printf '#!/bin/sh
#SBATCH --job-name=nbody_mpi_pure
#SBATCH --output=pure_mpi_%d_node.out
#SBATCH --nodes=%d
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=caslake
#SBATCH --ntasks-per-node=24
#SBATCH --account=mpcs51087
#SBATCH --exclusive

module load openmpi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
make main_mpi && mpirun --bind-to none ./main_mpi 102400 10 0' $node $node> $slurm_pure
        echo "$slurm_pure created!"
        sbatch -A mpcs51087 $slurm_pure > /dev/null
        echo "$slurm_pure submitted!"


        
        printf -v slurm_mix 'slurm_mix_mpi_%d_node' $node
        printf '#!/bin/sh
#SBATCH --job-name=nbody_mpi_mix
#SBATCH --output=mix_mpi_%d_node.out
#SBATCH --nodes=%d
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=24
#SBATCH --partition=caslake
#SBATCH --ntasks-per-node=1
#SBATCH --account=mpcs51087
#SBATCH --exclusive


module load openmpi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
make main_mpi && mpirun --bind-to none ./main_mpi 102400 10 0' $node $node> $slurm_mix
        echo "$slurm_mix created!"
        sbatch -A mpcs51087 $slurm_mix > /dev/null
        echo "$slurm_mix submitted!"

    done
