#! /bin/bash
module load openmpi
make mpi

case "$1" in
    "strong")
        for node in 1 2 4 8 16 32 64
            do
                printf -v filename 'slurm_batch_%d' $node
                printf '#!/bin/bash
#SBATCH --job-name=cg_mpi
#SBATCH --output=cg_mpi_%d.out
#SBATCH --nodes=%d
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=caslake
#SBATCH --ntasks-per-node=1
module load openmpi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun --bind-to none ./main 2000\n' $node $node > $filename
                echo "$filename created"

                sbatch -A mpcs51087 $filename > /dev/null

                echo "$filename submitted"
                rm $filename
            done
    ;;



    "weak")
        for node in 1 2 4 8 16 32 64
        do
            n=`echo "scale=6; sqrt($node)*600" | bc`
            printf -v filename 'slurm_batch_%d_weak' $node
            printf '#!/bin/bash
#SBATCH --job-name=cg_mpi
#SBATCH --output=cg_mpi_%d_weak.out
#SBATCH --nodes=%d
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=caslake
#SBATCH --ntasks-per-node=1
module load openmpi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun --bind-to none ./main %.0f\n' $node $node $n> $filename
            echo "$filename created"
            sbatch -A mpcs51087 $filename > /dev/null
            echo "$filename submitted"
            rm $filename

        done
    ;;

esac



# make test && sbatch -A mpcs51087 testbatch
