#! /bin/bash
module load openmpi

case "$1" in
    "strong") 
        echo "run strong"
        if make hybrid; then 
            sbatch -A mpcs51087 batchfile1
            sbatch -A mpcs51087 batchfile4
            sbatch -A mpcs51087 batchfile16
        fi
    ;; 
    "strong_no_omp")
        if make hybrid; then 
            sbatch -A mpcs51087 batchfile1_no_omp
            # sbatch -A mpcs51087 batchfile4_no_omp
            # sbatch -A mpcs51087 batchfile16_no_omp
        fi
    ;;
    "weak")
        if make hybrid; then
            sbatch -A mpcs51087 batchfile1_weak
            sbatch -A mpcs51087 batchfile4_weak
            sbatch -A mpcs51087 batchfile16_weak
        fi
    ;;
    "weak_no_omp")
        if make hybrid; then
            sbatch -A mpcs51087 batchfile1_weak_no_omp
            sbatch -A mpcs51087 batchfile4_weak_no_omp
            sbatch -A mpcs51087 batchfile16_weak_no_omp
        fi
    ;;
    "first_order")
        if make first_order; then
            sbatch -A mpcs51087 batchfile_first_order
        fi
    ;;
    "second_order")
        if make second_order; then
            sbatch -A mpcs51087 batchfile_second_order
        fi
    ;;
    "test")
        echo "run test"
        make test && sbatch -A mpcs51087 test_batchfile
    ;;
    "demo")
        echo "generate movie"
        if make hybrid; then 
            sbatch -A mpcs51087 batchfile_demo
        fi
    ;;
esac