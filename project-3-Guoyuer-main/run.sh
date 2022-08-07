module load cuda
nvcc -O2 main.cu -o main_cuda -lm -arch=sm_60
gcc -O2 main.c -o main_cpu -lm -std=c99

n_thread=128
for n_ray in 100000000 200000000 400000000 800000000 1600000000
do
    ./main_cuda $n_ray 1000 $n_thread
    ./main_cpu $n_ray 1000
done

