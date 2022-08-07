n_particle=100
n_iter=1000
save=1

export OMP_NUM_THREADS=8
make && ./main $n_particle $n_iter $save

python make_animation.py $n_particle