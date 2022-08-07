#!/bin/bash

thread_num=(1 2 4 8 12 16 24 32)
#thread_num=(1)
if [ ! -f "omp_advection" ]; then
  make
else
  echo "omp_advection found"
fi

if [ ! -f "../milestone-1/advection" ]; then
  cd ../milestone-1 && make
else
  echo "regular advection found"
fi

echo "=====strong scaling 1====="
for nT in "${thread_num[@]}"; do
  tmp_arr=()
  for i in {1..3}; do
    if [ $nT == 1 ]; then
      cmd="../milestone-1/advection 3200 400 1.0 1.0e3 5.0e-7 2.85e-7"
    else
      cmd="./omp_advection 3200 400 1.0 1.0e3 5.0e-7 2.85e-7 $nT"
    fi

    rate=$($cmd | tail -1 | cut -d' ' -f3)
    tmp_arr[i]=$rate
  done
  sorted=($(printf '%s\n' "${tmp_arr[@]}" | sort))
  median=${sorted[1]}
  echo "$nT $median"
done

echo "=====strong scaling 1====="
for nT in "${thread_num[@]}"; do
  tmp_arr=()
  for i in {1..3}; do
    if [ $nT == 1 ]; then
      cmd="../milestone-1/advection 200 400 1.0 1.0e3 5.0e-7 2.85e-7"
    else
      cmd="./omp_advection 200 400 1.0 1.0e3 5.0e-7 2.85e-7 $nT"
    fi
    rate=$($cmd | tail -1 | cut -d' ' -f3)
    tmp_arr[i]=$rate
  done
  sorted=($(printf '%s\n' "${tmp_arr[@]}" | sort))
  median=${sorted[1]}
  echo "$nT $median"
done

echo "=====weak scaling benchmark====="

for nT in "${thread_num[@]}"; do
  tmp_arr=()
  for i in {1..3}; do
    N=$((800 * nT))
    if [ $nT == 1 ]; then
      cmd="../milestone-1/advection $N 400 1.0 1.0e3 5.0e-7 2.85e-7"
    else
      cmd="./omp_advection $N 400 1.0 1.0e3 5.0e-7 2.85e-7 $nT"
    fi
    rate=$($cmd | tail -1 | cut -d' ' -f3)
    tmp_arr[i]=$rate
  done
  sorted=($(printf '%s\n' "${tmp_arr[@]}" | sort))
  median=${sorted[1]}
  echo "$nT $median"
done
