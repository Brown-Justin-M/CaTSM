#!/bin/bash
for i in 128 256
do
    for j in 5.0e-6 2.0e-6 1.0e-6 5.0e-7 2.0e-7
    do
        python gen_params.py $i $j
        python gen_restart.py
        ../../CaTSM
        python check.py res_$i >> converge.dat
    done
done
