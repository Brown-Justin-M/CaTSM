#!/bin/bash
for i in 32 64 128
do
    for j in 5.0e-3 2.0e-3 1.0e-3 5.0e-4 2.0e-4 1.0e-4 5.0e-5 2.0e-5
    do
        # Note, to use the gen_restart_p.py tool, first install the 
        # f90nml, numpy, and netCDF4 python modules
        python gen_params.py $i $j
        python gen_restart_p.py
        ../../ROME
        python check.py >> converge.dat
    done
done