import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("nx", type=int)
parser.add_argument("dt", type=float)

args = parser.parse_args()

nx = args.nx
nz = args.nx

params = "&grid\nx_modes = %i\nz_modes = %i\n" % (nx, nz)

params = params + """x_length = 2.0
z_length = 2.0
/

&parameters
params%max_time = 1.0
"""

params = params + "params%%start_dt = %e\nparams%%max_dt = %e\n" % (args.dt, args.dt)

params = params + """params%limit_mult_dt = 1000.0
params%diff_vel = 0.0
params%strat_t = 0.0
params%buoy_t = 0.0
params%diff_t = 0.0
params%strat_c = 0.0
params%buoy_c = 0.0
params%diff_c = 0.0
/

&io
"""

params = params + "controls%%steps_data = %i" % (int(1.0 / args.dt))

params = params + """controls%steps_restart = 100
controls%steps_diag = 1
controls%steps_profile = 100
controls%file_input = "infile.nc"
/
"""

file = open("parameters", "w")

file.write(params)

file.close()
