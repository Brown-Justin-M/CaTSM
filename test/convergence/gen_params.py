import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("nz", type=int)
parser.add_argument("dt", type=float)
args = parser.parse_args()

nz = args.nz
dt = args.dt

steps = int(0.01 / dt)

params = "&grid\nx_modes = 2\nz_modes = %i\nx_length = 1.0\nz_length = 2.0\n/\n\n" % (nz)

params = params + "&parameters\nparams%%max_time = 0.01\nparams%%start_dt = %e\nparams%%max_dt = %e\n" % (dt,dt)

params = params + """params%diff_vel = 0.0
params%diff_t = 0.0e-8
params%buoy = 0.0
params%diff_c = 0.0e-8
params%diff_vel = 0.0e-8
params%rho_ref = 25.0
params%gamma = 1.4
params%cv = 714.2857143
params%t_ref = 300.0
params%p_ref = 1.0e4
params%eos = "ideal"
params%reflection = .false.
params%implicit = .false.
/

&io
"""

params = params + "controls%%steps_data = %i\n" % (int(steps / 10))

params = params + """controls%steps_restart = 100
controls%steps_diag = 1
controls%steps_profile = 100
controls%file_input = "infile.nc"
controls%file_data = "data01.nc"
controls%file_diag = "diag01.csv"
controls%file_profile = "profiles01.nc"
controls%file_restart = "restart01.nc"
/
"""

file = open("parameters", "w")

file.write(params)

file.close()
