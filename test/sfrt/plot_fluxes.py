import numpy as np
import matplotlib.pyplot as plt
from paddi_utils.data import Diagnostic
import netCDF4 as nc

# d = np.genfromtxt("diag.csv", delimiter=",", names=True)

d = nc.Dataset("data.nc", "r")

nz = len(d["z"]) // 8

fluxes = []

for i in range(len(d["time"])):
    meanT = np.mean(d["T"][i,nz:3*nz], axis=(1,2))
    fluxes.append(np.mean(d["W"][i,nz:3*nz,:,:] * (d["T"][i,nz:3*nz,:,:] - meanT[:,np.newaxis,np.newaxis])))

plt.plot(d["time"], fluxes)

d = np.genfromtxt("/Users/justinbrown/Dropbox/Research/DDC/Compressible/sims/ICADDI/SFRT_problem_2/OUT01")

plt.plot(d[:,1], d[:,9]/7.5e6)

plt.show()

