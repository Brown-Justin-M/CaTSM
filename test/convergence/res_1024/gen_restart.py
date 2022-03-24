import f90nml
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

file = open("parameters", "r")
nl = f90nml.read(file)

grid = nl["grid"]
params = nl["parameters"]["params"]
io = nl["io"]["controls"]

nx = grid["x_modes"]
ny = 0
nz = grid["z_modes"]

d = nc.Dataset("infile.nc", "w", format="NETCDF3_64BIT")

xid = d.createDimension("x", 3 * nx)
yid = d.createDimension("y", max(3 * ny, 1))
zid = d.createDimension("z", 3 * nz)

dims = ("z", "y", "x")

d.createVariable("x", "f8", "x")
d.createVariable("y", "f8", "y")
d.createVariable("z", "f8", "z")

d.createVariable("RHO", "f8", dims)
d.createVariable("ETA", "f8", dims)
d.createVariable("C", "f8", dims)
d.createVariable("U", "f8", dims)
d.createVariable("V", "f8", dims)
d.createVariable("W", "f8", dims)

x = np.linspace(0, grid["x_length"], 3 * nx + 1)[:-1]
d["x"][:] = x
x = np.array(d["ETA"][:]) * 0.0 + x[np.newaxis,np.newaxis,:]
z = np.linspace(0, grid["z_length"], 3 * nz + 1)[:-1]
d["z"][:] = z
z = np.array(d["ETA"][:]) * 0.0 + z[:,np.newaxis,np.newaxis]

mask = 0.5 * np.tanh((z - 0.25*grid["z_length"]) / 0.001) - 0.5 * np.tanh((z - 0.75*grid["z_length"]) / 0.001)

# plt.plot(mask[:,0,0])
# plt.show()

d["W"][:] = 0.2 + 0.002 * mask

d["RHO"][:] = params["rho_ref"]

pres = 0.0*d["W"][:] + 1e4 + 0.85 * mask
# pres[:,:,np.logical_or(z<0.5,z>1.5)] = 1.0

gam_m1 = params["gamma"] - 1.0
rtmp = params["cv"] * gam_m1 * params["t_ref"]
# d["E"][:] = params["cv"] * params["gamma"] * np.log((pres * params["p_ref"]**(params["gamma"] - 1.0))**(1.0/params["gamma"]) / rtmp / d["rho"][:])
d["ETA"][:] = params["cv"] * params["gamma"] * np.log((pres * params["p_ref"]**gam_m1) ** (1/params["gamma"]) / rtmp / d["RHO"][:])

d["C"][:] = 0.0
d["V"][:] = 0.0
d["U"][:] = 0.0

d.close()