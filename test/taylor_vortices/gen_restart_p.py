import f90nml
import netCDF4 as nc
import numpy as np

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

d.createVariable("ETA", "f8", dims)
d.createVariable("RHO", "f8", dims)
d.createVariable("C", "f8", dims)
d.createVariable("U", "f8", dims)
d.createVariable("V", "f8", dims)
d.createVariable("W", "f8", dims)

x = np.linspace(0, grid["x_length"], 3 * nx + 1)[:-1]
x = np.array(d["ETA"][:]) * 0.0 + x[np.newaxis,np.newaxis,:]
z = np.linspace(0, grid["z_length"], 3 * nz + 1)[:-1]
z = np.array(d["ETA"][:]) * 0.0 + z[:,np.newaxis,np.newaxis]

d["ETA"][:] = 0.0
d["RHO"][:] = params["rho_ref"]
d["C"][:] = 0.0
d["U"][:] = -np.cos(np.pi * x) * np.sin(np.pi * z)
d["V"][:] = 0.0
d["W"][:] = np.sin(np.pi * x) * np.cos(np.pi * z)

d.close()
