import f90nml
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

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

mask = (np.tanh(((z - grid["z_length"] / 4.0) - 0.05 * np.cos(2.0 * np.pi * x / grid["x_length"]) + 0.001 * np.cos(10.0 * np.pi * x / grid["x_length"])) / 0.01) / 2.0 + 0.5)
mask[len(d["z"])//2:,:,:] = mask[len(d["z"])//2-1::-1,:,:]

rho = params["rho_ref"] + 3.8613e-05 * mask
rhog = rho*params["buoy"]*2*np.arctan(np.sin(2*np.pi*z/grid["z_length"])/0.01)/np.pi
gint = integrate.cumtrapz(rhog,z,axis=0,initial=0)
pres = params["p_ref"] - gint + 5.8

# plt.plot(pres[:,0,0])
# plt.show()


temp = params["t_ref"] + 0.0035 - 0.01 * mask
rho = pres / ((params["gamma"] - 1.0) * params["cv"]) / temp
rhog = rho*params["buoy"]*2*np.arctan(np.sin(2*np.pi*z/grid["z_length"])/0.01)/np.pi
gint = integrate.cumtrapz(rhog,z,axis=0,initial=0)
pres = params["p_ref"] - gint + 5.8
rho = pres / ((params["gamma"] - 1.0) * params["cv"]) / temp
print(np.max(pres), np.min(pres))
en = params["cv"] * params["gamma"] * np.log(temp / params["T_ref"] / (pres / params["p_ref"]) ** ((params["gamma"] - 1) / params["gamma"]))

mask = np.cos(2*np.pi*x / grid["x_length"]) * np.exp(-(z - grid["z_length"]/4.0)**2/0.1**2)
mask[len(d["z"])//2:,:,:] = -mask[len(d["z"])//2-1::-1,:,:]

d["RHO"][:] = rho
d["ETA"][:] = en
d["C"][:] = rho

d["U"][:] = 0.0
d["V"][:] = 0.0
d["W"][:] = 0.0e-2*mask

d.close()