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

mask = np.exp(-((z - grid["z_length"] / 4.0)**2 + (x - grid["x_length"] / 2.0)**2) / 0.01**2)
mask[len(d["z"])//2:,:,:] = mask[len(d["z"])//2-1::-1,:,:]

rho0 = params["rho_ref"]
alpha = params["beta_t"]
beta = params["beta_c"]
rho = 0.001 * mask
rho = rho - np.mean(rho) + params["rho_ref"]
rhog = rho*params["buoy"]*2*np.arctan(np.sin(2*np.pi*z/grid["z_length"])/0.01)/np.pi
gint = integrate.cumtrapz(rhog,z,axis=0,initial=0)
pres = params["p_ref"] - gint

salt = 0.001 * mask / beta / params["rho_ref"]

# temp = -0.001 * mask / alpha / params["rho_ref"]
temp = -0.0 * mask / alpha / params["rho_ref"]

en = 1.0 / rho0 / params["alr"] * (alpha * temp - beta * salt)

# plt.plot(temp[:,0,64],z[:,0,64])
# plt.show()

rho = 1 / (-pres / (params["rho_ref"]**2 * params["cs"]**2) + 1.0 / params["rho_ref"] + params["alr"] * en)

print(np.max(np.max(rho, axis=(1,2)) - np.min(rho, axis=(1,2))))

d["RHO"][:] = rho
d["ETA"][:] = en
d["C"][:] = salt

d["U"][:] = 0.0
d["V"][:] = 0.0
d["W"][:] = 0.0

d.close()
