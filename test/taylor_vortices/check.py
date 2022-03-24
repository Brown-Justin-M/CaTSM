import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

d = nc.Dataset("data.nc", "r")

x = d["x"][:]
x = x[np.newaxis,:] + 0.0 * d["T"][0,:,0,:]
z = d["z"][:]
z = z[:,np.newaxis] + 0.0 * d["T"][0,:,0,:]

i = -1
dv = (d["x"][1] - d["x"][0]) * (d["z"][1] - d["z"][0])

analytic_u = -np.cos(np.pi * x) * np.sin(np.pi * z)
analytic_w = np.sin(np.pi * x) * np.cos(np.pi * z)
val = np.sum(np.abs(d["U"][i,:,0,:] - analytic_u) + np.abs(d["W"][i,:,0,:] - analytic_w)) / np.sum(analytic_u == analytic_u) / 2.0

# print(np.max(analytic_u))
# print(np.max(d["U"][i,:,0,:]))

# val = np.sum(analytic_u**2 + analytic_w**2) * dv
# print(val)
# print(np.sum(d["U"][i,:,0,:]**2 + d["W"][i,:,0,:]**2) * dv)
# val = np.abs(val - np.sum(d["U"][i,:,0,:]**2 + d["W"][i,:,0,:]**2) * dv)

fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

pc1 = ax1.pcolormesh(x, z, analytic_u**2 + analytic_w**2)
cb = plt.colorbar(pc1, ax=ax1)
pc2 = ax2.pcolormesh(x, z, d["U"][i,:,0,:]**2 + d["W"][i,:,0,:]**2 - analytic_u**2 - analytic_w**2)
cb = plt.colorbar(pc2, ax=ax2)

print(len(d["x"]) // 3, d["dt"][i], val)
plt.show()
