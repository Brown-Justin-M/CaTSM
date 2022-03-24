import numpy as np
import matplotlib.pyplot as plt
from fluidfoam import readmesh, readvector, readscalar
import netCDF4 as nc

plt.style.use("paper")

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(6.5, 4.0))

time = 50

d = nc.Dataset("res_128/data.nc", "r")

i = np.argmax(d["time"][:] > time - 0.1)
print(d["time"][i])

nz = len(d["z"])//2
rhoMean = np.mean(d["T"][i,:nz,0,:] - 300)
rhoDiff = (np.max(d["T"][i,:nz,0,:]) - np.min(d["T"][i,:nz,0,:])) / 4.0
print(np.mean(d["T"][i,:nz,0,:]), np.min(d["T"][i,:nz,0,:]))
pc = ax3.pcolor(d["x"], d["z"][:nz], d["T"][i,:nz,0,:] - 300, vmin=rhoMean - rhoDiff, vmax=rhoMean + rhoDiff)
pc.set_rasterized(True)
print(np.max(d["T"][i]), np.min(d["T"][i]))
ax3.contour(d["x"][:], d["z"][:nz], d["T"][i,:nz,0,:], levels=[300.0], colors="w")
ax3.set_xlim(0, 0.25)
ax3.set_ylim(0, 1.0)


d = nc.Dataset("res_32/data.nc", "r")

i = np.argmax(d["time"][:] > time - 0.1)
print(d["time"][i])

if (i == 0): i = -1

nz = len(d["z"])//2
# rhoMean = np.mean(d["Temp"][i,:nz,0,:]) - 300
# rhoDiff = (np.max(d["Temp"][i,:nz,0,:]) - np.min(d["Temp"][i,:nz,0,:])) / 4.0
print(np.mean(d["T"][i,:nz,0,:]), np.min(d["T"][i,:nz,0,:]))
pc = ax2.pcolor(d["x"], d["z"][:nz], d["T"][i,:nz,0,:] - 300, vmin=rhoMean - rhoDiff, vmax=rhoMean + rhoDiff)
pc.set_rasterized(True)
ax2.contour(d["x"][:], d["z"][:nz], d["T"][i,:nz,0,:], levels=[300.0], colors="w")
ax2.set_xlim(0, 0.25)
ax2.set_ylim(0, 1.0)




x, y, z = readmesh("/Users/justinbrown/Dropbox/Research/DDC/Compressible/data/RT_validation/CRT_instability", False)
T = readscalar("/Users/justinbrown/Dropbox/Research/DDC/Compressible/data/RT_validation/CRT_instability", str(time), "T", False)
p = readscalar("/Users/justinbrown/Dropbox/Research/DDC/Compressible/data/RT_validation/CRT_instability", str(time), "p", False)
rho = p / T / (8.316382 / 0.029107337)
print(np.mean(T), np.max(T))

# rhoMean = np.mean(T)
# rhoDiff = (np.max(T) - np.min(T)) / 4.0
print(np.max(T), np.min(T))
pc = ax1.tripcolor(x, y + 0.5, T - 300, vmin=rhoMean - rhoDiff, vmax=rhoMean + rhoDiff)
pc.set_rasterized(True)
ax1.tricontour(x, y + 0.5, T, [300.0], zorder=100, colors=["w"])

ax1.set_xlabel("$x$ [m]")
ax2.set_xlabel("$x$ [m]")
ax3.set_xlabel("$x$ [m]")
ax1.set_ylabel("$z$ [m]")

ax1.text(0.01, 0.02, "OpenFOAM")
ax2.text(0.01, 0.02, "Spectral ($N_x=32$)")
ax3.text(0.01, 0.02, "Spectral ($N_x=128$)")

plt.tight_layout()
cb = fig.colorbar(pc, ax=[ax1, ax2, ax3])
cb.set_label("$T - T_0$ [K]")

plt.savefig("rt_compare.pdf", dpi=300)
plt.show()
