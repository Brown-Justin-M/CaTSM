import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

d = nc.Dataset("restart.nc", "r")

data = d["W"][:] + 1.0j * d["WI"][:]

res = np.fft.irfft2(data,axes=(2,1))

print(np.max(np.real(data)))

plt.pcolor(res[0].T)
plt.show()
