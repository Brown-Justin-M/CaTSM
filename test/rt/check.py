import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

d = np.genfromtxt("diag.csv", delimiter=",", names=True)

pe0 = d["pe"][5]
te0 = np.mean(d["ie"][-10:] + d["pe"][-10:] + d["ke"][-10:])
plt.plot(d["time"], d["ke"], label="KE")
plt.plot(d["time"], d["pe"] - pe0, label="PE")
plt.plot(d["time"], d["ie"] - te0 + pe0, label="IE")
plt.plot(d["time"], d["ie"] + d["pe"] + d["ke"] - te0, label="TE", ls="--", c="k")

plt.legend()
plt.show()
