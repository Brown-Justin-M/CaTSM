import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

plt.style.use("paper")

fig, ax1 = plt.subplots(1, 1, figsize=(6.5, 4.5))
fig, ax = plt.subplots(1, 1, figsize=(6.5, 4.5))

d = nc.Dataset("res_128/data.nc", "r")
nz = len(d["z"])
rhoMean = np.trapz(d["RHO"][0,:nz//2,0,:], d["x"], axis=1)
rhoGH = rhoMean * 9.81 * d["z"][:nz//2]
rhoMean = np.sort(rhoMean)
rhoGHSort = rhoMean * 9.81 * d["z"][:nz//2]
APE = np.abs(np.trapz(rhoGH, d["z"][:nz//2]) - np.trapz(rhoGHSort, d["z"][:nz//2]))
# APE = 1.0

print(APE)

dA = (d["z"][1] - d["z"][0]) * (d["x"][1] - d["x"][0])

for i in [32, 64, 128]:
    d = np.genfromtxt("res_%i/diag.csv" % i, delimiter=",", names=True)

    n = np.min([2000, len(d["time"]) - 2]) // 2

    eoidx = np.argmax(d["time"] > 1.0)
    eidx = np.argmax(d["time"] > 50.0)
    if eidx == 0: eidx = len(d["time"]) - 1
    sidx = np.argmax(d["time"] > 20.0)
    if sidx == 0: sidx = np.max(eidx - 2 * n, 0)

    E = d["ke"] + d["pe"] + d["ie"]
    E0 = np.mean(E[eoidx:eoidx + n])
    E -= E0
    pe0 = d["pe"][10]
    ax1.plot(d["time"], E)
    # ax1.plot(d["time"], d["ke"])
    # ax1.plot(d["time"], d["ie"] - E0 + pe0)
    # ax1.plot(d["time"], d["pe"] - pe0)
    print(E[eidx], E[sidx])
    ax1.scatter([d["time"][sidx], d["time"][eidx]], [np.mean(E[sidx:sidx+n]), np.mean(E[eidx-n:eidx])], zorder=10, edgecolors="k")
    # ax1.plot(d["time"], d["ke"] + d["pe"] - pe0)
    dt = np.mean(d["time"][eidx-n:eidx]) - np.mean(d["time"][sidx:sidx+n])
    # plt.scatter([i], [np.abs(((np.mean(E[eidx-n:eidx]) - np.mean(E[sidx:sidx+n])) / APE)) / dt], c="C0")
    plt.scatter([i], [(np.abs(np.mean(E[eidx-n:eidx]) - np.mean(E[sidx:sidx+n]))) / APE], c="C0")

    if i == 128: control = np.abs((E[eidx] - E[sidx] / APE))

nx = np.array([32., 64., 128., 256., 512.])
# plt.plot(nx, (128**2 * control)*nx**-2, ls="--", label="2nd order", c="C1")
# plt.plot(nx, (128**3 * control)*nx**-3, ls="--", label="3rd order", c="C2")
# plt.plot(nx, (128**4 * control)*nx**-4, ls="--", label="4th order", c="C3")

plt.yscale("log")
plt.xscale("log")
ax.set_xticks([32, 64, 128, 256, 512])
ax.set_xticks([], minor=True)
ax.set_xticklabels(["$32$", "$64$", "$128$", "$256$", "$512$"])

plt.xlabel("$N_x$")
plt.ylabel("$\overline{\mathrm{E}'} / \mathrm{APE}$")
plt.legend()
plt.tight_layout()
plt.savefig("rt_energy.pdf")
plt.show()
