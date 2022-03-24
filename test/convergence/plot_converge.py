import numpy as np
import netCDF4 as nc
import f90nml
from scipy.optimize import root

import matplotlib.pyplot as plt

fig, ((axp, axu), (axrho, axcon)) = plt.subplots(2, 2, figsize=(6, 5))

dir = "res_1024"
file = open(dir + "/parameters", "r")
nl = f90nml.read(file)

grid = nl["grid"]
params = nl["parameters"]["params"]

d = nc.Dataset(dir + "/data01.nc","r")

pL = 1.0e4
pR = 1.0e4 + 0.85
uL = 0.2
uR = 0.202
rhoL = 25
rhoR = 25

gamma = 1.4

aL = np.sqrt(gamma * pL / rhoL)
aR = np.sqrt(gamma * pR / rhoR)
AL = 2 / rhoL / (gamma + 1)
BL = pL * (gamma - 1) / (gamma + 1)
AR = 2 / rhoR / (gamma + 1)
BR = pR * (gamma - 1) / (gamma + 1)

def pStarEqn(pStar):
    return -(pStar - pL) * np.sqrt(AL / (pStar + BL)) - (pStar - pR) * np.sqrt(AR / (pStar + BR)) + uL - uR

r = root(pStarEqn, pL)
pStar = r.x[0]
uStar = uL - (pStar - pL) * np.sqrt(AL / (pStar + BL))
rhoStarL = rhoL * (pL * (gamma - 1) + pStar * (gamma + 1)) / (pStar * (gamma - 1) + pL * (gamma + 1))
rhoStarR = rhoR * (pR * (gamma - 1) + pStar * (gamma + 1)) / (pStar * (gamma - 1) + pR * (gamma + 1))

vL = uL - aL * np.sqrt((gamma + 1) * pStar / 2 / pL / gamma + (gamma - 1) / 2 / gamma)
vR = uR + aR * np.sqrt((gamma + 1) * pStar / 2 / pR / gamma + (gamma - 1) / 2 / gamma)

for dir in ["res_256", "res_1024"]:
    d = nc.Dataset(dir + "/data01.nc","r")
    n = len(d["z"])

    label = "$N_z=%i$" % (n//3)
    if "cn" in file:
        label = label + ", RK2CN"
    axp.plot(d["z"][:], d["P"][-1,:,0,0])
    axu.plot(d["z"][:], d["W"][-1,:,0,0], label=label)

    axrho.plot(d["z"][:], d["RHO"][-1,:,0,0], label=label)

n = len(d["z"])
t = np.max(d["time"])
x = np.linspace(0, 2.0, n + 1)[:-1]
analytic = 0.0 * x
analytic[x < 0.5] = pL
analytic[x >= 0.5] = pR

axp.plot(x, analytic, c="k", ls="--")

analytic[np.logical_and(x > vL * t + 0.5, x <= vR * t + 0.5)] = pStar

axp.plot(x, analytic, c="k", ls=":")

analytic[x <= 0.5] = rhoL
analytic[x > 0.5] = rhoR

axrho.plot(x, analytic, c="k", ls="--")

analytic[np.logical_and(x <= 0.5, x > vL * t + 0.5)] = rhoStarL
analytic[np.logical_and(x > 0.5, x <= vR * t + 0.5)] = rhoStarR

axrho.plot(x, analytic, c="k", ls=":", label="Analytic")

axrho.set_xlim(0, 1.0)

analytic[x <= 0.5] = uL
analytic[x > 0.5] = uR

axu.plot(x, analytic, c="k", ls="--", label="Initial")

analytic[np.logical_and(x > vL * t + 0.5, x <= vR * t + 0.5)] = uStar

axu.plot(x, analytic, c="k", ls=":", label="Analytic")

d1 = np.genfromtxt("converge.dat")

axu.legend()

axp.set_xlabel("$z$")
axu.set_xlabel("$z$")
axrho.set_xlabel("$z$")

axp.set_ylabel("$p$")
axu.set_ylabel("$u$")
axrho.set_ylabel(r"$\rho$")

axp.set_xlim(0,1.0)
axu.set_xlim(0,1.0)
axrho.set_xlim(0,1.0)

# p = np.polyfit(np.log(dts), np.log(vals), 1)
# print(p)

# dts = np.array(dts)
# vals = np.array(vals)
# res = np.array(res)s
d = d1[np.logical_and(d1[:,1] <= 2.0e-6, d1[:,0] == 128)]
# axcon.scatter(d[:,1], d[:,3], c="C0", label="N_z=256")
d = d1[np.logical_and(d1[:,1] <= 2.0e-6, d1[:,0] == 256)]
axcon.scatter(d[:,1], d[:,3], c="C0", label="$N_z=256$")
dts = np.array([2.0e-7, 2.0e-6])
# # axcon.scatter(dts[0], vals[0], c="C0")
# # axcon.scatter(dts[1], vals[1], c="C1")
# # axcon.scatter(dts[2], vals[2], c="C2")
# dts.sort()
# fit = np.exp(np.polyval([1,0], np.log(dts)))
# fit /= np.max(fit)
# fit *= np.max(d[:,3])
# axcon.plot(dts, fit, c="g", ls="--", label="1st order")
# fit = np.exp(np.polyval([2,0], np.log(dts)))
# fit /= np.max(fit)
# fit *= np.max(d[:,3])
# axcon.plot(dts, fit, c="b", ls="--", label="2nd order")\
fit = np.exp(np.polyval([3,0], np.log(dts)))
fit /= np.max(fit)
fit *= np.max(d[:,3])
axcon.plot(dts, fit, c="r", ls="--", label="3rd order")
axcon.set_xscale("log")
axcon.set_yscale("log")
axcon.set_xlabel("$dt$")
axcon.set_ylabel("$l_1$")
axcon.set_ylim(1e-11, None)
# axcon.set_ylim(np.min(np.abs(vals[vals>0]))*0.9, np.max(vals)*1.1)
axcon.legend()
plt.tight_layout()
plt.savefig("converge.pdf")
plt.show()
