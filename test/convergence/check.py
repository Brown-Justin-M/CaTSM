import numpy as np
import netCDF4 as nc
import f90nml
from scipy.optimize import root
import argparse

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("ref")
args = parser.parse_args()

fig, ((axp, axu), (axrho, axcon)) = plt.subplots(2, 2, figsize=(6, 5))

dir = "."
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

n = len(d["z"])
t = np.max(d["time"])
x = np.linspace(0, 2.0, n + 1)[:-1]
analytic = 0.0 * x
analytic[x < 0.5] = pL
analytic[x >= 0.5] = pR

axp.plot(x, analytic, c="k", ls="--")

analytic[np.logical_and(x > vL * t + 0.5, x <= vR * t + 0.5)] = pStar

axp.plot(x, analytic, c="k")

analytic[x <= 0.5] = rhoL
analytic[x > 0.5] = rhoR

axrho.plot(x, analytic, c="k", ls="--")

analytic[np.logical_and(x <= 0.5, x > vL * t + 0.5)] = rhoStarL
analytic[np.logical_and(x > 0.5, x <= vR * t + 0.5)] = rhoStarR

axrho.plot(x, analytic, c="k", label="Analytic")

axrho.set_xlim(0, 1.0)

analytic[x <= 0.5] = uL
analytic[x > 0.5] = uR

axu.plot(x, analytic, c="k", ls="--", label="Initial")

analytic[np.logical_and(x > vL * t + 0.5, x <= vR * t + 0.5)] = uStar

axu.plot(x, analytic, c="k", label="Analytic")

dts = []
vals = []
res = []

d = nc.Dataset(dir + "/data01.nc","r")
dref = nc.Dataset(args.ref + "/data01.nc","r")

ratio = len(dref["z"]) // len(d["z"])
ref = dref["W"][-1,::ratio,0,0]

x = d["z"][:]
n = len(x)
analytic = 0.0 * x
analytic[x <= 0.5] = uL
analytic[x > 0.5] = uR
analytic[np.logical_and(x > vL * d["time"][-1] + 0.5, x <= vR * d["time"][-1] + 0.5)] = uStar

print(n // 3, params["max_dt"], np.mean(np.abs(d["W"][-1,:n//2,0,0] - analytic[:n//2])), np.sum(np.abs(d["W"][-1,:n//2,0,0] - ref[:n//2])) / (n // 2), np.abs(d["W"][-1,n//4,0,0] - analytic[n//4]), np.abs(d["W"][-1,n//4,0,0] - ref[n//4]))

