from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import params_dp as p

fname = "run/"

timeSteps = (p.params03['endTime']) / p.params03['deltaT']
dumpInter = p.params03['dumpFreq'] / p.params03['deltaT']
print(f"Files stored at intervals of {dumpInter}")

Nz = p.grid_params['Nr']
Ny = p.grid_params['Ny']
Nx = p.grid_params['Nx']

Lx_m = p.run_config['Lx_m']
Ly_m = p.run_config['Ly_m']
H_m = p.domain_params['H']

x = np.linspace(0,Lx_m,Nx)
y = np.linspace(0,Ly_m,Ny)
z = np.linspace(0,H,Nz)

i = dumpIter * 2
T_ocean = np.fromfile(fname + "T." + f"{i:010d}" + ".001.001.data", dtype="<f4")  # try big-endian double
T_ocean = T_ocean.reshape((Nx,Ny,Nz), order = 'F')

olormap = 'coolwarm'
# ---------------------------
# Plot along the fjord (x vs z at plume center y=mp)
# ---------------------------
plt.figure(figsize=(8,4))
plt.contourf(np.arange(Nx), np.arange(Nz), T[:, mp, :].T, levels=50, cmap=colormap)
plt.gca().invert_yaxis()  # depth increases downward
plt.xlabel("x (along fjord)")
plt.ylabel("z (vertical level)")
plt.title("Temperature along the fjord (plume center)")
plt.colorbar(label="Temperature (°C)")
plt.show()
plt.savefig(f'T_alongFjord_{i}.png', dpi = 300)
# ---------------------------
# Plot across the fjord (y vs z at a given x, e.g., x=10)
# ---------------------------
x_index = 10
plt.figure(figsize=(6,4))
plt.contourf(np.arange(Ny), np.arange(Nz), T[x_index, :, :].T, levels=50, cmap=colormap)
plt.gca().invert_yaxis()
plt.xlabel("y (across fjord)")
plt.ylabel("z (vertical level)")
plt.title(f"Temperature across the fjord at x={x_index}")
plt.colorbar(label="Temperature (°C)")
plt.show()

# ---------------------------
# Plan view (x vs y at a given z level, e.g., near surface z=5)
# ---------------------------
z_index = 5
plt.figure(figsize=(10,4))
plt.imshow(T[:, :, z_index].T, origin='lower', aspect='auto', cmap=colormap)
plt.xlabel("x (along fjord)")
plt.ylabel("y (across fjord)")
plt.title(f"Plan view of temperature at z={z_index}")
plt.colorbar(label="Temperature (°C)")
plt.show()

