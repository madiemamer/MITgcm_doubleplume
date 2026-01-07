from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import params_dp as p

fname = "run/"

totalIters = (p.params03['endTime']) / p.params03['deltaT']
dumpInter = p.params03['dumpFreq'] / p.params03['deltaT']
print(f"Files stored at intervals of {dumpInter}")

Nz = p.grid_params['Nr']; 
Ny = p.grid_params['Ny']; y_mp = int(Ny/2)
Nx = p.grid_params['Nx']; x_mp = int(Nx/2)

Lx_m = p.run_config['Lx_m']
Ly_m = p.run_config['Ly_m']
H = p.domain_params['H']

x = np.linspace(0,Lx_m,Nx)/1000
y = np.linspace(0,Ly_m,Ny)/1000
z = np.linspace(0,H,Nz)

i = int(totalIters - dumpInter)

dynDiag_inst = np.fromfile(fname + "dynDiag_inst." + f"{i:010d}" + ".001.001.data", dtype="<f4")  # try big-endian double
dynDiag_inst = dynDiag_inst.reshape((5,50,11,50))#,order = 'F')

T = dynDiag_inst[0,:,:,:]
S = dynDiag_inst[1,:,:,:]
U = dynDiag_inst[2,:,:,:]
W = dynDiag_inst[3,:,:,:]
V = dynDiag_inst[4,:,:,:]

plumeDiag_inst = np.fromfile(fname + "plumeDiag_inst." + f"{i:010d}" + ".001.001.data", dtype="<f4")  # try big-endian double
plumeDiag_inst = plumeDiag_inst.reshape((5,50,11,50))#,order = 'F')

Qi = plumeDiag_inst[0,:,:,:]
idx = np.where(Qi !=0)
Qi = Qi[idx]
Mi = plumeDiag_inst[1,:,:,:]
Mi = Mi[idx]
Hi = plumeDiag_inst[2,:,:,:]
Hi = Hi[idx]
Bi = plumeDiag_inst[3,:,:,:]
Bi = Bi[idx]
Ci = plumeDiag_inst[4,:,:,:]
Ci = Ci[idx]

ui = Mi/Qi 
bi = np.sqrt(2 * Qi/(np.pi * ui))
Ti = Hi / (Qi * 1031 * 4000) - 273.15 
Si = Bi / Qi 
ci = Ci / Qi 

## Plotting Plume Profiles
fig, ax = plt.subplots(1,5,figsize = (25,7))

ax[0].plot(ui, z); ax[0].set_xlabel('Velocity [m s$^{-1}$]')
ax[1].plot(bi, z); ax[1].set_xlabel('Width [m]')
ax[2].plot(Ti, z); ax[2].set_xlabel('Temperature [C]')
ax[3].plot(Si, z); ax[3].set_xlabel('Salinity [ppt]')
ax[4].plot(ci, z); ax[4].set_xlabel('$O_2$ concentration [kg m$^{-3}$]')
plt.savefig(f'plumeDiag_{i}.png', dpi = 300)
plt.show()
colormap = 'coolwarm'

fig, ax = plt.subplots(1,3,figsize = (25,5))

ax[0].contourf(x, z, T[:,y_mp,:], levels = 50, cmap = colormap)
ax[0].set_xlabel('Along Fjord [km]')
ax[0].set_ylabel('Depth [m]')
c = ax[1].contourf(y, z, T[x_mp,:,:].T, levels = 50, cmap = colormap)
ax[1].set_xlabel('Across Fjord [km]')
ax[1].set_ylabel('Depth [m]')

ax[2].contourf(x, y, T[:,:,0].T, levels = 50, cmap = colormap)
ax[2].set_xlabel('Along Fjord [km]')
ax[2].set_ylabel('Across Fjord [km]')
plt.colorbar(c)
plt.savefig(f'T_{i}.png', dpi = 300)
plt.show()

fig, ax = plt.subplots(1,3,figsize = (25,5))

ax[0].contourf(x, z, W[:,y_mp,:], levels = 50, cmap = colormap)
ax[0].set_xlabel('Along Fjord [km]')
ax[0].set_ylabel('Depth [m]')
c = ax[1].contourf(y, z, W[x_mp,:,:].T, levels = 50, cmap = colormap)
ax[1].set_xlabel('Across Fjord [km]')
ax[1].set_ylabel('Depth [m]')

ax[2].contourf(x, y, W[:,:,0].T, levels = 50, cmap = colormap)
ax[2].set_xlabel('Along Fjord [km]')
ax[2].set_ylabel('Across Fjord [km]')
plt.colorbar(c)
plt.savefig(f'W_{i}.png', dpi = 300)
plt.show()

# ---------------------------
# Plot along the fjord (x vs z at plume center y=mp)
# ---------------------------
# plt.figure(figsize=(8,4))
# plt.contourf(x, z, T_ocean[:, 5, :].T, levels=50, cmap=colormap)
# plt.gca().invert_yaxis()  # depth increases downward
# plt.xlabel("x (along fjord)")
# plt.ylabel("z (vertical level)")
# plt.title("Temperature along the fjord (plume center)")
# plt.colorbar(label="Temperature (°C)")
# plt.show()
# plt.savefig(f'T_alongFjord_{i}.png', dpi = 300)
# ---------------------------
# Plot across the fjord (y vs z at a given x, e.g., x=10)
# ---------------------------
# x_index = 10
# plt.figure(figsize=(6,4))
# plt.contourf(np.arange(Ny), np.arange(Nz), T[x_index, :, :].T, levels=50, cmap=colormap)
# plt.gca().invert_yaxis()
# plt.xlabel("y (across fjord)")
# plt.ylabel("z (vertical level)")
# plt.title(f"Temperature across the fjord at x={x_index}")
# plt.colorbar(label="Temperature (°C)")
# plt.show()

# # ---------------------------
# # Plan view (x vs y at a given z level, e.g., near surface z=5)
# # ---------------------------
# z_index = 5
# plt.figure(figsize=(10,4))
# plt.imshow(T[:, :, z_index].T, origin='lower', aspect='auto', cmap=colormap)
# plt.xlabel("x (along fjord)")
# plt.ylabel("y (across fjord)")
# plt.title(f"Plan view of temperature at z={z_index}")
# plt.colorbar(label="Temperature (°C)")
# plt.show()

