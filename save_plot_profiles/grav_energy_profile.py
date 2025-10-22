import numpy as np
import yt
import pandas as pd
import unyt
from sys import argv

G = unyt.gravitational_constant_cgs.value

rin=2e8

rad_list = []
grav_list = []

f = argv[1]
ds = yt.load(f)
#r_cc  rho0  rhoh0  p0  gamma1bar tempbar
df = pd.read_csv(f"{f}/BaseCC_0", header=0, sep='\\s+')
df_fc = pd.read_csv(f"{f}/BaseFC_0", header=0, sep='\\s+')
rad = df.r_cc.to_numpy(dtype=np.float64)
rad_edge = df_fc.r_edge.to_numpy(dtype=np.float64)
rho = df.rho0.to_numpy(dtype=np.float64)
dr = rad_edge[1] - rad_edge[0]

#start summing enclosed mass and potential grav energy
rad_list.append(rad[0])
m_enc = 4./3. * np.pi * rho[0] * rad[0]**3
grav_list.append(-4 * np.pi * G * m_enc * rad[0] * rho[0]* dr)

for i in range(1, len(rad)):
    if rad[i] < rin:
        #calc mass from upper half of preceding cell to lower hafl of current cell.
        term1 = rho[i-1] * (rad_edge[i]**3 - rad[i-1]**3)
        term2 = rho[i] * (rad[i]**3 - rad_edge[i]**3)
        m_enc += 4./3. * np.pi * (term1 + term2)

        dr = rad_edge[i+1] - rad_edge[i]
        rad_list.append(rad[i])
        grav_list.append(-4 * np.pi * G * m_enc * rad[i] * rho[i] * dr)
    else:
        break
print(m_enc)
save_arr = np.array((rad_list, grav_list)).T #transpose to make it row x col
np.save(f"{ds.basename}_grav_energy_profile.npy", save_arr)
