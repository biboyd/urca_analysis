import numpy as np
import yt
import pandas as pd
import unyt
from os import listdir

G = unyt.gravitational_constant_cgs.value

rconv=7e7 #cm 

t_list = []
grav_list = []
mass_list = []

for f in np.sort(listdir("plotfiles/")):
    if '00' in f:
        ds = yt.load(f"plotfiles/{f}")
        #r_cc  rho0  rhoh0  p0  gamma1bar tempbar
        df = pd.read_csv(f"plotfiles/{f}/BaseCC_0", header=0, sep='\\s+')
        df_fc = pd.read_csv(f"plotfiles/{f}/BaseFC_0", header=0, sep='\\s+')
        rad = df.r_cc.to_numpy(dtype=np.float64)
        rad_edge = df_fc.r_edge.to_numpy(dtype=np.float64)
        rho = df.rho0.to_numpy(dtype=np.float64)
        dr = rad_edge[1] - rad_edge[0]

        #start summing enclosed mass and potential grav energy
        m_enc = 4./3. * np.pi * rho[0] * rad[0]**3
        tot_grav_energy = -4 * np.pi * G * m_enc * rad[0] * rho[0]* dr
            
        for i in range(1, len(rad)):
            if rad[i] < rconv:
                #calc mass from upper half of preceding cell to lower hafl of current cell.
                term1 = rho[i-1] * (rad_edge[i]**3 - rad[i-1]**3)
                term2 = rho[i] * (rad[i]**3 - rad_edge[i]**3)
                m_enc += 4./3. * np.pi * (term1 + term2)

                dr = rad_edge[i+1] - rad_edge[i]
                tot_grav_energy -= 4 * np.pi * G * m_enc * rad[i] * rho[i] * dr
            else:
                break
        grav_list.append(tot_grav_energy)    
        mass_list.append(m_enc)    
        t_list.append(ds.current_time)

save_arr = np.array((t_list, grav_list, mass_list))
np.save("grav_conv_energy_overtime.npy", save_arr)
