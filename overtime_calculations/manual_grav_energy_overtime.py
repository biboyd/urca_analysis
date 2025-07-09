import numpy as np
import yt
import pandas as pd
import unyt
from os import listdir

G = unyt.gravitational_constant_cgs.value

rin=2e8

t_list = []
grav_list = []

for f in np.sort(listdir("plotfiles/")):
    if '000' in f:
        ds = yt.load(f"plotfiles/{f}")
        prof = yt.create_profile(ds.all_data(), 'radius', [('boxlib', 'rho')], n_bins=300, extrema={'radius':(0, rin)}, logs={'radius':False})
        df = prof.to_dataframe()
        
        rad = df.radius.to_numpy(dtype=np.float64)
        dr = rad[11] - rad[10]
        rad_edge = rad - dr
        rho = df.rho.to_numpy(dtype=np.float64)

        #start summing enclosed mass and potential grav energy
        m_enc = 4./3. * np.pi * rho[0] * rad[0]**3
                    
        dr = rad_edge[1] - rad_edge[0]
        tot_grav_energy = -4 * np.pi * G * m_enc * rad[0] * rho[0]* dr
            
        for i in range(1, len(rad)-1):
            if rad[i] < rin:
                #calc mass from upper half of preceding cell to lower hafl of current cell.
                term1 = rho[i-1] * (rad_edge[i]**3 - rad[i-1]**3)
                term2 = rho[i] * (rad[i]**3 - rad_edge[i]**3)
                m_enc += 4./3. * np.pi * (term1 + term2)

                dr = rad_edge[i+1] - rad_edge[i]
                tot_grav_energy -= 4 * np.pi * G * m_enc * rad[i] * rho[i] * dr
            else:
                break
        grav_list.append(tot_grav_energy)    
        t_list.append(ds.current_time)

save_arr = np.array((t_list, grav_list))
np.save("manual_grav_energy_overtime.npy", save_arr)
