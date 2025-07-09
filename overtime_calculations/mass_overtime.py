import numpy as np
import yt
import pynucastro as pyna
import unyt
from os import listdir

rin=2.5e8

t_list = []
mass_list = []

for f in np.sort(listdir("plotfiles/")):
    if '00' in f:
        ds = yt.load(f"plotfiles/{f}")
        sph = ds.sphere(ds.domain_center, (rin, 'cm'))
        tot_mass = sph.sum('mass')
            
        mass_list.append(tot_mass.value)    
        t_list.append(ds.current_time)

save_arr = np.array((t_list, mass_list))
np.save("mass_overtime.npy", save_arr)
