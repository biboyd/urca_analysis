import numpy as np
import yt
import pynucastro as pyna
import unyt
from os import listdir

t_list = []
mass_list = []

for f in np.sort(listdir("plotfiles/")):
    if '00' in f:
        ds = yt.load(f"plotfiles/{f}")
        alldata = ds.all_data()
        tot_mass = alldata.sum('mass')
            
        mass_list.append(tot_mass.value)    
        t_list.append(ds.current_time)

save_arr = np.array((t_list, mass_list))
np.save("mass_all_overtime.npy", save_arr)
