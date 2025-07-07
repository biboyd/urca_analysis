import numpy as np
import yt
import pynucastro as pyna
import unyt
from os import listdir

G = unyt.gravitational_constant_cgs

ds = yt.load("plotfiles/plt0008400/")

@yt.derived_field(name=("gas", "mass"), units='g', sampling_type='local')
def _mass(field, data):
    # sum 1/2
    mass =data['boxlib', 'rho']*data['index', 'cell_volume'] * unyt.g/unyt.cm**3
    return mass

@yt.derived_field(name=("gas", "tot_internal_energy"), sampling_type='local', units='erg')
def _energy(field, data):
    return data['boxlib', 'internal_energy']*data['gas', 'mass'] * unyt.erg/unyt.g

rconv=7e7

t_list = []
int_list = []
mass_list = []

for f in np.sort(listdir("nuloss_plotfiles/")):
    if '00' in f:
        ds = yt.load(f"nuloss_plotfiles/{f}")
        sph = ds.sphere(ds.domain_center, (rconv, 'cm'))
        tot_int_energy = sph.sum('tot_internal_energy')
            
        int_list.append(tot_int_energy.value)    
        mass_list.append(sph.sum('mass').value)    
        t_list.append(ds.current_time)

save_arr = np.array((t_list, int_list, mass_list))
np.save("int_conv_energy_overtime.npy", save_arr)
