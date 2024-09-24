import yt
import unyt
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import pandas as pd


# fields to add and stuff.
def _invert_ecap_losses(field, data):
    return -1.0 * data[ 'A23_electron_capture_nu_loss'] * unyt.erg/unyt.g/unyt.s

def _invert_beta_losses(field, data):
    return -1.0 * data[ 'A23_beta_decay_nu_loss'] * unyt.erg/unyt.g/unyt.s

def _tot_nu_loss(field, data):
    return  (data['thermal_nu_loss'] -1.0 * data[ 'A23_electron_capture_nu_loss']  -1.0 * data[ 'A23_beta_decay_nu_loss'] )* unyt.erg/unyt.g/unyt.s


def _tot_energy_loss(field, data):
    return data[('boxlib', 'tot_nu_loss')]*data[('gas', 'mass')]

def _mass(field, data):
    return data[('boxlib', 'density')] * data[('boxlib', 'volume')]



# okay grab all the files I can
yt.set_log_level(40)

if len(argv) == 2:
    if argv[1] == "sw":
        rad = 378.4 #default sw radius
    elif argv[1] == "my" or argv[1] == "myfix":
        rad = None#482.4 #default myfix radius        
    else:
        rad = float(argv[1])
        
else:
    rad = None#482.4 #default myfix radius
    
plots_dir="plotfiles"
all_files = listdir(plots_dir)
output = np.zeros((len(all_files), 3), dtype=np.float64) -1.

#previous data
try:
    old_energy_time = np.load("nu_loss_over_time.npy")
    check_old = True
except FileNotFoundError:
    check_old = False
    print("Could not find old file, will be calculating for each plot-file")

for j, file in enumerate(all_files):
    if file[:3] == "plt" and file[-3:] == "_nu":
        
        ds = yt.load(f"{plots_dir}/{file}")
                 
        # add fields we need.
        ds.add_field(name=('gas', 'mass'), 
                        function=_mass, take_log=True, 
                        dimensions='mass', sampling_type='local', 
                        force_override=True)
        
        
        ds.add_field(
            name=("boxlib", "nu_ecap_losses"),
            function=_invert_ecap_losses,
            take_log=True,
            units='erg/g/s',
            display_name="$\\nu_{ecap}$ energy losses",
            sampling_type="local", force_override=True)
        
        ds.add_field(
            name=("boxlib", "nu_beta_losses"),
            function=_invert_beta_losses,
            take_log=True,
            units='erg/g/s',
            display_name="$\\nu_{beta}$ energy losses",
            sampling_type="local", force_override=True)
        
        ds.add_field(
            name=("boxlib", "tot_nu_loss"),
            function=_tot_nu_loss,
            take_log=True,
            units='erg/g/s',
            display_name="Energy Loss to Neutrino Emissions",
            sampling_type="local", force_override=True)
        
        
        ds.add_field(
            name=("boxlib", "tot_energy_loss"),
            function=_tot_energy_loss,
            take_log=True,
            units='erg/s',
            display_name="Total Power Loss to Neutrino Emissions",
            sampling_type="local", force_override=True)

        sim_time = ds.current_time
        timestep = float(ds.basename.removeprefix("plt").removesuffix("_nu"))
        
        #ignore early times
        if sim_time < ds.quan(200, 's'):
            continue
        else:
            # check if already calculated
            if check_old and timestep in old_energy_time[:, 2]:
                idx = np.argmin(np.abs(timestep - old_energy_time[:, 2]))
                output[j, :] = old_energy_time[idx, :]
                
            else:
                sph = ds.sphere(ds.domain_center, (800, 'km'))
                tot_energy_loss = sph.sum(('boxlib', 'tot_energy_loss'))
                #print(outflow_mean, inflow_mean, r, timescale, ds.basename)
                output[j, :] = [sim_time, tot_energy_loss, timestep]
                
#sort array by sim time
sorted_output = output[ np.argsort(output[:, 0]) ]
sorted_output = sorted_output[ np.where(sorted_output[:, 0] > 0) ]

#plot it all out
plt.figure(1)
plt.plot(sorted_output[:, 0], sorted_output[:, 1])

plt.ylabel("Energy Loss to Neutrino Emission (erg/s)")
plt.xlabel("Simulation Time (s)")

plt.savefig("nu_loss_over_time.png")
np.save("nu_loss_over_time.npy", sorted_output)