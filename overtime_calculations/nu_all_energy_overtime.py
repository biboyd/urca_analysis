import yt
import unyt
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from os import listdir

full_nu_list = [('boxlib', 'A21_beta_decay_nu_loss'), ('boxlib', 'A21_electron_capture_nu_loss'),
                ('boxlib', 'A23_beta_decay_nu_loss'), ('boxlib', 'A23_electron_capture_nu_loss'),
                ('boxlib', 'A25_beta_decay_nu_loss'), ('boxlib', 'A25_electron_capture_nu_loss')]

# fields to add and stuff.
def _tot_nu_loss(field, data):
    tot=0. 
    for nuloss in curr_nu_list:
        tot -= data[nuloss]
    tot+=data['thermal_nu_loss']
    tot *= unyt.erg/unyt.g/unyt.s
    return  tot*data[('gas', 'mass')]

def _mass(field, data):
    return data[('boxlib', 'rho')] * data[('boxlib', 'volume')] * unyt.g/unyt.cm**3



# okay grab all the files I can
yt.set_log_level(40)

if len(argv) == 2:
        rad = float(argv[1])
        
else:
    rad = 2e3 # km
    
plots_dir="nuloss_plotfiles/"
all_files = listdir(plots_dir)
output = np.zeros((len(all_files), 3), dtype=np.float64) -1.

#previous data
try:
    old_energy_time = np.load("nu_loss_overtime.npy")
    check_old = True
except FileNotFoundError:
    check_old = False
    print("Could not find old file, will be calculating for each plot-file")

for j, file in enumerate(np.sort(all_files)):
    if file:
        
        ds = yt.load(f"{plots_dir}/{file}")
                 
        # add fields we need.
        ds.add_field(name=('gas', 'mass'), 
                        function=_mass, take_log=True, 
                        dimensions='mass', sampling_type='local', 
                        force_override=True)
        
        curr_nu_list = []
        for fld in full_nu_list:
            if fld in ds.field_list:
                curr_nu_list.append(fld)

        ds.add_field(
            name=("boxlib", "tot_nu_loss"),
            function=_tot_nu_loss,
            take_log=True,
            units='erg/s',
            display_name="Energy Loss to Neutrino Emissions",
            sampling_type="local", force_override=True)
        
        sim_time = ds.current_time
        timestep = float(ds.basename.removeprefix("nu_loss.plt"))
        
        # check if already calculated
        if check_old and timestep in old_energy_time[:, 2]:
           idx = np.argmin(np.abs(timestep - old_energy_time[:, 2]))
           output[j, :] = old_energy_time[idx, :]
                
        else:
           sph = ds.sphere(ds.domain_center, (rad, 'km'))
           tot_energy_loss = sph.sum(('boxlib', 'tot_nu_loss'))
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

plt.savefig("nu_loss_overtime.png")
np.save("nu_loss_overtime.npy", sorted_output)
