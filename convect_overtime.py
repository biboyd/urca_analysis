import yt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from os import listdir
import pandas as pd

def _magvel(field, data):
    return np.sqrt(data['boxlib', 'velx']**2 + data['boxlib', 'vely']**2 + data['boxlib', 'velz']**2)

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
magvel_out = np.zeros((len(all_files), 3), dtype=np.float64) -1.
conv_rad_out = np.zeros((len(all_files), 3), dtype=np.float64) -1.

#previous data
try:
    old_conv_time = np.load("conv_time_over_time.npy")
    old_conv_zone = np.load("conv_zone_over_time.npy")
    old_magvel = np.load("magvel_over_time.npy")
    calc_all = False
    
except FileNotFoundError:
    calc_all = True
    print("Could not find old files, will be calculating for each file")

for j, file in enumerate(all_files):
    if file[:3] == "plt" and not file[-3:] == "_nu" and not file[-5:] == "_conv":
        
        ds = yt.load(f"{plots_dir}/{file}", hint='maestro')
        sim_time = ds.current_time
        timestep = float(ds.basename.removeprefix("plt"))
        if ('boxlib', 'magvel') not in ds.field_list:
            ds.add_field(name=('boxlib', 'magvel'), function=_magvel, sampling_type='local')
        
        #ignore early times
        if sim_time < ds.quan(20, 's'):
            continue
        else:
            # check if already calculated
            if not calc_all and timestep in old_conv_time[:, 2] and timestep in old_conv_zone[:, 2] and timestep in old_magvel[:, 2]:
                idx = np.argmin(np.abs(timestep - old_conv_time[:, 2]))
                output[j, :] = old_conv_time[idx, :]
                
                idx = np.argmin(np.abs(timestep - old_conv_zone[:, 2]))
                conv_rad_out[j, :] = old_conv_zone[idx, :]
                
                idx = np.argmin(np.abs(timestep - old_magvel[:, 2]))
                magvel_out[j, :] = old_magvel[idx, :]
            
            else:
                # calc convect zone
                if rad is None:
                    try:
                        df = pd.read_csv(f"profiles/{file}.csv", index_col=0).dropna()
                    except FileNotFoundError:
                        continue
                    i = np.argmin(np.abs(df['X(c12)'].to_numpy() - (0.39974)))
                    
                    r = ds.quan(df['radius'].iloc[i], 'cm').in_units('km')
        
                else:
                    r = ds.quan(rad, 'km')
                    
                try:
                    sph = ds.sphere(ds.domain_center, r)
                except yt.utilities.exceptions.YTSphereTooSmall:
                    continue
                mean_magvel = sph.mean("magvel").in_units('km/s')
                # now take actual avg inflow /outflow speed
                inflow_conv_reg = sph.cut_region(["obj[('gas', 'radial_velocity')] < 0."])
                outflow_conv_reg = sph.cut_region(["obj[('gas', 'radial_velocity')] > 0."])
        
                print(ds.basename, r)
                inflow_mean = inflow_conv_reg.mean("radial_velocity").in_units('km/s')
                outflow_mean = outflow_conv_reg.mean("radial_velocity").in_units('km/s')
        
                timescale =  (r / outflow_mean) - (r / inflow_mean)
        
                #print(outflow_mean, inflow_mean, r, timescale, ds.basename)
                output[j, :] = [sim_time, timescale.value, timestep]
                magvel_out[j, :] = [sim_time, mean_magvel.value, float(ds.basename.removeprefix("plt"))]
                conv_rad_out[j, :] = [sim_time, r.value, float(ds.basename.removeprefix("plt"))]
        

#sort array by sim time
sorted_output = output[ np.argsort(output[:, 0]) ]
sorted_output = sorted_output[ np.where(sorted_output[:, 0] > 0) ]

#sort array by sim time
sorted_magvel = magvel_out[ np.argsort(magvel_out[:, 0]) ]
sorted_magvel = sorted_magvel[ np.where(sorted_magvel[:, 0] > 0) ]

#sort array by sim time
sorted_conv = conv_rad_out[ np.argsort(conv_rad_out[:, 0]) ]
sorted_conv = sorted_conv[ np.where(sorted_conv[:, 0] > 0) ]

#plot it all out
plt.figure(1)
plt.plot(sorted_output[:, 0], sorted_output[:, 1])

plt.ylabel("Convection Timescale (s)")
plt.xlabel("simulation time (s)")

plt.savefig("conv_time_over_time.png")
np.save("conv_time_over_time.npy", sorted_output)

#and second plot
plt.figure(2)
plt.plot(sorted_magvel[:, 0], sorted_magvel[:, 1])

plt.ylabel("avg magnitude velocity (km/s)")
plt.xlabel("simulation time (s)")

plt.savefig("magvel_over_time.png")
np.save("magvel_over_time.npy", sorted_magvel)

# and third plot
plt.figure(3)
plt.plot(sorted_conv[:, 0], sorted_conv[:, 1])

plt.ylabel("Convection Zone Radius (km)")
plt.xlabel("simulation time (s)")

plt.savefig("conv_zone_over_time.png")
np.save("conv_zone_over_time.npy", sorted_conv)
