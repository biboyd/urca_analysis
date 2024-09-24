import yt
import unyt as u
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from os import listdir
import pandas as pd

yt.enable_parallelism()

def _Ye(field, data):
    # sum 1/2
    Ye=0.5*(data['boxlib', 'X(c12)']+data['boxlib', 'X(o16)']+data['boxlib', 'X(he4)']+data['boxlib', 'X(ne20)'])
    
    #sum ones
    Ye+=(data['boxlib', 'X(h1)'])
    
    #A=23 nuc
    Ye+= (10.*data['boxlib', 'X(ne23)'] + 11.*data['boxlib', 'X(na23)'] + 12.*data['boxlib', 'X(mg23)'])/23.
    
    return Ye

def _mass(field, data):
    return data[('boxlib', 'rho')] * data[('gas', 'volume')]  * u.g / u.cm**3

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
output = np.zeros((len(all_files), 6), dtype=np.float64) -1.

#previous data
try:
    old_Ye_time = np.load("Ye_over_time.npy")
    calc_all = False
    
except FileNotFoundError:
    calc_all = True
    print("Could not find old files, will be calculating for each file")

for j, file in enumerate(all_files):
    if file[:3] == "plt" and not file[-3:] == "_nu":
        
        print(file)
        ds = yt.load(f"{plots_dir}/{file}", hint='maestro')
        sim_time = ds.current_time
        timestep = float(ds.basename.removeprefix("plt"))
        
        ds.add_field(name=("boxlib", "Ye"),
                    function=_Ye,
                    take_log=False,
                    units = "dimensionless",
                    sampling_type="local")
        
        ds.add_field(name=("gas", "mass"),
                    function=_mass,
                    take_log=True,
                    units = "g",
                    sampling_type="local")
        #ignore early times
        if sim_time < ds.quan(20, 's'):
            continue
        else:
            # check if already calculated
            if not calc_all and timestep in old_Ye_time[:, -1]:
                idx = np.argmin(np.abs(timestep - old_Ye_time[:, -1]))
                output[j, :] = old_Ye_time[idx, :]
            
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
                    conv_sph = ds.sphere(ds.domain_center, r)
                except yt.utilities.exceptions.YTSphereTooSmall:
                    continue
                    
                sph800 = ds.sphere(ds.domain_center, (800, 'km'))    
                sph1000 = ds.sphere(ds.domain_center, (1000, 'km'))    
                sph2000 = ds.sphere(ds.domain_center, (2000, 'km'))    
                
                conv_Ye = conv_sph.mean("Ye", weight=('gas', 'mass'))
                Ye800 = sph800.mean("Ye", weight=('gas', 'mass'))
                Ye1000 = sph1000.mean("Ye", weight=('gas', 'mass'))
                Ye2000 = sph2000.mean("Ye", weight=('gas', 'mass'))
        
                #print(outflow_mean, inflow_mean, r, timescale, ds.basename)
                output[j, :] = [sim_time, conv_Ye, Ye800, Ye1000, Ye2000, timestep]
        

#sort array by sim time
sorted_output = output[ np.argsort(output[:, 0]) ]
sorted_output = sorted_output[ np.where(sorted_output[:, 0] > 0) ]

#plot it all out
plt.figure(1)
plt.plot(sorted_output[:, 0], (sorted_output[:, 1] -0.5)/1e-5, label="Ye in Conv zone")
plt.plot(sorted_output[:, 0], (sorted_output[:, 2] -0.5)/1e-5, label="Ye in 800 km")
plt.plot(sorted_output[:, 0], (sorted_output[:, 3] -0.5)/1e-5, label="Ye in 1000 km")
plt.plot(sorted_output[:, 0], (sorted_output[:, 4] -0.5)/1e-5, label="Ye in 2000 km")

plt.ylabel("Electron Fraction (Ye - 0.5)/1e-5")
plt.xlabel("simulation time (s)")
plt.legend()
plt.savefig("Ye_over_time.png")
np.save("Ye_over_time.npy", sorted_output)