import yt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from os import listdir
import pandas as pd

def calc_radius(ds, file):
    # calc convect zone
    try:
        df = pd.read_csv(f"profiles/{file}.csv", index_col=0).dropna()
    except FileNotFoundError:
        print(f"no profile found for {file}")
        return -1

    ## construct c12 grad using center diff.
    c12_arr = df['X(c12)'].to_numpy()
    rad_arr = df['radius'].to_numpy()

    c12_grad_arr = (c12_arr[1:] - c12_arr[:-1])/(rad_arr[1:] - rad_arr[:-1])
    avg_rad_arr = (rad_arr[1:] + rad_arr[:-1])/2.
    i = np.argmax(c12_grad_arr) 

    r = ds.quan(avg_rad_arr[i], 'cm').in_units('km')
    return r
   

def calc_magvel(ds, r):
    try:
        sph = ds.sphere(ds.domain_center, r)
    except yt.utilities.exceptions.YTSphereTooSmall:
        return -1
    
    mean_magvel = sph.mean("velocity_magnitude").in_units('km/s')
    return mean_magvel.value
    
def calc_vrms(ds, r):
    try:
        sph = ds.sphere(ds.domain_center, r)
    except yt.utilities.exceptions.YTSphereTooSmall:
        return -1
    
    vrms = np.sqrt(sph.std("velx").in_units('km/s')**2 + 
                   sph.std("vely").in_units('km/s')**2 + 
                   sph.std("velz").in_units('km/s')**2)
    return vrms.value
    

def calc_mass(ds, r):
    try:
        sph = ds.sphere(ds.domain_center, r)
    except yt.utilities.exceptions.YTSphereTooSmall:
        return -1

    mass = sph.sum("mass").in_units('Msun')
    return mass.value

# okay grab all the files I can
yt.set_log_level(40)

outdir = './'
if len(argv) == 1:
    plots_dir = "plotfiles"
else:
    plots_dir = argv[1]
    if len(argv) == 3:
        outdir = argv[2]

all_files = listdir(plots_dir)
conv_time_out = np.zeros((len(all_files), 3), dtype=np.float64) -1.
conv_zone_out = np.zeros((len(all_files), 3), dtype=np.float64) -1.
magvel_out = np.zeros((len(all_files), 3), dtype=np.float64) -1.
vrms_out = np.zeros((len(all_files), 3), dtype=np.float64) -1.
conv_mass_out = np.zeros((len(all_files), 3), dtype=np.float64) -1.
new_arrs = [conv_time_out, conv_zone_out, magvel_out, vrms_out, conv_mass_out]

#previous data
try:
    old_conv_time = np.load("conv_time_over_time.npy")
    old_conv_zone = np.load("conv_zone_over_time.npy")
    old_magvel = np.load("magvel_over_time.npy")
    old_vrms = np.load("vrms_over_time.npy")
    old_conv_mass = np.load("conv_mass_over_time.npy")
    calc_all = False
    
    old_arrs = [old_conv_time, old_conv_zone, old_magvel, old_vrms, old_conv_mass]
    
except FileNotFoundError:
    calc_all = True
    print("Could not find old files, will be calculating for each file")
    
quantity_arr = ['time', 'radius', 'magvel', 'vrms', 'mass']

for j, file in enumerate(all_files):
    if file[:3] == "plt" and len(file) == 10:
        
        ds = yt.load(f"{plots_dir}/{file}", hint='maestro')
        sim_time = ds.current_time
        timestep = float(ds.basename.removeprefix("plt"))
        
        #ignore early times
        if sim_time < ds.quan(20, 's'):
            continue
        else:
            # check if already calculated
            if not calc_all:
                for old_arr, new_arr, quantity in zip(old_arrs, new_arrs, quantity_arr):
                    if timestep in old_arr:
                        idx = np.argmin(np.abs(timestep - old_arr[:, 2]))
                        new_arr[j, :] = old_arr[idx, :]
                        
                    elif quantity == 'time':
                        new_arr[j, :] = [sim_time, timescale.value, timestep]
                        
                    elif quantity == 'radius':
                        new_arr[j, :] = [sim_time, calc_radius(ds, file), timestep]
                        
                    elif quantity == 'magvel':
                        r = conv_zone_out[j, 1]
                        new_arr[j, :] = [sim_time, calc_magvel(ds, r), timestep]
                        
                    elif quantity == 'vrms':
                        new_arr[j, :] = [sim_time, calc_vrms(ds, r), timestep]
                        
                    elif quantity == 'mass':
                        new_arr[j, :] = [sim_time, calc_mass(ds, r), timestep]
                        
                        
            else:
                # calc convect zone
                try:
                    df = pd.read_csv(f"profiles/{file}.csv", index_col=0).dropna()
                except FileNotFoundError:
                    print(f"no profile found for {file}")
                    continue

                ## construct c12 grad using center diff.
                c12_arr = df['X(c12)'].to_numpy()
                rad_arr = df['radius'].to_numpy()

                c12_grad_arr = (c12_arr[1:] - c12_arr[:-1])/(rad_arr[1:] - rad_arr[:-1])
                avg_rad_arr = (rad_arr[1:] + rad_arr[:-1])/2.

                i = np.argmax(c12_grad_arr) 

                r = ds.quan(avg_rad_arr[i], 'cm').in_units('km')

                    
                try:
                    sph = ds.sphere(ds.domain_center, r)
                except yt.utilities.exceptions.YTSphereTooSmall:
                    continue
                mean_magvel = sph.mean("velocity_magnitude").in_units('km/s')
                vrms = np.sqrt(sph.std("velx").in_units('km/s')**2 + 
                              sph.std("vely").in_units('km/s')**2 + 
                              sph.std("velz").in_units('km/s')**2)
                mass = sph.sum("mass").in_units('Msun')
        
                print(ds.basename, r)
        
                timescale =  2. * r / mean_magvel
        
                #print(outflow_mean, inflow_mean, r, timescale, ds.basename)
                conv_time_out[j, :] = [sim_time, timescale.value, timestep]
                magvel_out[j, :]    = [sim_time, mean_magvel.value,  timestep]
                vrms_out[j, :]    = [sim_time, vrms.value,  timestep]
                conv_zone_out[j, :] = [sim_time, r.value,  timestep]
                conv_mass_out[j, :] = [sim_time, mass.value,  timestep]
        

#sort time array by sim time
sorted_conv_time = conv_time_out[ np.argsort(conv_time_out[:, 0]) ]
sorted_conv_time = sorted_conv_time[ np.where(sorted_conv_time[:, 0] > 0) ]

#sort vel array by sim time
sorted_magvel = magvel_out[ np.argsort(magvel_out[:, 0]) ]
sorted_magvel = sorted_magvel[ np.where(sorted_magvel[:, 0] > 0) ]

#sort vel array by sim time
sorted_vrms = vrms_out[ np.argsort(vrms_out[:, 0]) ]
sorted_vrms = sorted_vrms[ np.where(sorted_vrms[:, 0] > 0) ]

#sort radius array by sim time
sorted_conv_zone = conv_zone_out[ np.argsort(conv_zone_out[:, 0]) ]
sorted_conv_zone = sorted_conv_zone[ np.where(sorted_conv_zone[:, 0] > 0) ]

#sort mass array by sim time
sorted_conv_mass = conv_mass_out[ np.argsort(conv_mass_out[:, 0]) ]
sorted_conv_mass = sorted_conv_mass[ np.where(sorted_conv_mass[:, 0] > 0) ]


#plot it all out
# first time
plt.figure(1)
plt.plot(sorted_conv_time[:, 0], sorted_conv_time[:, 1])

plt.ylabel("Convection Timescale (s)")
plt.xlabel("simulation time (s)")

plt.savefig(f"{outdir}/conv_time_over_time.png")
np.save(f"{outdir}/conv_time_over_time.npy", sorted_conv_time)

#and second plot. magvel
plt.figure(2)
plt.plot(sorted_magvel[:, 0], sorted_magvel[:, 1])

plt.ylabel("avg magnitude velocity (km/s)")
plt.xlabel("simulation time (s)")

plt.savefig(f"{outdir}/magvel_over_time.png")
np.save(f"{outdir}/magvel_over_time.npy", sorted_magvel)

#and third plot. vrms
plt.figure(3)
plt.plot(sorted_vrms[:, 0], sorted_vrms[:, 1])

plt.ylabel("Vrms (km/s)")
plt.xlabel("simulation time (s)")

plt.savefig(f"{outdir}/vrms_over_time.png")
np.save(f"{outdir}/vrms_over_time.npy", sorted_vrms)

# and fourth plot. conv size
plt.figure(4)
plt.plot(sorted_conv_zone[:, 0], sorted_conv_zone[:, 1])

plt.ylabel("Convection Zone Radius (km)")
plt.xlabel("simulation time (s)")

plt.savefig(f"{outdir}/conv_zone_over_time.png")
np.save(f"{outdir}/conv_zone_over_time.npy", sorted_conv_zone)

# and fifth conv mass
plt.figure(5)
plt.plot(sorted_conv_mass[:, 0], sorted_conv_mass[:, 1])

plt.ylabel("Mass of Convection Zone $(M_{\\odot})$")
plt.xlabel("simulation time (s)")

plt.savefig(f"{outdir}/conv_mass_over_time.png")
np.save(f"{outdir}/conv_mass_over_time.npy", sorted_conv_mass)
