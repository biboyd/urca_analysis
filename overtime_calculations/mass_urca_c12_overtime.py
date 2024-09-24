import yt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from os import listdir
import pandas as pd

def _urca_mass(field, data):
    return (data[('boxlib', 'X(na23)')] + data[('boxlib', 'X(ne23)')]) * data[('gas', 'mass')]

def _c12_mass(field, data):
    return (data[('boxlib', 'X(c12)')]) * data[('gas', 'mass')]

def _urca_mass_rate(field, data):
    return data[('gas', 'mass')] * (data[('boxlib', 'omegadot(ne23)')] + data[('boxlib', 'omegadot(na23)')])

def _c12_mass_rate(field, data):
    return data[('gas', 'mass')] * (data[('boxlib', 'omegadot(c12)')])

# okay grab all the files I can
yt.set_log_level(40)

if len(argv) == 1:
    plots_dir = "small_plotfiles"
else:
    plots_dir = argv[1]

all_files = listdir(plots_dir)

conv_c12mass_out = np.zeros((len(all_files), 3), dtype=np.float64) -1. 
conv_c12rate_out = np.zeros((len(all_files), 3), dtype=np.float64) -1. 
conv_urca_mass_out = np.zeros((len(all_files), 3), dtype=np.float64) -1. 
conv_urca_rate_out = np.zeros((len(all_files), 3), dtype=np.float64) -1. 
tot_c12mass_out = np.zeros((len(all_files), 3), dtype=np.float64) -1. 
tot_c12rate_out = np.zeros((len(all_files), 3), dtype=np.float64) -1. 
tot_urca_mass_out = np.zeros((len(all_files), 3), dtype=np.float64) -1. 
tot_urca_rate_out = np.zeros((len(all_files), 3), dtype=np.float64) -1. 
all_arr_out = [conv_c12mass_out, conv_c12rate_out, conv_urca_mass_out, conv_urca_rate_out,
               tot_c12mass_out, tot_c12rate_out, tot_urca_mass_out, tot_urca_rate_out] 

#previous data
try:
    old_conv_c12mass = np.load("conv_c12mass_over_time.npy")
    old_conv_c12rate = np.load("conv_c12rate_over_time.npy")
    old_conv_urca_mass = np.load("conv_urca_mass_over_time.npy")
    old_conv_urca_rate = np.load("conv_urca_rate_over_time.npy")
    
    old_tot_c12mass = np.load("tot_c12mass_over_time.npy")
    old_tot_c12rate = np.load("tot_c12rate_over_time.npy")
    old_tot_urca_mass = np.load("tot_urca_mass_over_time.npy")
    old_tot_urca_rate = np.load("tot_urca_rate_over_time.npy")
    old_all_arr = [old_conv_c12mass, old_conv_c12rate, old_conv_urca_mass, old_conv_urca_rate,
                   old_tot_c12mass, old_tot_c12rate, old_tot_urca_mass, old_tot_urca_rate] 
    calc_all = False
    
except FileNotFoundError:
    calc_all = True
    print("Could not find old files, will be calculating for each file")

for j, file in enumerate(all_files):
    if file[:3] == "plt" and len(file) == 10:
        
        ds = yt.load(f"{plots_dir}/{file}", hint='maestro')
        sim_time = ds.current_time
        timestep = float(ds.basename.removeprefix("plt"))
        
        #ignore early times
        if sim_time < ds.quan(20, 's'):
            continue
        else:
            ds.add_field(('gas', 'urca_mass'), _urca_mass, 'local', units='Msun', force_override=True)
            ds.add_field(('gas', 'c12_mass'), _c12_mass, 'local', units='Msun', force_override=True)
            ds.add_field(('gas', 'urca_mass_rate'), _urca_mass_rate, 'local', units='Msun/hr', force_override=True)
            ds.add_field(('gas', 'c12_mass_rate'), _c12_mass_rate, 'local', units='Msun/hr', force_override=True)
            
            if not calc_all:
                # check if already calculated
                t_inall_arr = []
                for arr in old_all_arr:
                    t_inall_arr.append(timestep in arr)
                
                if np.all(t_inall_arr):
                    for old_arr, arr_out in zip(old_all_arr, all_arr_out):
                        idx = np.argmin(np.abs(timestep - old_arr[:, 2]))
                        arr_out[j, :] = old_arr[idx, :]

            else:
                # calc convect zone
                try:
                    conv_rad_arr = np.load("conv_zone_over_time.npy")
                    idx = np.argmin(np.abs(timestep - conv_rad_arr[:, 2]))
                except FileNotFoundError:
                    print(f"conv zone file found for {file}")
                    continue

                r = ds.quan(conv_rad_arr[idx, 1], 'km')

                    
                try:
                    conv_sph = ds.sphere(ds.domain_center, r)
                except yt.utilities.exceptions.YTSphereTooSmall:
                    continue
                
                #full wd sphere
                tot_sph = ds.sphere(ds.domain_center, (1.6e3, 'km'))
                
                conv_c12_mass  = conv_sph.sum(('gas', 'c12_mass'))
                conv_c12_rate  = conv_sph.sum(('gas', 'c12_mass_rate'))
        
                conv_urca_mass = conv_sph.sum(('gas', 'urca_mass'))
                conv_urca_rate = conv_sph.sum(('gas', 'urca_mass_rate'))
                
                tot_c12_mass  = tot_sph.sum(('gas', 'c12_mass'))
                tot_c12_rate  = tot_sph.sum(('gas', 'c12_mass_rate'))
                
                tot_urca_mass = tot_sph.sum(('gas', 'urca_mass'))
                tot_urca_rate = tot_sph.sum(('gas', 'urca_mass_rate'))
                
                print(ds.basename)
        
                #print(outflow_mean, inflow_mean, r, timescale, ds.basename)
                conv_c12mass_out[j, :]   = [sim_time, conv_c12_mass.value,  timestep]
                conv_c12rate_out[j, :]   = [sim_time, conv_c12_rate.value,  timestep]
                conv_urca_mass_out[j, :] = [sim_time, conv_urca_mass.value, timestep]
                conv_urca_rate_out[j, :] = [sim_time, conv_urca_rate.value, timestep]
                tot_c12mass_out[j, :]    = [sim_time, tot_c12_mass.value,   timestep]
                tot_c12rate_out[j, :]    = [sim_time, tot_c12_rate.value,   timestep]
                tot_urca_mass_out[j, :]  = [sim_time, tot_urca_mass.value,  timestep]
                tot_urca_rate_out[j, :]  = [sim_time, tot_urca_rate.value,  timestep]
                
        

# sort the arrays
all_names_out = ["conv_c12mass", "conv_c12rate", "conv_urca_mass", "conv_urca_rate",
                 "tot_c12mass", "tot_c12rate", "tot_urca_mass", "tot_urca_rate"] 

for arr_out, name in zip(all_arr_out, all_names_out):
    #sort and remove empty
    arr_out = arr_out[ np.argsort(arr_out[:, 0]) ]
    arr_out = arr_out[ np.where(arr_out[:, 0] > 0) ]
    
    if 'mass' in name:
        units = 'Msun'
    elif 'rate' in name:
        units = 'Msun/hr'
    else:
        raise RuntimeError(f"Need to add units to name, {name}, that was added")
        
    fig, ax  = plt.subplots(1,1)
    ax.plot(arr_out[:, 0], arr_out[:, 1])
    ax.set_ylabel(f"{name} ({units})")
    ax.set_xlabel("simulation time (s)")
    fig.savefig(f"{name}_over_time.png", bbox_inches='tight')
    np.save(f"{name}_over_time.npy", arr_out)