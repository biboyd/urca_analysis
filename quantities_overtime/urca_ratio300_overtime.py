from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt

from sys import argv


yt.set_log_level(40)
time_arr = []
ratio_arr = []
timestep_arr = []

if len(argv) == 1:
    plot_dir = "plotfiles"
else:
    plot_dir = argv[1]

try:
    old_ratio_file = "ratio300_over_time.npy"
    old_ratio_arr = np.load(old_ratio_file)
    calc_all = False
    
except FileNotFoundError:
    calc_all = True
    print("Could not find old files, will be calculating for each file")

for f in listdir(plot_dir):
    if f[:3] == "plt" and len(f) == 10:
        try:
            ds = yt.load(f"{plot_dir}/{f}", hint='amrex')

            time = ds.current_time.value
            timestep = float(ds.basename.removeprefix("plt"))
        except FileNotFoundError:
            continue

    # check if already calculated
    if not calc_all and timestep in old_ratio_arr[:, 2]: 
        idx = np.argmin(np.abs(timestep - old_ratio_arr[:, 2]))
        time_arr.append(old_ratio_arr[idx, 0])
        ratio_arr.append(old_ratio_arr[idx, 1])
        timestep_arr.append(old_ratio_arr[idx, 2])

    else:
        sph = ds.sphere(ds.domain_center, (300, 'km'))
        ratio = sph.mean("X(na23)")/sph.mean("X(ne23)")
        
        time_arr.append(time)
        ratio_arr.append(ratio)
        timestep_arr.append(timestep)

output =  np.array([time_arr, ratio_arr, timestep_arr]).swapaxes(0, 1)
sorted_output = output[ np.argsort(output[:, 0]) ]
np.save("ratio300_over_time.npy", sorted_output)
fig, ax = plt.subplots(1, 1)

ax.plot(time_arr, ratio_arr, 'o')

#ax.set_ylim(0.07, 0.09)
#ax.set_xlim(4000)
ax.set_xlabel('time (s)')
ax.set_title("average $X({}^{23} Na)$/$X({}^{23} Ne)$ \nfrom 0 to 300 km in radius")

ax.grid()
fig.tight_layout()


fig.savefig(f"ratio300_over_time.png")
