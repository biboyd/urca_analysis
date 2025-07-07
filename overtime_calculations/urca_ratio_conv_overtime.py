from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt

from sys import argv


yt.set_log_level(40)

# grab all the files
if len(argv) == 1:
    plot_dir = "plotfiles"
else:
    plot_dir = argv[1]
all_files = listdir(plot_dir)

urca_ratio_out = np.zeros((len(all_files), 3), dtype=np.float64) -1.

try:
    old_ratio_arr = np.load("urca_ratio_over_time.npy")
    calc_all = False
    
except FileNotFoundError:
    calc_all = True
    print("Could not find old files, will be calculating for each file")

# grab convect zone array
try:
    conv_zone_arr = np.load("conv_zone_over_time.npy")
except FileNotFoundError:
    print("Could not find the conv zone file. Make sure to run convect_overtim.py first.")
    raise FileNotFoundError("conv_zone_over_time.npy")

for j, f in enumerate(all_files):
    if f[:3] == "plt" and len(f) == 10:
        try:
            ds = yt.load(f"{plot_dir}/{f}", hint='amrex')

            time = ds.current_time.value
            timestep = float(ds.basename.removeprefix("plt"))
        except FileNotFoundError:
            continue

        #ignore initial time stuff
        if time > 20.:
            # check if entry is in Rconv array
            assert timestep in conv_zone_arr[:, 2], f"{timestep} not in conv zone array"

            # check if already calculated
            if not calc_all and timestep in old_ratio_arr[:, 2]: 
                idx = np.argmin(np.abs(timestep - old_ratio_arr[:, 2]))
                urca_ratio_out[j, :] = old_ratio_arr[idx, :]

            else:
                idx = np.argmin(np.abs(timestep - conv_zone_arr[:, 2]))
                radius = conv_zone_arr[idx, 1]
                
                sph = ds.sphere(ds.domain_center, (radius, 'km'))
                curr_ratio = sph.mean("rhoX(Ne23)")/sph.mean("rhoX(Na23)")
                
                urca_ratio_out[j, :] = [time, curr_ratio, timestep]

sorted_urca_ratio = urca_ratio_out[ np.argsort(urca_ratio_out[:, 0]) ]
sorted_urca_ratio = sorted_urca_ratio[ np.where(sorted_urca_ratio[:, 0] > 0) ]

np.save("urca_ratio_over_time.npy", sorted_urca_ratio)
fig, ax = plt.subplots(1, 1)

ax.plot(sorted_urca_ratio[:, 0],sorted_urca_ratio[:, 1], 'o')

#ax.set_ylim(0.07, 0.09)
#ax.set_xlim(4000)
ax.set_xlabel('time (s)')
ax.set_title("average dens $X({}^{23} Ne)$/$X({}^{l3} Na)$ \nin Conv zone")

ax.grid()
fig.tight_layout()


fig.savefig(f"urca_ratio_over_time.png")
