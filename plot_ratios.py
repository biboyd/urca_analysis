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

try:
    old_ratio_file = "ratio_over_time.npy"
    old_ratio_arr = np.load(old_ratio_file)
    calc_all = False
    
except FileNotFoundError:
    calc_all = True
    print("Could not find old files, will be calculating for each file")

for f in listdir("plotfiles"):
    if f[:3] == "plt" and not f[-3:] == "_nu" and not f[-5:] == "_conv":
        fname = f.replace(".csv", '').replace('profiles/', '')
        try:
            ds = yt.load(f"plotfiles/{f}", hint='amrex')

            time = ds.current_time.value
            timestep = float(ds.basename.removeprefix("plt"))
        except FileNotFoundError:
            continue
            # check if already calculated

    if not calc_all and timestep in old_ratio_arr[:, 2]: 
        idx = np.argmin(np.abs(time - old_ratio_arr[:, 0]))
        time_arr.append(old_ratio_arr[idx, 0])
        ratio_arr.append(old_ratio_arr[idx, 1])
        timestep_arr.append(old_ratio_arr[idx, 2])
        continue

    """
    df = pd.read_csv(f"profiles/{f}",index_col=0).dropna()

    rad = df.radius.to_numpy(np.float64)/1e5 #km
    i = np.argmin(np.abs(rad - 300))
    ratio = np.mean(df["X(na23)"].iloc[:i])/np.mean(df["X(ne23)"].iloc[:i])
    """
    sph = ds.sphere(ds.domain_center, (300, 'km'))
    ratio = sph.mean("X(na23)")/sph.mean("X(ne23)")
    

    time_arr.append(time)
    ratio_arr.append(ratio)
    timestep_arr.append(timestep)

output =  np.array([time_arr, ratio_arr, timestep_arr]).swapaxes(0, 1)
sorted_output = output[ np.argsort(output[:, 0]) ]
np.save("ratio_over_time.npy", sorted_output)
fig, ax = plt.subplots(1, 1)

ax.plot(time_arr, ratio_arr, 'o')

#ax.set_ylim(0.07, 0.09)
#ax.set_xlim(4000)
ax.set_xlabel('time (s)')
ax.set_title("average $X({}^{23} Na)$/$X({}^{23} Ne)$ \nfrom 0 to 300 km in radius")

ax.grid()
fig.tight_layout()


fig.savefig(f"ratio_over_time.png")
