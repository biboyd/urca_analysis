from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt

from sys import argv


yt.set_log_level(40)
time_arr = []
Ye_arr = []
for f in listdir("profiles"):
    fname = f.replace(".csv", '').replace('profiles/', '')
    try:
        ds = yt.load(f"plotfiles/{fname}", hint='amrex')
    
        time = ds.current_time
    except FileNotFoundError:
        continue

    df = pd.read_csv(f"profiles/{f}",index_col=0).dropna()

    rad = df.radius.to_numpy(np.float64)/1e5 #km
    i = np.argmin(np.abs(rad - 300))
    Ye = ( np.mean(df["Ye"].iloc[:i]) - 0.5 )/1e-5

    time_arr.append(time)
    Ye_arr.append(Ye)
fig, ax = plt.subplots(1, 1)

ax.plot(time_arr, Ye_arr, 'o')

ax.set_xlabel('time (s)')
ax.set_ylabel("average Ye from 0 to 300 km in radius")

ax.grid()
fig.tight_layout()


fig.savefig(f"avg_Ye_over_time.png")
