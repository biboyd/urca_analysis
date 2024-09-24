from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt

from sys import argv


yt.set_log_level(40)
time_arr = []
na23_arr = []
ne23_arr = []

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
    na23_curr = np.mean(df["X(na23)"].iloc[:i])
    ne23_curr = np.mean(df["X(ne23)"].iloc[:i])

    time_arr.append(time)
    na23_arr.append(na23_curr)
    ne23_arr.append(ne23_curr)

fig, ax = plt.subplots(1, 1)

ax.plot(time_arr, na23_arr, 'o', label= "$X({}^{23} Na)$")
ax.plot(time_arr, ne23_arr, 'o', label="$X({}^{23} Na)$")

ax.set_yscale('log')
ax.set_xlabel('time (s)')
ax.set_title("average Mass Fraction \nfrom 0 to 300 km in radius")
ax.legend()
ax.grid()
fig.tight_layout()


fig.savefig(f"na_ne23_over_time.png")
