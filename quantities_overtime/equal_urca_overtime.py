import yt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import listdir


# okay grab all the files I can
yt.set_log_level(40)

all_files = listdir("plotfiles")
output = np.zeros((len(all_files), 2), dtype=np.float64) -1.

for j, file in enumerate(all_files):
    if "Conv" in file:
        continue
    else:
        curr_t = yt.load(f"plotfiles/{file}", hint='amrex').current_time
        
    try:
        curr_df = pd.read_csv(f"profiles/{file}.csv", index_col=0)
    except FileNotFoundError:
        continue
    i_min = np.argmin(np.abs(curr_df['X(na23)'] - curr_df['X(ne23)']))
    output[j, :] = [curr_t, curr_df['radius'].iloc[i_min]]

#clean up numpy array
sorted_output = output[ np.argsort(output[:, 0]) ]
sorted_output = sorted_output[ np.where(sorted_output[:, 0] > 0) ]

#plot it all out
plt.plot(sorted_output[:, 0], sorted_output[:, 1]/1.e5)
np.save("equal_over_time.npy", sorted_output) 

plt.title("Average radius where X(na23) = X(ne23)")
plt.ylabel("Radius (km)")
plt.xlabel("time (s)")

plt.savefig("equal_over_time.png")
