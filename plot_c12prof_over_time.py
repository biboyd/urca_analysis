from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt

from matplotlib import colormaps
from sys import argv


cmap = colormaps['viridis']
max_tstep=163536
yt.set_log_level(40)
fig, ax = plt.subplots(1, 1)

for f in listdir("profiles"):
    if "all_fields" not in f:
        pass
    else:

        nocsv = f.replace("_all_fields.csv", '')
        noprofile = nocsv.replace('profiles/', '')
        basename = noprofile.replace('plt', '')

        curr_color = cmap(float(basename)/max_tstep)
        df = pd.read_csv(f"profiles/{f}",index_col=0).dropna()
        if 'X(c12)' in df.columns:
            ax.plot(df['radius']/1e5, df['X(c12)'], color=curr_color)
        else:
            pass

ax.set_ylabel("$X({}^{12} C)$")
ax.set_xlabel('Radius (km)')
ax.set_title("Carbon Radially Averaged Profile")
ax.set_xlim(0., 600.)
ax.grid()
fig.tight_layout()


fig.savefig(f"c12_prof_overtime.png")
