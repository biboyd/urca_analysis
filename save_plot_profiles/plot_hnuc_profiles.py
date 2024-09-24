from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt

from sys import argv


outdir='plots_profiles'
f = argv[1]
if len(argv) == 3 :
    pfname = argv[2]
    
df = pd.read_csv(f"{f}",index_col=0).dropna()

rad = df.radius.to_numpy(np.float64)/1e5 #km
hnuc = df['Hnuc'].to_numpy(np.float64)

fig, ax = plt.subplots(1, 1)

ax.plot(rad, hnuc)
    
#ax.set_yscale('log')
bot, top = ax.set_ylim(-5e11, 2.1e12)
ax.set_xlim(0, 600)


ax.set_ylabel('Avg. Specific Energy Generation (erg/g/s)')
ax.set_xlabel("Radius (km)")

ax.grid()
fig.tight_layout()

fname = f.replace(".csv", '').replace('profiles/', '')

try:
    if len(argv) == 3:
        ds = yt.load(f"{pfname}", hint='amrex')
    else:
        ds = yt.load(f"plotfiles/{fname}", hint='amrex')
    
    time_str = f"t = {ds.current_time:0.2f} s"
except FileNotFoundError:
    time_str = "t = - s"

ax.set_title(time_str)
fig.savefig(f"{outdir}/{fname}_hnuc_profile.png", dpi=300)
plt.close(fig)
