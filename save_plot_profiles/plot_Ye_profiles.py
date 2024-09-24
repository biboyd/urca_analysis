from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yt

from sys import argv

indir = 'profiles'
outdir='plots_profiles'

for f in listdir(indir):
    # check if have plotfile
    fname = f.replace(".csv", '').replace('profiles/', '')

    try:
        ds = yt.load(f"plotfiles/{fname}", hint='amrex')
    
    except FileNotFoundError:
        continue

    time_str = f"t = {ds.current_time:0.2f} s"
    df = pd.read_csv(f"{indir}/{f}",index_col=0).dropna()

    rad = df.radius.to_numpy(np.float64)/1e5 #km
    Ye= df['Ye'].to_numpy(np.float64)

    fig, ax = plt.subplots(1, 1)

    ax.plot(rad, (Ye - 0.5)/1e-5)
    
    #ax.set_yscale('log')
    bot, top = ax.set_ylim(-5.5, -1)
    ax.set_xlim(0, 600)


    ax.set_ylabel('Radially Avg. (Ye - 0.5)/1e-5 ')
    ax.set_xlabel("Radius (km)")

    ax.grid()
    fig.tight_layout()

    ax.set_title(time_str)
    fig.savefig(f"{outdir}/{fname}_Ye_profile.png", dpi=300)
    plt.close(fig)
