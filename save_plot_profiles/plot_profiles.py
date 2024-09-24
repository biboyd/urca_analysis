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
Xna23 = df['X(na23)'].to_numpy(np.float64)
Xne23 = df['X(ne23)'].to_numpy(np.float64)

fig, ax = plt.subplots(1, 1)

ax.plot(rad, Xna23, label = '${}^{23}\mathrm{Na}$')
ax.plot(rad, Xne23, label = '${}^{23}\mathrm{Ne}$')
    
equal = rad[np.argmin(abs(Xna23 - Xne23))]
#ax.set_yscale('log')
bot, top = ax.set_ylim(0, 6e-4)
ax.set_xlim(0, 600)

ax.vlines(equal, bot, top, color='k', linestyle='--', alpha=0.7, label='Equal Mass Fractions', zorder=-1)

ax.set_ylabel("Mass Fraction")
ax.set_xlabel("Radius (km)")

ax.legend()
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
fig.savefig(f"{outdir}/{fname}_urca_profile.png", dpi=300)
plt.close(fig)

# plot the againt rho
fig, ax = plt.subplots(1, 1)

rho = df.rho.to_numpy(np.float64)
ax.plot(rho, Xna23, label = '${}^{23}\mathrm{Na}$')
ax.plot(rho, Xne23, label = '${}^{23}\mathrm{Ne}$')
    
equal = rho[np.argmin(abs(Xna23 - Xne23))]
#ax.set_yscale('log')
bot, top = ax.set_ylim(0, 1e-3)

ax.set_ylabel("Mass Fraction")
ax.set_xlabel("Density ($g/cm^3$)")
    
ax.vlines(equal, bot, top, color='k', linestyle='--', alpha=0.7, label='Equal Mass Fractions', zorder=-1)

ax.invert_xaxis()
ax.legend()
fig.tight_layout()

ax.set_title(time_str)
fig.savefig(f"{outdir}/{fname}_urca_density_profile.png", dpi=300)
plt.close(fig)

