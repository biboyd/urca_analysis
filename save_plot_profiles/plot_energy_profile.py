import yt
import unyt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

@yt.derived_field(name=("gas", "mass"), sampling_type='local')
def _mass(field, data):
    return data["boxlib", "rho"]*unyt.g/unyt.cm**3 * data['boxlib', 'volume']

@yt.derived_field(name=("gas", "energy_rate"), sampling_type='local')
def _energy_rate(field, data):
    return data[('boxlib', 'Hnuc')] * data[('gas', 'mass')] * unyt.erg/unyt.g/unyt.s

ds = yt.load(argv[1])

# calc and save energy rate profile
N_bins = 100
tot_nuc_prof = yt.create_profile(ds.all_data(), 'radius', ('gas', "energy_rate"), weight_field=None, logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})
df = tot_nuc_prof.to_dataframe()
df.to_csv(f"{ds.basename}_tot_energy_loss.csv")

# plot
fig, ax = plt.subplots(1, 1)
binsize = (df['radius'].iloc[1] - df['radius'].iloc[0])/1e5 #in km
ax.plot(df['radius']/1e5, df['energy_rate']/binsize, c='tab:red',label="Nuclear Reactions")

plt.grid()
ax.set_xlim(0, 600)

ax.set_title("Total Energy Generation Rate per Radial Bin")
ax.set_ylabel("$\\dot{E}(r)$ ($10^{41} \\; \\mathrm{erg}~\\mathrm{s}^{-1}~\\mathrm{km}^{-1}$)")
ax.set_xlabel("Radius (km)")

# annotate plot with key reactions
#c12 burn
text_fs = 13
text_bbox = dict(fc="white", ec="white", alpha=0.5, pad=0.)
ax.annotate("   $\\mathrm{{}^{12} C}$ \n burning",(50, 1.5),
            fontsize=text_fs, color='k', 
            bbox=text_bbox, zorder=2)
#na23 capture
ax.annotate("      $\\mathrm{{}^{23} Na}$\n $e^{-}-\\mathrm{capture}$",(189, 0.25),
            fontsize=text_fs, color='k', 
            bbox=text_bbox, zorder=2)
#ne23 decay
ax.annotate("    $\\mathrm{{}^{23} Ne}$\n $\\beta - \\mathrm{decay}$",(450, 1.25),
            fontsize=text_fs, color='k', 
            bbox=text_bbox, zorder=2)

fig.savefig(f"{ds.basename}_tot_energy_loss.png", dpi=300)
