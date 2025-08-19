import yt
import unyt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

@yt.derived_field(name=("gas", "mass"), units='g', sampling_type='local')
def _mass(field, data):
    return data["boxlib", "density"]* data['boxlib', 'volume']

@yt.derived_field(name=("gas", "therm_energy_rate"), units='erg/s', sampling_type='local')
def _energy_rate(field, data):
    return data[('boxlib', 'therm_term_h')] * data[('gas', 'mass')] * unyt.erg/unyt.g/unyt.s

@yt.derived_field(name=("gas", "pres_energy_rate"), units='erg/s', sampling_type='local')
def _energy_rate(field, data):
    return data[('boxlib', 'pres_term')] * data[('gas', 'mass')] * unyt.erg/unyt.g/unyt.s


def plot_therm_energy(df):
    # plot therm energy
    fig, ax = plt.subplots(1, 1)
    binsize = (df['radius'].iloc[1] - df['radius'].iloc[0])/1e5 #in km
    ax.plot(df['radius']/1e5, df['therm_energy_rate']/binsize, c='tab:red',label="Thermal component")
    ax.plot(df['radius']/1e5, df['pres_energy_rate']/binsize, c='tab:red',label="pressure component")
    
    plt.grid()
    
    ax.set_title("Total Energy Generation Rate per Radial Bin")
    ax.set_ylabel("$\\dot{E}(r)$ ($10^{41} \\; \\mathrm{erg}~\\mathrm{s}^{-1}~\\mathrm{km}^{-1}$)")
    ax.set_xlabel("Radius (km)")
    
    # annotate plot with key reactions
    #c12 burn
    """text_fs = 13
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
                bbox=text_bbox, zorder=2)"""
    
    fig.savefig(f"{ds.basename}_therm_and_pres_energy.png", dpi=300)


ds = yt.load(argv[1])

# calc and save energy rate profile
N_bins = 100
fields = [('gas', 'pres_energy_rate'), ('gas', 'therm_energy_rate')]
tot_nuc_prof = yt.create_profile(ds.all_data(), 'radius', fields, weight_field=None, logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})
df = tot_nuc_prof.to_dataframe()
df.to_csv(f"eos_profiles/{ds.basename}_energy.csv")

# quick plots
plot_therm_energy(df)
