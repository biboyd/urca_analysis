import yt
import unyt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

@yt.derived_field(name=("gas", "mass"), units='g', sampling_type='local')
def _mass(field, data):
    return data["boxlib", "rho"]*unyt.g/unyt.cm**3 * data['boxlib', 'volume']

@yt.derived_field(name=("gas", "energy_rate"), units='erg/s', sampling_type='local')
def _energy_rate(field, data):
    return data[('boxlib', 'Hnuc')] * data[('gas', 'mass')] * unyt.erg/unyt.g/unyt.s


@yt.derived_field(name=("gas", "kinetic_energy"), units='erg', sampling_type='local')
def _energy_rate(field, data):
    return 0.5 * (data[('boxlib', 'velx')]**2 + data[('boxlib', 'vely')]**2 + data[('boxlib', 'velz')]**2)   * data[('gas', 'mass')] * unyt.cm**2/unyt.s**2

def plot_nuc_energy(df):
    # plot nuc energy
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
    
    fig.savefig(f"{ds.basename}_tot_nuc_energy.png", dpi=300)


def plot_kin_energy(df):
    # plot
    fig, (axt, axb) = plt.subplots(2, 1)
    binsize = (df['radius'].iloc[1] - df['radius'].iloc[0])/1e5 #in km
    axt.plot(df['radius']/1e5, df['kinetic_energy']/binsize)
    axb.plot(df['radius']/1e5, np.cumsum(df['kinetic_energy']))
    
    plt.grid()
    axt.set_xlim(0, 600)
    axb.set_xlim(0, 600)
    
    axt.set_title("Kinetic Energy ")
    axt.set_xlabel("Radius (km)")
    axt.set_ylabel("KE per Radial Bin (erg/km)")
    
    axb.set_xlabel("Radius (km)")
    axb.set_ylabel("Cumulative KE Inside Radius (erg)")
    
    fig.tight_layout()
    fig.savefig(f"{ds.basename}_kin_energy.png", dpi=300)


ds = yt.load(argv[1])

# calc and save energy rate profile
N_bins = 100
fields = [('gas', 'energy_rate'), ('gas', 'kinetic_energy')]
tot_nuc_prof = yt.create_profile(ds.all_data(), 'radius', fields, weight_field=None, logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})
df = tot_nuc_prof.to_dataframe()
df.to_csv(f"{ds.basename}_nuc_and_kin_energy.csv")

# quick plots
plot_nuc_energy(df)
plot_kin_energy(df)
