import yt
import unyt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import pynucastro as pyna

m_p = 1.672621777e-24 #g
m_n = 1.674927351e-24 #g
m_e = 9.10938291e-28 #g

path_to_pfile = "/pscratch/sd/b/biboyd/gpuruns/bstate_on_urca_0.9Mconv_large"
ds = yt.load(f"{path_to_pfile}/plotfiles/plt0008400/")

# calc bind diffs
Na23_nuc = pyna.Nucleus('Na23')
delta = (Na23_nuc.Z/Na23_nuc.A - 0.5)*(m_p + m_e - m_n) *unyt.g 
Na23_bind = Na23_nuc.nucbind*unyt.MeV - delta*unyt.c**2

Ne23_nuc = pyna.Nucleus('Ne23')
delta = (Ne23_nuc.Z/Ne23_nuc.A - 0.5)*(m_p + m_e - m_n) *unyt.g 
Ne23_bind = Ne23_nuc.nucbind*unyt.MeV - delta*unyt.c**2


@yt.derived_field(name=("gas", "mass"), units='g', sampling_type='local')
def _mass(field, data):
    return data["boxlib", "rho"]* data['boxlib', 'volume'] * unyt.g/unyt.cm**3

@yt.derived_field(name=("gas", "urca_energy_rate"), units='erg/s', sampling_type='local')
def _energy_rate(field, data):
    return (data[('boxlib', 'A23_electron_capture_rate')]*(Ne23_bind-Na23_bind) + data[('boxlib', 'A23_beta_decay_rate')]*(Na23_bind-Ne23_bind)) * data[('gas', 'mass')] * unyt.avogadros_number_mks.value / unyt.g/ unyt.s 


def plot_urca_energy(df):
    # plot urca energy. no neutrino losses
    fig, ax = plt.subplots(1, 1)
    binsize = (df['radius'].iloc[1] - df['radius'].iloc[0])/1e5 #in km
    ax.plot(df['radius']/1e5, df['urca_energy_rate']/binsize, c='tab:red',label="Urca energy (no nu's)")
    
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
    
    fig.savefig(f"{ds.basename}_tot_nuc_energy.png", dpi=300)


ds = yt.load(argv[1])

# calc and save energy rate profile
N_bins = 100
fields = [('gas', 'urca_energy_rate')]
tot_nuc_prof = yt.create_profile(ds.all_data(), 'radius', fields, weight_field=None, logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})
df = tot_nuc_prof.to_dataframe()
df.to_csv(f"eos_profiles/{ds.basename}_urca_energy.csv")

# quick plots
plot_urca_energy(df)
