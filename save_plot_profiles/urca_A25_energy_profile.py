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
Mg25_nuc = pyna.Nucleus('Mg25')
delta = (Mg25_nuc.Z/Mg25_nuc.A - 0.5)*(m_p + m_e - m_n) *unyt.g 
Mg25_bind = Mg25_nuc.nucbind*unyt.MeV - delta*unyt.c**2

Na25_nuc = pyna.Nucleus('Na25')
delta = (Na25_nuc.Z/Na25_nuc.A - 0.5)*(m_p + m_e - m_n) *unyt.g 
Na25_bind = Na25_nuc.nucbind*unyt.MeV - delta*unyt.c**2


@yt.derived_field(name=("gas", "mass"), units='g', sampling_type='local')
def _mass(field, data):
    return data["boxlib", "rho"]* data['boxlib', 'volume'] * unyt.g/unyt.cm**3

@yt.derived_field(name=("gas", "urca_energy_rate"), units='erg/s', sampling_type='local')
def _energy_rate(field, data):
    return (data[('boxlib', 'A25_electron_capture_rate')]*(Na25_bind-Mg25_bind) + data[('boxlib', 'A25_beta_decay_rate')]*(Mg25_bind-Na25_bind)) * data[('gas', 'mass')] * unyt.avogadros_number_mks.value / unyt.g/ unyt.s 

A25_nu_list = [('boxlib', 'A25_beta_decay_nu_loss'), ('boxlib', 'A25_electron_capture_nu_loss')]

@yt.derived_field(name=("gas", "nuloss_rate"), units='erg/s', sampling_type='local')
def _tot_nu_loss(field, data):
    tot=0. 
    for nuloss in A25_nu_list:
        tot -= data[nuloss]
    tot+=data['thermal_nu_loss']
    tot *= unyt.erg/unyt.g/unyt.s
    return  tot*data[('gas', 'mass')]


# load data 
ds = yt.load(argv[1]) #make sure it's nuloss file

# calc and save energy rate profile
N_bins = 100
fields = [('gas', 'urca_energy_rate'), ('gas', 'nuloss_rate')]
tot_nuc_prof = yt.create_profile(ds.all_data(), 'radius', fields, weight_field=None, logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})
df = tot_nuc_prof.to_dataframe()
df.to_csv(f"eos_profiles/{ds.basename}_urca_A25_energy.csv")

