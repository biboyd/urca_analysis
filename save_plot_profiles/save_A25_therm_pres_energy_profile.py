import yt
import unyt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

@yt.derived_field(name=("gas", "mass"), units='g', sampling_type='local')
def _mass(field, data):
    return data["boxlib", "density"]* data['boxlib', 'volume']

@yt.derived_field(name=("gas", "therm(Na25)_energy_rate"), units='erg/s', sampling_type='local')
def _energy_rate(field, data):
    return data[('boxlib', 'therm(Na25)')] * data[('gas', 'mass')] * unyt.erg/unyt.g/unyt.s

@yt.derived_field(name=("gas", "pres(Na25)_energy_rate"), units='erg/s', sampling_type='local')
def _energy_rate(field, data):
    return data[('boxlib', 'pres(Na25)')] * data[('gas', 'mass')] * unyt.erg/unyt.g/unyt.s

@yt.derived_field(name=("gas", "therm(Mg25)_energy_rate"), units='erg/s', sampling_type='local')
def _energy_rate(field, data):
    return data[('boxlib', 'therm(Mg25)')] * data[('gas', 'mass')] * unyt.erg/unyt.g/unyt.s

@yt.derived_field(name=("gas", "pres(Mg25)_energy_rate"), units='erg/s', sampling_type='local')
def _energy_rate(field, data):
    return data[('boxlib', 'pres(Mg25)')] * data[('gas', 'mass')] * unyt.erg/unyt.g/unyt.s


ds = yt.load(argv[1])

# calc and save energy rate profile
N_bins = 100
fields = [('gas', 'pres(Na25)_energy_rate'), ('gas', 'therm(Na25)_energy_rate'),
          ('gas', 'pres(Mg25)_energy_rate'), ('gas', 'therm(Mg25)_energy_rate')]
tot_nuc_prof = yt.create_profile(ds.all_data(), 'radius', fields, weight_field=None, logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})
df = tot_nuc_prof.to_dataframe()
df.to_csv(f"eos_profiles/{ds.basename}_A25_energy.csv")

