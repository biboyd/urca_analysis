import yt
import unyt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

@yt.derived_field(name=("gas", "conv_work"), units='erg/s', sampling_type='local')
def _energy_rate(field, data):
    return data[('boxlib', 'eps_conv')] * data[('gas', 'cell_volume')] * unyt.erg/unyt.cm**3/unyt.s

ds = yt.load(argv[1])

# calc and save energy rate profile
N_bins = 100
fields = [('gas', 'conv_work')]
tot_nuc_prof = yt.create_profile(ds.all_data(), 'radius', fields, weight_field=None, logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})
df = tot_nuc_prof.to_dataframe()
df.to_csv(f"convwork_profiles/{ds.basename}_convwork.csv")

# avg this out
avg_prof = yt.create_profile(ds.all_data(), 'radius', ds.field_list, weight_field='volume', logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})
df = avg_prof.to_dataframe()
df.to_csv(f"convwork_profiles/{ds.basename}_profiles.csv")
