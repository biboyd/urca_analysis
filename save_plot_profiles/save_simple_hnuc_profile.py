import yt
import unyt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

ds = yt.load(argv[1])

# calc and save energy rate profile
N_bins = 100
fields = ['Hnuc']
tot_nuc_prof = yt.create_profile(ds.all_data(), 'radius', fields, logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})
df = tot_nuc_prof.to_dataframe()
df.to_csv(f"{ds.basename}_simple_hnuc_profile.csv")

