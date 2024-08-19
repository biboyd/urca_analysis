import yt
import unyt
import pandas as pd
from sys import argv

def _mass(field, data):
    return data["boxlib", "rho"]*unyt.g/unyt.cm**3 *data['boxlib', 'volume']
def save_profile(fname):
    ds = yt.load(fname, hint='maestro')
    #generate profiles
    #prof_fields=[('boxlib', 'X(ne23)'), ('boxlib', 'X(na23)'), ('boxlib', 'X(c12)'), ('boxlib', 'X(o16)'), ('boxlib', 'tfromp'), ('boxlib', 'vort'), ('boxlib', 'radial_velocity'), ('boxlib', 'rho'), ('boxlib', 'p0pluspi'), ('boxlib', 'Hnuc'), ('boxlib', 'entropy')]
    
    prof_fields = ds.field_list + [("gas", "mass"), ("gas", "radial_velocity"), ("gas", "tangential_velocity"), ("gas", "velocity_magnitude")]
    prof = yt.create_profile(ds.all_data(), 'radius', prof_fields, n_bins=600, extrema={'radius':(0, 8e7)}, logs={'radius':False})

    #save as csv for easy access later
    df = prof.to_dataframe(include_std=True).to_csv(f"profiles/{ds.basename}_all_fields.csv")
    
if __name__ == '__main__':
    #load and add mass field
    fname = argv[1]
    save_profile(fname)
