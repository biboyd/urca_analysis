import yt
import unyt
import pandas as pd
from sys import argv

@yt.derived_field(name=("gas", "velocity_x_squared"), sampling_type='local')
def velocity_x_squared(field, data):
    return data["velx"] ** 2.0

@yt.derived_field(name=("gas", "velocity_y_squared"), sampling_type='local')
def velocity_y_squared(field, data):
    return data["vely"] ** 2.0

@yt.derived_field(name=("gas", "velocity_z_squared"), sampling_type='local')
def velocity_z_squared(field, data):
    return data["velz"] ** 2.0

def save_profile(fname):
    ds = yt.load(fname)
    #generate profiles
    
    additional_fields = [("gas", "mass"),
                        ("gas", "radial_velocity"),
                        ("gas", "tangential_velocity"),
                        ("gas", "velocity_magnitude"),
                        ("gas", "velocity_x_squared"),
                        ("gas", "velocity_y_squared"),
                        ("gas", "velocity_z_squared")]

    prof_fields = ds.field_list + additional_fields
    prof = yt.create_profile(ds.all_data(), 'radius', prof_fields, n_bins=320, extrema={'radius':(0, 8e7)}, logs={'radius':False})

    #save as csv for easy access later
    df = prof.to_dataframe(include_std=True).to_csv(f"profiles/{ds.basename}_all_fields.csv")
    
if __name__ == '__main__':
    #load and add mass field
    fname = argv[1]
    save_profile(fname)

