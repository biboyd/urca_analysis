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

@yt.derived_field(name=("gas", "velocity_magnitude"), sampling_type='local')
def velocity_magnitude(field, data):
    return (data["velx"] ** 2.0 + data["vely"] ** 2.0 + data["velz"] ** 2.0)**0.5
    
#@yt.derived_field(name=("gas", "mass"), sampling_type='local')
#def _mass(field, data):
#    return data["boxlib", "rho"]*unyt.g/unyt.cm**3 *data['boxlib', 'volume']

    
def save_profile(fname):
    ds = yt.load(fname)

    
    @yt.derived_field(name=("gas", "radial_velocity"), sampling_type='local')
    def velocity_radial(field, data):
        return (data["velx"] * (data[('index', 'x')] - ds.domain_center[0]) + 
                data["vely"] *  (data[('index', 'y')] - ds.domain_center[1]) + 
                data["velz"]  *  (data[('index', 'z')] - ds.domain_center[2]))/data[('index', 'radius')]
        
    @yt.derived_field(name=("gas", "tangential_velocity"), sampling_type='local')
    def velocity_tangential(field, data):
        return (data["velocity_magnitude"]**2.0 - data["gas", "radial_velocity"]**2.0)**0.5
    
    #generate profiles
    
    additional_fields = [("gas", "mass"),
                        ("gas", "radial_velocity"),
                        ("gas", "tangential_velocity"),
                        ("gas", "velocity_magnitude"),
                        ("gas", "velocity_x_squared"),
                        ("gas", "velocity_y_squared"),
                        ("gas", "velocity_z_squared")]

    prof_fields = ds.field_list + additional_fields
    prof = yt.create_profile(ds.all_data(), 'radius', prof_fields, n_bins=100, extrema={'radius':(0, 1e8)}, logs={'radius':False})

    #save as csv for easy access later
    df = prof.to_dataframe(include_std=True).to_csv(f"profiles/{ds.basename}_all_fields.csv")
    
if __name__ == '__main__':
    #load and add mass field
    fname = argv[1]
    save_profile(fname)

