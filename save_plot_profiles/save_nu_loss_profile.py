import yt
import unyt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv


# define fields
def _mass(field, data):
    return data[('index', 'volume')]*data[('boxlib', 'rho')]* unyt.g/unyt.cm**3


def _tot_nu_loss(field, data):
    return data['specific_nu_loss'] * data[('gas', 'mass')]


def add_nuloss_fields(ds_nu, nu_field_list):
    spec_fields = ['specific_nu_loss']
    raw_fields = ['tot_nu_loss']

    # define our nu loss
    def _spec_nu_loss(field, data):
        tot = 0.
        for nuloss in nu_field_list:
            tot += -1.*data[nuloss]
        tot += data[('thermal_nu_loss')]

        return tot * unyt.erg/unyt.g/unyt.s

    # add mass and totals
    ds_nu.add_field(name=('gas', 'mass'), function=_mass, units='g', sampling_type='local')
    ds_nu.add_field(name=('gas', 'specific_nu_loss'), function=_spec_nu_loss, units='erg/g/s', sampling_type='local')
    ds_nu.add_field(name=('gas', 'tot_nu_loss_rate'), function=_tot_nu_loss, units='erg/s', sampling_type='local')

    # add individual mass multi fields
    for fld in nu_field_list:

        ds_nu.add_field(name=('gas', 'tot_{fld}'),
                        function=lambda field, data: -1 * data[fld] * data[('gas', 'mass')] * unyt.erg/unyt.g/unyt.s,
                        units='erg/s', sampling_type='local')

        raw_fields.append(f"tot_{fld}")

    return spec_fields, raw_fields


def save_nuloss_profile(ds):
    # grab relevant nuloss fields
    nu_field_list = []
    for fld in ds.field_list:
        # check field starts 'A' ends in 'nu_loss'
        name = fld[1]
        if 'A' == name[0] and '_nu_loss' == name[-8:]:
            nu_field_list.append(fld)

    # add fields
    spec_fields, raw_fields = add_nuloss_fields(ds, nu_field_list)

    # calc and save energy rate profile
    N_bins = 100
    fields = ds.field_list + spec_fields
    prof = yt.create_profile(ds.all_data(), 'radius', fields, logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})
    prof_sum = yt.create_profile(ds.all_data(), 'radius', raw_fields, weight=None, logs={'radius':False}, n_bins=N_bins, extrema={'radius':(0, 1e8)})

    # save profiles
    df = prof.to_dataframe()
    df.to_csv(f"nuloss_profiles/{ds.basename}_avg_profiles.csv")

    df_sum = prof_sum.to_dataframe()
    df_sum.to_csv(f"nuloss_profiles/{ds.basename}_sum_profiles.csv")


if __init__ == "__main__":
    fname = argv[1]
    ds = yt.load(fname)
    save_nuloss_profile(ds)
