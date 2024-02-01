import yt
import unyt
import pandas as pd
from sys import argv

def _mass(field, data):
    return data["boxlib", "rho"]*unyt.g/unyt.cm**3 *data['boxlib', 'volume']

def _c12_complement(field, data):
    return (0.39975 - data["boxlib", "X(c12)"])

def _A23_frac(field, data):

    #A=23 nuc
    A23_frac = (data['boxlib', 'X(ne23)'] - data['boxlib', 'X(na23)'])/ (data['boxlib', 'X(ne23)'] + data['boxlib', 'X(na23)'])
    
    return A23_frac

def _Ye23(field, data):
    # sum 1/2
    Ye=0.5*(data['boxlib', 'X(c12)']+data['boxlib', 'X(o16)']+data['boxlib', 'X(he4)']+data['boxlib', 'X(ne20)'])
    
    #sum ones
    Ye+=(data['boxlib', 'X(h1)'])
    
    #A=23 nuc
    Ye+= (10.*data['boxlib', 'X(ne23)'] + 11.*data['boxlib', 'X(na23)'] + 12.*data['boxlib', 'X(mg23)'])/23.
    
    return Ye

def _Ye21_23_25(field, data):
    # sum 1/2
    Ye=0.5*(data['boxlib', 'X(c12)']+data['boxlib', 'X(o16)']+data['boxlib', 'X(he4)']+data['boxlib', 'X(ne20)'])
    
    #sum ones
    Ye+=(data['boxlib', 'X(h1)'])
    
    #A=21 nuc
    Ye+= (9.*data['boxlib', 'X(f21)'] + 10.*data['boxlib', 'X(ne21)'])/21.
        
    #A=23 nuc
    Ye+= (10.*data['boxlib', 'X(ne23)'] + 11.*data['boxlib', 'X(na23)'] + 12.*data['boxlib', 'X(mg23)'])/23.
        
    #A=25 nuc
    Ye+= (11.*data['boxlib', 'X(na25)'] + 12.*data['boxlib', 'X(mg25)'])/25.
    
    return Ye

def _Ye_asymmetry(field, data):
    return data['boxlib', 'Ye'] -0.5

def _Tbar(field, data):
    return -1.*(data['boxlib', 'tfromp'] - data['boxlib', 'tpert'])

# mean molecular weight
def _mu23(field, data):
    muinv = 2*data['boxlib', 'X(h1)'] + data['boxlib', 'X(n)'] + 3./2. * data['boxlib', 'X(he4)'] + 7./12. * data['boxlib', 'X(c12)'] +  9./16. * data['boxlib', 'X(o16)'] + 11./20. * data['boxlib', 'X(ne20)'] + 11./23. * data['boxlib', 'X(ne23)'] + 12./23. * data['boxlib', 'X(na23)'] + 13./23. * data['boxlib', 'X(mg23)']
    
    return 1./(muinv)

# mean molecular weight
def _mu21_23_25(field, data):
    muinv = 2*data['boxlib', 'X(h1)'] + data['boxlib', 'X(n)'] + 3./2. * data['boxlib', 'X(he4)'] + 7./12. * data['boxlib', 'X(c12)'] +  9./16. * data['boxlib', 'X(o16)'] + 11./20. * data['boxlib', 'X(ne20)'] + 11./21. * data['boxlib', 'X(ne21)'] + 11./23. * data['boxlib', 'X(ne23)'] + 10./21. * data['boxlib', 'X(f21)'] + 12./23. * data['boxlib', 'X(na23)'] + 12./25. * data['boxlib', 'X(na25)'] + 13./23. * data['boxlib', 'X(mg23)'] + 13./25. * data['boxlib', 'X(mg25)']
    
    return 1./(muinv)


def save_profile(fname):
    ds = yt.load(fname, hint='amrex')
    # add additional fields
    ds.add_field(name=("gas", "mass"),
                function=_mass,
                take_log=False,
                display_name="$m$ ",
                units='g',
                sampling_type="local")
    """
    if ("boxlib", "X(ne21)") in ds.field_list:
        ds.add_field(
            name=("boxlib", "Ye"),
            function=_Ye21_23_25,
            take_log=False,
            sampling_type="local")
        ds.add_field(name=("boxlib", "mu"),
            function=_mu21_23_25,
            take_log=False,
            sampling_type="local")
    else:
        ds.add_field(
            name=("boxlib", "Ye"),
            function=_Ye23,
            take_log=False,
            sampling_type="local")
        ds.add_field(name=("boxlib", "mu"),
            function=_mu23,
            take_log=False,
            sampling_type="local")
    """
    def _radvel(field, data):
        return (data[('boxlib', 'velx')] * (data[('index', 'x')] - ds.domain_center[0]) + 
                data[('boxlib', 'vely')] *  (data[('index', 'y')] - ds.domain_center[1]) + 
                data[('boxlib', 'velz')]  *  (data[('index', 'z')] - ds.domain_center[2]))/data[('index', 'radius')]
    if ('boxlib', 'radial_velocity') not in ds.field_list:
        ds.add_field(name=('boxlib', 'radial_velocity'),
                    function=_radvel,
                    take_log=True,
                    sampling_type='local')
    #generate profiles
    prof_fields=[('boxlib', 'X(ne23)'), ('boxlib', 'X(na23)'), ('boxlib', 'X(c12)'), ('boxlib', 'X(o16)'), ('boxlib', 'tfromp'), ('boxlib', 'vort'), ('boxlib', 'radial_velocity'), ('boxlib', 'rho'), ('boxlib', 'p0pluspi'), ('boxlib', 'Hnuc'), ('boxlib', 'entropy')]
    
    prof = yt.create_profile(ds.all_data(), 'radius', prof_fields, n_bins=500, extrema={'radius':(0, 8e7)}, logs={'radius':False})

    #save as csv for easy access later
    df = prof.to_dataframe().to_csv(f"profiles/{ds.basename}.csv")
    
if __name__ == '__main__':
    #load and add mass field
    fname = argv[1]
    save_profile(fname)
