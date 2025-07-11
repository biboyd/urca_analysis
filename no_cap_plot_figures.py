import argparse
import yt
import unyt as u
#import pandas as pd
import numpy as np

# pass arguments
parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('field', type=str, help='field of slice to plot')
parser.add_argument('-w', '--width', type=float, help='width of slice in km.', default=1.2e3)
parser.add_argument('-g', '--uselog', type=str, help='whether to use log/symlog. defaults to whatever yt decides. if want to set declare True/False', default=None)
parser.add_argument('-t', '--linthresh', type=float, help='lin threshold to use in a symlog. defaults to no threshold', default=None)
parser.add_argument('-l', '--zlo', type=str, help='lower bound of plot. Default depends on field. genrally min', default=None)
parser.add_argument('-u', '--zup', type=str, help='upper bound of plot. Default depends on field. gnerally max', default=None)
parser.add_argument('-c', '--cmap', type=str, help='colormap of plot. Default depends on field', default=None)
parser.add_argument('-o', '--outdir', type=str, help='out directory of plot. Default depends on field', default=None)
parser.add_argument('-pg', '--plotgrid', help='plot gridlines', default=False, action='store_true')
parser.add_argument('-a', '--axis', type=str, help='x, y, or z axis to take slice of. defaults to x', default="z")
parser.add_argument('-flip', '--flip_axis', default=False, action='store_true', help='flip vertical axis of slice')
parser.add_argument('-cH', '--contour_Hnuc0', help='plot Hnuc=0 contours', default=False, action='store_true')
parser.add_argument('-cR', '--contour_rho', help='plot density contours', default=None, type=float, nargs="*")
parser.add_argument('-sph', '--sphere', help='plot sphere at radius (in km)', default=None, type=float, nargs="*")
parser.add_argument('-cU21', '--contour_Urca21', help='plot A21_frac=0 contours', default=False, action='store_true')
parser.add_argument('-cU23', '--contour_Urca23', help='plot A23_frac=0 contours', default=False, action='store_true')
parser.add_argument('-cU25', '--contour_Urca25', help='plot A25_frac=0 contours', default=False, action='store_true')
parser.add_argument('-mR', '--mixingregion', help='plot mixing region bounds', default=False, action='store_true')
parser.add_argument('-s', '--streamlines', help='plot streamlines', default=False, action='store_true')
parser.add_argument('-arr', '--arrows', help='plot velocity arrows', default=False, action='store_true')
parser.add_argument('-lc', '--linecolor', type=str, help='color for streamlines', default='black')
parser.add_argument('-cc', '--contourcolor', type=str, help='contour color', default='black')
parser.add_argument('-z', '--conv_zone', help='plot just convection zone', default=False, action='store_true')
parser.add_argument('-fs', '--fontsize', type=int, help='fontsize for plot', default=None)
parser.add_argument('-at', '--annotate_text', type=str, help='text to annotate', default=None)
parser.add_argument('-ap', '--annotate_position', type=float, nargs=2, help='postion to place annotated text', default=(0.,0.))


args = parser.parse_args()

# default dicts

default_uselog = {"magvel" : None,
                  "radial_velocity" : True,
                  "Hnuc" : None,
                  "MachNumber" : None,
                  "tfromp" : False,
                  "c12_complement" : False,
                  "X(Ne23)" : True,
                  "X(Na23)" : True,
                  "Ye" : False,
                  "Ye_asymmetry" : False,
                  "mu" : False,
                  "vort" : False,
                  "tpert" : False

                 }

default_linthresh = {"magvel" : None,
                     "radial_velocity" : 1e-1,
                     "MachNumber" : None,
                     "Hnuc" : 1e5,
                     "tfromp" : None,
                     "c12_complement" : None,
                     "X(Ne23)" : None,
                     "X(Na23)" : None,
                     "Ye" : None
                 }
                                  
default_zlim = {"magvel" : [1e-1, 1e2],
                "radial_velocity" : [-2e2, 2e2],
                "Hnuc" : [-2e11, 2e11],
                "MachNumber" : ['min', 'max'],
                "tfromp" : ['min', 'max'],
                "c12_complement" : ['min', 'max'],
                "Ye_asymmetry" : ['min', 'max'],
                "mu" : ['min', 'max'],
                "A21_frac" : [-1., 1.],
                "A23_frac" : [-1., 1.],
                "A25_frac" : [-1., 1.],
                "ad_excess" : [-0.5, 0.5],
                "ad_excess_led" : [-0.5, 0.5],
                "vort" : [0., 6.],
                "tpert" : [-5e5, 5e5]
                }

default_cmap = {"magvel" : "cividis",
                "radial_velocity" : "RdBu_r",
                "Hnuc" : 'PiYG',
                "MachNumber" : 'cividis',
                "tfromp" : 'magma',
                "c12_complement" : 'viridis',
                "X(Ne23)" : None,
                "X(Na23)" : None,
                "Ye_asymmetry" : 'viridis',
                "rho_Ye" : "plasma",
                "mu" : 'magma',
                "eta" : 'plasma',
                "A21_frac" : 'PuOr',
                "A23_frac" : 'PuOr',
                "A25_frac" : 'PuOr',
                "ad_excess" : "RdBu",
                "ad_excess_led" : "RdBu",
                "tpert" : "magma",
                "tot_nu_loss" : "inferno"
                }

default_outdir = {"magvel": "plots_magvel/",
                  "radial_velocity" : "plots_radvel/",
                  "Hnuc" : "plots_Hnuc/",
                  "MachNumber" : "plots_machNumber/",
                  "tfromp" : "plots_temp/",
                  "c12_complement" : "plots_c12_complement/",
                  "X(C12)" : "plots_c12_frac/",
                  "X(C13)" : "plots_c13_frac/",
                  "X(N13)" : "plots_n13_frac/",
                  "X(O16)" : "plots_o16_frac/",
                  "X(Ne20)" : "plots_ne20_frac/",
                  "X(Ne21)" : "plots_ne21_frac/",
                  "X(F21)" : "plots_f21_frac/",
                  "X(Ne23)" : "plots_ne23_frac/",
                  "X(Na23)" : "plots_na23_frac/",
                  "X(Mg25)" : "plots_mg25_frac/",
                  "X(Na25)" : "plots_na25_frac/",
                  "omegadot(C12)"  : "plots_omegadot_c12/",
                  "omegadot(O16)"  : "plots_omegadot_o16/",
                  "omegadot(Ne21)" : "plots_omegadot_ne21/",
                  "omegadot(F21)"  : "plots_omegadot_f21/",
                  "omegadot(Ne23)" : "plots_omegadot_ne23/",
                  "omegadot(Na23)" : "plots_omegadot_na23/",
                  "omegadot(Mg25)" : "plots_omegadot_mg25/",
                  "omegadot(Na25)" : "plots_omegadot_na25/",
                  "rhopert" : "plots_rhopert/",
                  "Pi" : "plots_Pi/",
                  "rho" : "plots_density/",
                  "Ye" : "plots_Ye/",
                  "Ye_asymmetry" : "plots_Ye_asymmetry/",
                  "rho_Ye" : "plots_rho_Ye/",
                  "mu" : "plots_mu/",
                  "eta" : "plots_eta/",
                  "A21_frac" : "plots_A21_frac/",
                  "A23_frac" : "plots_A23_frac/",
                  "A25_frac" : "plots_A25_frac/",
                  "A21_ratio" : "plots_A21_ratio/",
                  "A23_ratio" : "plots_A23_ratio/",
                  "A25_ratio" : "plots_A25_ratio/",
                  "vort" : "plots_vorticity/",
                  "ad_excess" : "plots_ad_excess/",
                  "ad_excess_led" : "plots_ad_excess_led/",
                  "tpert" : "plots_tpert/",
                  "tot_nu_loss" : "plots_nu_loss/"
                 }

# derived field functions

def _UrcaActivity(field, data):
    #A=23 nuc
    return data['boxlib', 'omegadot(Ne23)']

def _c12_complement(field, data):
    return (1.0 - data["boxlib", "X(C12)"])

def _A21_frac(field, data):

    #A=21 nuc
    return (data['boxlib', 'X(F21)'] - data['boxlib', 'X(Ne21)'])/ (data['boxlib', 'X(F21)'] + data['boxlib', 'X(Ne21)'])

def _A21_ratio(field, data):

    #A=23 nuc
    return data['boxlib', 'X(F21)']/data['boxlib', 'X(Ne21)']
   
def _A23_ratio(field, data):

    #A=23 nuc
    return data['boxlib', 'X(Ne23)']/data['boxlib', 'X(Na23)']

def _A23_frac(field, data):

    #A=23 nuc
    return(data['boxlib', 'X(Ne23)'] - data['boxlib', 'X(Na23)'])/ (data['boxlib', 'X(Ne23)'] + data['boxlib', 'X(Na23)'])
     
def _A25_ratio(field, data):

    #A=25 nuc
    return data['boxlib', 'X(Na25)'] / data['boxlib', 'X(Mg25)']
   
def _A25_frac(field, data):

    #A=25 nuc
    return (data['boxlib', 'X(Na25)'] - data['boxlib', 'X(Mg25)'])/ (data['boxlib', 'X(Na25)'] + data['boxlib', 'X(Mg25)'])

def _A23_tot(field, data):

    #A=23 nuc
    return(data['boxlib', 'X(Ne23)'] + data['boxlib', 'X(Na23)'])

def _mytan_vel(field, data):
        return np.sqrt(data[('boxlib', 'magvel')]**2 -  data[('boxlib', 'radial_velocity')]**2)
# electron fraction if just A=23
def _Ye23(field, data):
    # sum 1/2
    Ye=0.5*(data['boxlib', 'X(c12)']+data['boxlib', 'X(o16)']+data['boxlib', 'X(he4)']+data['boxlib', 'X(ne20)'])
    
    #sum ones
    Ye+=(data['boxlib', 'X(h1)'])
    
    #A=23 nuc
    Ye+= (10.*data['boxlib', 'X(ne23)'] + 11.*data['boxlib', 'X(na23)'] + 12.*data['boxlib', 'X(mg23)'])/23.
    
    return Ye

# electron fraction with A=21,23,25 nuclei
def _Ye21_23_25(field, data):
    # sum 1/2
    Ye=0.5*(data['boxlib', 'X(C12)']+data['boxlib', 'X(O16)']+data['boxlib', 'X(He4)']+data['boxlib', 'X(Ne20)'])
    
    #sum ones
    Ye+=(data['boxlib', 'X(H1)'])
    
    #A=21 nuc
    Ye+= (9.*data['boxlib', 'X(F21)'] + 10.*data['boxlib', 'X(Ne21)'])/21.
        
    #A=23 nuc
    Ye+= (10.*data['boxlib', 'X(Ne23)'] + 11.*data['boxlib', 'X(Na23)'] + 12.*data['boxlib', 'X(Mg23)'])/23.
        
    #A=25 nuc
    Ye+= (11.*data['boxlib', 'X(Na25)'] + 12.*data['boxlib', 'X(Mg25)'])/25.
    
    return Ye

def _Ye_asymmetry(field, data):
    return data['boxlib', 'Ye'] - 0.5

def _rho_Ye(field, data):
    return data['boxlib', 'rho'] * data['boxlib', 'Ye']

# neutronization
def _eta(field, data):
    return 1- 2*data['boxlib', 'Ye']

def _Tbar(field, data):
    return -1.*(data['boxlib', 'tfromp'] - data['boxlib', 'tpert'])

# mean molecular weight with only A=23 nuclei
def _mu23(field, data):
    muinv = 2*data['boxlib', 'X(H1)'] + data['boxlib', 'X(N)'] + 3./2. * data['boxlib', 'X(He4)'] + 7./12. * data['boxlib', 'X(C12)'] +  9./16. * data['boxlib', 'X(O16)'] + 11./20. * data['boxlib', 'X(Ne20)'] + 11./23. * data['boxlib', 'X(Ne23)'] + 12./23. * data['boxlib', 'X(Na23)'] + 13./23. * data['boxlib', 'X(Mg23)']
    
    return 1./(muinv)

# mean molecular weight with A=21,23,25 nuclei
def _mu21_23_25(field, data):
    muinv = 2*data['boxlib', 'X(H1)'] + data['boxlib', 'X(N)'] + 3./2. * data['boxlib', 'X(He4)'] + 7./12. * data['boxlib', 'X(C12)'] +  9./16. * data['boxlib', 'X(O16)'] + 11./20. * data['boxlib', 'X(Ne20)'] + 11./21. * data['boxlib', 'X(Ne21)'] + 11./23. * data['boxlib', 'X(Ne23)'] + 10./21. * data['boxlib', 'X(F21)'] + 12./23. * data['boxlib', 'X(Na23)'] + 12./25. * data['boxlib', 'X(Na25)'] + 13./23. * data['boxlib', 'X(Mg23)'] + 13./25. * data['boxlib', 'X(Mg25)']
    
    return 1./(muinv)

def _mass(field, data):
    return data[('boxlib', 'rho')] * data[('gas', 'volume')] * u.g/u.cm**3

# convection related fields
def _ad_excess(field, data):
    return  data[('boxlib', 'conv_actual')] - data[('boxlib', 'conv_adiabatic')]

def _ad_excess_led(field, data):
    return  data[('boxlib', 'conv_actual')] - data[('boxlib', 'conv_ledoux')]

# nu loss rate stuff
def _invert_ecap_losses(field, data):
    return -1.0 * data[ 'A23_electron_capture_nu_loss'] * u.erg/u.g/u.s

def _invert_beta_losses(field, data):
    return -1.0 * data[ 'A23_beta_decay_nu_loss'] * u.erg/u.g/u.s

def _tot_nu_loss(field, data):
    return  (data['thermal_nu_loss'] -1.0 * data[ 'A23_electron_capture_nu_loss']  -1.0 * data[ 'A23_beta_decay_nu_loss'] )* u.erg/u.g/u.s

def _rel_tpert(field, data):
    return data[('boxlib', 'tpert')]/data[('boxlib', 'tfromp')]

def _w0_magnitude(field, data):
    return np.sqrt(data[('boxlib', 'w0x')]**2 + data[('boxlib', 'w0y')]**2 + data[('boxlib', 'w0z')]**2) * u.cm/u.s

def plot_slice(ds, slice_field, args):

    width, uselog, linthresh, zlo, zup, cmap, outdir = load_defaults(slice_field, args)
    field = slice_field

    #treat a21/23/25 differently b/c of possible contours
    if slice_field == "A21_frac" or args.contour_Urca21:
        ds.add_field(
            name=("boxlib", "A21_frac"),
            function=_A21_frac,
            take_log=False,
            dimensions='dimensionless',
            display_name="$\\frac{\\mathrm{X({}^{21}F)} - \\mathrm{X({}^{21}Ne)}}{\\mathrm{{}^{21}X(F)} + \\mathrm{{}^{21}X(Ne)}}$",
            sampling_type="local")
        
    if slice_field == "A23_frac" or args.contour_Urca23:
        ds.add_field(
            name=("boxlib", "A23_frac"),
            function=_A23_frac,
            take_log=False,
            dimensions='dimensionless',
            display_name="$\\frac{\\mathrm{X({}^{23}Ne)} - \\mathrm{X({}^{23}Na)}}{\\mathrm{{}^{23}X(Ne)} + \\mathrm{{}^{23}X(Na)}}$",
            sampling_type="local")

    if slice_field == "w0_magnitude":
        ds.add_field(
            name=("boxlib", "w0_magnitude"),
            function=_w0_magnitude,
            take_log=False,
            units='km/s',
            display_name="$ | w_0 | $",
            sampling_type="local")
        
        
    if slice_field == "A25_frac" or args.contour_Urca25:
        ds.add_field(
            name=("boxlib", "A25_frac"),
            function=_A25_frac,
            take_log=False,
            dimensions='dimensionless',
            display_name="$\\frac{\\mathrm{X({}^{25}Na)} - \\mathrm{X({}^{25}Mg)}}{\\mathrm{{}^{25}X(Na)} + \\mathrm{{}^{25}X(Mg)}}$",
            sampling_type="local")
        
    if slice_field == "A23_tot":
        ds.add_field(
            name=("boxlib", "A23_tot"),
            function=_A23_tot,
            take_log=False,
            dimensions='dimensionless',
            display_name="$\\mathrm{X({}^{23}Ne)} + \\mathrm{X({}^{23}Na)}$",
            sampling_type="local")
         
    if slice_field == "A21_ratio":
        ds.add_field(
            name=("boxlib", "A21_ratio"),
            function=_A21_ratio,
            take_log=False,
            dimensions='dimensionless',
            display_name="$\\mathrm{X({}^{21}F)} / \\mathrm{X({}^{21}Ne)}$",
            sampling_type="local")

    if slice_field == "A23_ratio":
        ds.add_field(
            name=("boxlib", "A23_ratio"),
            function=_A23_ratio,
            take_log=False,
            dimensions='dimensionless',
            display_name="$\\mathrm{X({}^{23}Ne)} / \\mathrm{X({}^{23}Na)}$",
            sampling_type="local")
        
    if slice_field == "A25_ratio":
        ds.add_field(
            name=("boxlib", "A25_ratio"),
            function=_A25_ratio,
            take_log=False,
            dimensions='dimensionless',
            display_name="$\\mathrm{X({}^{25}Na)} / \\mathrm{X({}^{25}Mg)}$",
            sampling_type="local")

    if slice_field == "UrcaActivity":
        ds.add_field(
            name=("boxlib", "UrcaActivity"),
            function=_UrcaActivity,
            take_log=False,
            sampling_type="local")
        
    #add fields to dataset if needed
    if slice_field == "c12_complement":
        ds.add_field(
            name=("boxlib", "c12_complement"),
            function=_c12_complement,
            take_log=False,
            units = "dimensionless",
            display_name="$1 - \\mathrm{X}({}^{12}\\mathrm{C})$ ",
            sampling_type="local")
    if slice_field == "mytan_vel":
        ds.add_field(('boxlib', 'mytan_vel'), _mytan_vel, units='km/s', sampling_type='local', force_override=True)

    # check if include A=21 urca pair in there. include those in calc.
    elif slice_field == "Ye" or slice_field == "Ye_asymmetry" or slice_field == "eta" or slice_field == "rho_Ye":
        if ("boxlib", "X(ne21)") in ds.field_list:
            ds.add_field(
                name=("boxlib", "Ye"),
                function=_Ye21_23_25,
                units = "dimensionless",
                take_log=False,
                sampling_type="local")
        else:
            ds.add_field(
                name=("boxlib", "Ye"),
                function=_Ye23,
                take_log=False,
                units = "dimensionless",
                sampling_type="local")
    
        ds.add_field(
            name=("boxlib", "Ye_asymmetry"),
            function=_Ye_asymmetry,
            take_log=False,
            display_name="$\\mathrm{Y_e} - 0.5$",
            units = "dimensionless",
            sampling_type="local")
        
        ds.add_field(
            name=("boxlib", "rho_Ye"),
            function=_rho_Ye,
            take_log=True,
            display_name="$\\rho \\mathrm{Y_e}$",
            units = "dimensionless",
            sampling_type="local")
        
        ds.add_field(
            name=("boxlib", "eta"),
            function=_eta,
            take_log=False,
            units = "dimensionless",
            display_name="$\\eta$",
            sampling_type="local")
        
    elif slice_field == "Tbar":
        ds.add_field(
        name=("boxlib", "Tbar"),
        function=_Tbar,
        take_log=False,
        sampling_type="local")

    elif slice_field == "mu":
        if ("boxlib", "X(ne21)") in ds.field_list:
            ds.add_field(
                name=("boxlib", "mu"),
                function=_mu21_23_25,
                take_log=False,
                display_name="$\\mu$ ",
                units = "dimensionless",
                sampling_type="local")
        else:
            ds.add_field(
                name=("boxlib", "mu"),
                function=_mu23,
                take_log=False,
                display_name="$\\mu$ ",
                units = "dimensionless",
                sampling_type="local")
        
    elif slice_field == "ad_excess":
        ds.add_field(name=("boxlib", "ad_excess"),
            function=_ad_excess,
            take_log=False,
            display_name="Schwarzschild $\\Delta \\nabla $",
            sampling_type="local")
        
    elif slice_field == "ad_excess_led":
        ds.add_field(name=("boxlib", "ad_excess_led"),
            function=_ad_excess_led,
            take_log=False,
            display_name="Ledoux $\\Delta \\nabla$",
            sampling_type="local")
        
    elif slice_field == "nu_ecap_losses":
        ds.add_field(
            name=("boxlib", "nu_ecap_losses"),
            function=_invert_ecap_losses,
            take_log=True,
            units='erg/g/s',
            display_name="$\\nu_{ecap}$ energy losses",
            sampling_type="local", force_override=True)

    elif slice_field == "nu_beta_losses":
        ds.add_field(
            name=("boxlib", "nu_beta_losses"),
            function=_invert_beta_losses,
            take_log=True,
            units='erg/g/s',
            display_name="$\\nu_{beta}$ energy losses",
            sampling_type="local", force_override=True)

    elif slice_field == "tot_nu_loss":
        ds.add_field(
            name=("boxlib", "tot_nu_loss"),
            function=_tot_nu_loss,
            take_log=True,
            units='erg/g/s',
            display_name="Energy Loss to Neutrino Emissions",
            sampling_type="local", force_override=True)

    elif slice_field == "rel_tpert":
        ds.add_field(
            name=("boxlib", "rel_tpert"),
            function=_rel_tpert,
            take_log=False,
            display_name="$\\frac{T - \\bar{T}}{T}$",
            sampling_type="local", force_override=True)
        
    if args.axis in ('x', 'y', 'z'):
        if args.conv_zone:
            zone = np.load("conv_zone_over_time.npy")
            r = zone[np.argmin(np.abs(ds.current_time.value - zone[:, 0])), 1]
            dat_source = ds.sphere(ds.domain_center, (r, 'km'))
        else:
            dat_source = ds.all_data()
        
        s = yt.SlicePlot(ds, args.axis, field, width = width, data_source=dat_source)
        s.flip_vertical()
    else:
        raise ValueError(f"axis argument given ({args.a}) invalid. must be 'x', 'y', or 'z' ")
    
    if args.flip_axis:
        s.flip_vertical()

    # set velocity to default as km/s
    if field == ('boxlib', 'radial_velocity') or field == ('boxlib', 'magvel') or field == 'radial_velocity' or field == 'velocity_magnitude':
        s.set_unit(field, "km/s")
    
    s.set_zlim(field, zlo, zup)

    #set colorbar and scale
    if uselog is not None or linthresh is not None:
        s.set_log(field, uselog, linthresh=linthresh)
    
    if cmap is not None:
        s.set_cmap(field, cmap)
    
    # field specific changes. Mostly lable/units related
    if field == ('boxlib', 'Hnuc'):
        s.set_colorbar_label(field, "$\\dot{\\epsilon}_{\\mathrm{nuc}} \\; (\\frac{\\mathrm{erg}}{g \\cdot s})$")
        

    if args.plotgrid:
        s.annotate_grids()
    
    if args.contour_Hnuc0:
        s.annotate_contour(('boxlib', 'Hnuc'), levels=1, factor=1, take_log=False,
                       clim=(0.0,0.0), plot_args={'colors' : args.contourcolor})
    
    if args.contour_rho is not None:
        for curr_dens in args.contour_rho:
            s.annotate_contour(('boxlib', 'rho'), levels=1, factor=1, take_log=False,
                       clim=(curr_dens, curr_dens), plot_args={'colors' : args.contourcolor})
    
    if args.sphere is not None:
        for anno_sphere in args.sphere:
            s.annotate_sphere(ds.domain_center, (anno_sphere, 'km'))
        
    if args.contour_Urca21:
        s.annotate_contour(('boxlib', 'A21_frac'), levels=1, factor=1, take_log=False,
                       clim=(0.0,0.0), plot_args={'colors' : args.contourcolor, 'linestyles' : '--'})
        
    if args.contour_Urca23:
        s.annotate_contour(('boxlib', 'A23_frac'), levels=1, factor=1, take_log=False,
                       clim=(0.0,0.0), plot_args={'colors' : args.contourcolor, 'linestyles' : '--'})
        
    if args.contour_Urca25:
        s.annotate_contour(('boxlib', 'A25_frac'), levels=1, factor=1, take_log=False,
                       clim=(0.0,0.0), plot_args={'colors' : args.contourcolor, 'linestyles' : '--'})
    
    if args.mixingregion:
        mix_inner, mix_outer = get_mixing_region(ds)
        
        s.annotate_sphere(ds.domain_center, mix_inner, circle_args={"color": "black", "linestyle" : "-."})
        s.annotate_sphere(ds.domain_center, mix_outer, circle_args={"color": "black", "linestyle" : "-."})
        print(mix_inner)
        
    if args.streamlines or args.arrows:
        if args.axis == 'x':
            vel_fld1 = ('boxlib', 'vely')
            vel_fld2 = ('boxlib', 'velz')
        elif args.axis == 'y':
            vel_fld1 = ('boxlib', 'velz')
            vel_fld2 = ('boxlib', 'velx')
        elif args.axis == 'z':
            vel_fld1 = ('boxlib', 'velx')
            vel_fld2 = ('boxlib', 'vely')
            
        if args.streamlines:
            s.annotate_streamlines(vel_fld1, vel_fld2,density=0.8, color=args.linecolor)
        elif args.arrows:
            s.annotate_quiver(vel_fld1, vel_fld2, color=args.linecolor)
        
    # annotate text
    if args.annotate_text is not None:
        s.annotate_text(args.annotate_position, args.annotate_text, coord_system='plot')
    # edit the fontsize
    if args.fontsize is not None:
        s.set_font_size(args.fontsize)
        
    # add simulation time and save
    s.annotate_timestamp(draw_inset_box=True)
    return s, outdir 

def load_defaults(slice_field, args):

    # everything default to None unless specified or included in default dictionaries
    width = (args.width, 'km')
    
    # alias the velocit magnitude field for same defaults
    if slice_field == 'velocity_magnitude':
        slice_field = 'magvel'

    # set log defaults
    if args.uselog is None and slice_field in default_uselog.keys():
        uselog = default_uselog[slice_field]
    else:
        if args.uselog is None:
            uselog = None
        elif args.uselog == 'True':
            uselog = True
        elif args.uselog == 'False':
            uselog = False
        else:
            raise ValueError(f"uselog must be set to either True or False. {args.uselog} if not accepted")
    
    # if linthresh is set then also log must be set to true
    if uselog is False:
        linthresh = None #set to None b/c linthresh overrides uselog in yt
    elif args.linthresh is None and slice_field in default_linthresh.keys():
        linthresh = default_linthresh[slice_field]
    else:
        linthresh = args.linthresh
    
    # set lower/upper bound. set to min/max if not defined otherwise
    if args.zlo is None:
        if slice_field in default_zlim.keys():
            zlo = default_zlim[slice_field][0]
        else:
            zlo = 'min'
    elif args.zlo == 'min':
        zlo = 'min'
    else:
        zlo = float(args.zlo)

    if args.zup is None:
        if slice_field in default_zlim.keys():
            zup = default_zlim[slice_field][1]
        else:
            zup = 'max'
    elif args.zup == 'max':
        zup = 'max'
    else:
        zup = float(args.zup)
    
    # set colormap
    if args.cmap is None and slice_field in default_cmap.keys():
        cmap = default_cmap[slice_field]
    else:
        cmap = args.cmap
    
    # set out directory
    if args.outdir is None and slice_field in default_outdir.keys():
        outdir = default_outdir[slice_field]
    else:
        outdir = args.outdir

    return width, uselog, linthresh, zlo, zup, cmap, outdir

def get_mixing_region(ds):
    # first find if profile has been made
    profile_fname = f"profiles/{ds.basename}.csv"
    prof = pd.read_csv(profile_fname, index_col=0)
    
    na23_range = np.max(prof['X(Na23)']) - np.min(prof['X(Na23)'])
    ne23_range = np.max(prof['X(Ne23)']) - np.min(prof['X(Ne23)'])
    # estimate mixing edge as cutttoff for fraction of min/max
    cutoff_frac=0.9
    frac_na23 = np.max(prof['X(Na23)']) - cutoff_frac * na23_range
    frac_ne23 = np.max(prof['X(Ne23)']) - cutoff_frac * ne23_range
    
    # return nearest index to edge. ignore innermost region
    prof = prof[ prof['radius'] > 1e7 ]

    idx_min = np.argmin(np.abs(prof['X(Na23)'] - frac_na23))
    idx_max = np.argmin(np.abs(prof['X(Ne23)'] - frac_ne23))
    
    return (prof['radius'].iloc[idx_min], 'cm'), (prof['radius'].iloc[idx_max], 'cm')
        
if __name__=="__main__":
    try:
        ds = yt.load(args.infile, hint='maestro')
    except yt.utilities.exceptions.YTUnidentifiedDataType:
        ds = yt.load(args.infile)
    
    s, outdir = plot_slice(ds, args.field, args)
    s.save(outdir)
