#!/usr/bin/env python
import yt
import unyt
from yt.units import dimensions
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.api import Scene, create_volume_source, Camera, ColorTransferFunction
import numpy as np
import argparse

from matplotlib import colormaps
Blues = colormaps['Blues']
Reds = colormaps['Reds']

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-f', '--frac_isotope', type=str, default="na23", help='isotope field fraction to plot')
parser.add_argument('-rup', '--rup', type=float, default=1.0e8, help='Maximum radius (cm). Default 1.0e8.')
parser.add_argument('-zoom', '--zoom', type=float, default=1.0, help='Camera zoom factor. Default 1.0.')
parser.add_argument('-cpos', '--camera_position', type=float, nargs=3, default=[1., 0., 0.], help='3-D Camera position in fractions of maximum radius (--rup).')
parser.add_argument('-cnorth', '--camera_north', type=float, nargs=3, default=[0., 0., 1.], help='Camera north vector (direction of up).')
parser.add_argument('-rmin', '--frac_minimum', type=float, default=1e-7, help='Minimum frac for transfer function. (Default is 1e-7).')
parser.add_argument('-rmax', '--frac_maximum', type=float, default=1e-4, help='Maximum frac for transfer function. (Default is 1e-4).')
parser.add_argument('-sig', '--frac_sigma', type=float, default=0.01, help='frac transfer function width parameter. (Default is 0.01).')
parser.add_argument('-a', '--angle', type=float, default=0., help='angle in which to rotate about the north vector.')
parser.add_argument('-n', '--num_layers', type=int, default=4, help='Number of layers in transf func (Default is 4).')
parser.add_argument('-dd', '--drawdomain', action='store_true', help='If supplied, draw the boundaries of the domain.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-da', '--drawaxes', action='store_true', help='If supplied, draw an axes triad.')
parser.add_argument('-alpha_ones', '--alpha_ones', default=True, action='store_true', help='If supplied, set the transfer function values to ones.')
parser.add_argument('-res', '--resolution', type=int, default=2048, help='Resolution for output plot.')
parser.add_argument('-u', '--urca_rho', type=float, default=None, help='Plot density source (intended to be at urca shell 1.9e9 for myfix2048).')
parser.add_argument('-o', '--outprefix', type=str, default="", help='prefix to put at front of file')
parser.add_argument('-ptf', '--plot_tfunction', action='store_true',  help='plot the transferfunction files')
parser.add_argument('-dry', '--dry_run', action='store_true', help='Plot only the transfer functions and quit.')
args = parser.parse_args()

yt.enable_parallelism()


# Open Dataset
ds = yt.load(args.infile, hint='maestro')
field = ('boxlib', f"X({args.frac_isotope})")
core = ds.sphere(ds.domain_center, (args.rup, 'cm'))
# Create Scene
sc = Scene()

# Create Sources
so_frac = create_volume_source(core, field)
so_urca_shell = create_volume_source(core, 'rho')

# Assign Transfer Functions to Sources
# give buffer on bounds
bounds = np.array([args.frac_minimum/2, 2*args.frac_maximum])
log_min = np.log10(args.frac_minimum)
log_max = np.log10(args.frac_maximum)

if args.alpha_ones:
    alphavec = np.ones(args.num_layers)
else:
    alphavec = np.logspace(-3, 0, num=args.num_layers, endpoint=True)

tfh = TransferFunctionHelper(ds)
tfh.set_field(field)
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(bounds)
tfh.build_transfer_function()
tfh.tf.add_layers(args.num_layers,
               w= args.frac_sigma,
               mi=log_min,
               ma=log_max,
               alpha = alphavec,
               colormap = 'plasma')
#if args.plot_tfunction:
#    tfh.plot(f"{args.outprefix}{ds.basename}_tfun_ratio.png")#, profile_field=('boxlib', 'urca_ratio'))
             
so_frac.transfer_function = tfh.tf

if args.urca_rho is not None:
    tfh = TransferFunctionHelper(ds)
    tfh.set_field(('boxlib', 'rho'))
    tfh.set_log(True)
    tfh.grey_opacity = False
    tfh.set_bounds((1.e9, 4.5e9))
    tfh.build_transfer_function()
    tfh.tf.add_gaussian(np.log10(args.urca_rho), (2 * args.frac_sigma)**2, [1., 1., 1., 0.1]) # should give a white shell
    if args.plot_tfunction:
        tfh.plot(f"{args.outprefix}{ds.basename}_tfun_urca_shell.png") 
    so_urca_shell.transfer_function = tfh.tf

if args.dry_run:
    exit()

# Add sources to scene
sc.add_source(so_frac)
if args.urca_rho is not None:
    sc.add_source(so_urca_shell)

# Add camera to scene
sc.add_camera(ds, lens_type="perspective")

# Set camera properties
sc.camera.focus = ds.domain_center
sc.camera.resolution = args.resolution
sc.camera.north_vector = unyt.unyt_array(args.camera_north, 'cm')
sc.camera.position = ds.domain_center + unyt.unyt_array(args.camera_position, 'cm') * args.rup
sc.camera.set_width((2*args.rup, 'cm'))
# Annotate domain - draw boundaries
if args.drawdomain:
    sc.annotate_domain(ds, color=[1, 1, 1, 0.01])

# Annotate by drawing grids
if args.drawgrids:
    sc.annotate_grids(ds, alpha=0.01)

# Annotate by drawing axes triad
if True: #args.drawaxes:
    sc.annotate_axes(alpha=0.01) 

sc.camera.yaw(args.angle, rot_center=ds.domain_center)
sc.render()
sc.save(f"{args.outprefix}{ds.basename}_rendering_{args.frac_isotope}_frac.png", sigma_clip=4, render=False)
sc.save_annotated(f"{args.outprefix}{ds.basename}_annotated_rendering_{args.frac_isotope}_frac.png", sigma_clip=4,  render=False)
