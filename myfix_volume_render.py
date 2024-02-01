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
parser.add_argument('-rup', '--rup', type=float, default=1.0e8, help='Maximum radius (cm). Default 1.0e8.')
parser.add_argument('-zoom', '--zoom', type=float, default=1.0, help='Camera zoom factor. Default 1.0.')
parser.add_argument('-cpos', '--camera_position', type=float, nargs=3, help='3-D Camera position in fractions of maximum radius (--rup).')
parser.add_argument('-cnorth', '--camera_north', type=float, nargs=3, help='Camera north vector (direction of up).')
parser.add_argument('-vmin', '--velocity_minimum', type=float, default=1.0e-2, help='Minimum velocity for transfer function. (Default is 1.0e-2 km/s).')
parser.add_argument('-vmax', '--velocity_maximum', type=float, default=1.0e2, help='Maximum velocity for transfer function. (Default is 1.0e2 km/s).')
parser.add_argument('-v', '--velocity_center', type=float, default=2., help='center velocity for transfer function in log10. (Default is 2 (1e2 km/s)).')
parser.add_argument('-vsig', '--velocity_sigma', type=float, default=0.05, help='Velocity transfer function width parameter. (Default is 0.05).')
parser.add_argument('-a', '--angle', type=float, default=0., help='angle in which to rotate about the north vector.')
parser.add_argument('-n', '--num_layers', type=int, default=5, help='Number of layers for each of +/- velocity. (Default is 5).')
parser.add_argument('-dd', '--drawdomain', action='store_true', help='If supplied, draw the boundaries of the domain.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-da', '--drawaxes', action='store_true', help='If supplied, draw an axes triad.')
parser.add_argument('-alpha_ones', '--alpha_ones', action='store_true', help='If supplied, set the transfer function values to ones.')
parser.add_argument('-res', '--resolution', type=int, default=2048, help='Resolution for output plot.')
parser.add_argument('-u', '--urca_rho', type=float, default=None, help='Plot density source (intended to be at urca shell 1.9e9 for myfix2048).')
parser.add_argument('-o', '--outprefix', type=str, default="", help='prefix to put at front of file')
parser.add_argument('-ptf', '--plot_tfunction', action='store_true',  help='plot the transferfunction files')
parser.add_argument('-dry', '--dry_run', action='store_true', help='Plot only the transfer functions and quit.')
args = parser.parse_args()

yt.enable_parallelism()

# Hack: because rendering likes log fields ...
## create positive_radial_velocity and negative_radial_velocity fields.
## must do this before opening dataset
def _pos_radial_velocity(field, data):
    return np.maximum(data[('boxlib','radial_velocity')], yt.YTQuantity(1.0e-99, 'km/s'))
def _neg_radial_velocity(field, data):
    return np.maximum(-data[('boxlib','radial_velocity')], yt.YTQuantity(1.0e-99, 'km/s'))

# Open Dataset
ds = yt.load(args.infile, hint='maestro')

ds.add_field(name=("boxlib", "pos_radial_velocity"),
            function=_pos_radial_velocity,
            take_log=True,
            display_name="Outflow Velocity",
            units = "km/s",
            sampling_type="local")
    
ds.add_field(name=("boxlib", "neg_radial_velocity"),
            function=_neg_radial_velocity,
            take_log=True,
            display_name="Inflow Velocity",
            units = "km/s",
            sampling_type="local")

core = ds.sphere(ds.domain_center, (args.rup, 'cm'))
# Create Scene
sc = Scene()

# Create Sources
so_pos_vrad = create_volume_source(core, 'pos_radial_velocity')
so_neg_vrad = create_volume_source(core, 'neg_radial_velocity')
so_urca = create_volume_source(core, 'rho')

mag_vel_bounds = np.array([np.min([args.velocity_minimum, args.velocity_center])/2., np.max([args.velocity_maximum, args.velocity_center])*2.])
mag_vel_sigma  = args.velocity_sigma
log_min = np.log10(args.velocity_minimum)
log_max = np.log10(args.velocity_maximum)

nlayers = args.num_layers
if args.alpha_ones:
    alphavec = np.ones(nlayers)
else:
    alphavec = np.logspace(-3, 0, num=nlayers, endpoint=True)

# positive velocity
tfh = TransferFunctionHelper(ds)
tfh.set_field('pos_radial_velocity')
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(mag_vel_bounds)
tfh.build_transfer_function()
tfh.tf.add_gaussian(np.log10(args.velocity_center), mag_vel_sigma**2, Reds(1.))
tfh.tf.add_gaussian(np.log10(args.velocity_center/2.), mag_vel_sigma**2, Reds(0.5, alpha=0.1))
tfh.tf.add_gaussian(np.log10(args.velocity_center/4.), mag_vel_sigma**2, Reds(0.25, alpha=0.01))

if args.plot_tfunction:
    tfh.plot(f"{args.outprefix}{ds.basename}_tfun_pos_vrad.png")#, profile_field=('boxlib', 'pos_radial_velocity'))
so_pos_vrad.transfer_function = tfh.tf

# negative velocity
tfh = TransferFunctionHelper(ds)
tfh.set_field('neg_radial_velocity')
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(mag_vel_bounds)
tfh.build_transfer_function()
tfh.tf.add_gaussian(np.log10(args.velocity_center), mag_vel_sigma**2, Blues(1.))
tfh.tf.add_gaussian(np.log10(args.velocity_center/2.), mag_vel_sigma**2, Blues(.5, alpha=0.1))
tfh.tf.add_gaussian(np.log10(args.velocity_center/4.), mag_vel_sigma**2, Blues(.25, alpha=0.01))

if args.plot_tfunction:
    tfh.plot(f"{args.outprefix}{ds.basename}_tfun_neg_vrad.png")#, profile_field=('neg_radial_velocity'))
so_neg_vrad.transfer_function = tfh.tf

if args.urca_rho is not None:
    tfh = TransferFunctionHelper(ds)
    tfh.set_field(('boxlib', 'rho'))
    tfh.set_log(True)
    tfh.grey_opacity = False
    tfh.set_bounds((1.e9, 4.5e9))
    tfh.build_transfer_function()
    tfh.tf.add_gaussian(np.log10(args.urca_rho), (mag_vel_sigma/10.)**2, [1., 1., 1., 0.1]) # should give a white shell
    if args.plot_tfunction:
        tfh.plot(f"{args.outprefix}{ds.basename}_tfun_urca_shell.png") 
    so_urca.transfer_function = tfh.tf

if args.dry_run:
    exit()

# Add sources to scene
sc.add_source(so_pos_vrad)
sc.add_source(so_neg_vrad)
if args.urca_rho is not None:
    sc.add_source(so_urca)

# Add camera to scene
sc.add_camera(ds, lens_type="perspective")

# Set camera properties
sc.camera.focus = ds.domain_center
sc.camera.resolution = args.resolution
sc.camera.north_vector = unyt.unyt_array(args.camera_north, 'cm')
sc.camera.position = ds.domain_center + unyt.unyt_array(args.camera_position, 'cm') * args.rup
sc.camera.set_width((3*args.rup, 'cm'))
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
sc.save(f"{args.outprefix}{ds.basename}_rendering_rad-vel.png", sigma_clip=4, render=False)
sc.save_annotated(f"{args.outprefix}{ds.basename}_annotated_rendering_rad-vel.png", sigma_clip=4,  render=False, label_fmt="%.2d")
