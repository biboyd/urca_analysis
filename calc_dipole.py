import yt
import numpy as np
from sys import argv

#import dataset
ds = yt.load(argv[1])

cx = ds.domain_center[0]
cy = ds.domain_center[1]
cz = ds.domain_center[2]

def _radvel_alongx(field, data):
    return data[('boxlib', 'radial_velocity')]*(data[('index', 'x')] - cx)/data[('index', 'radius')]
def _radvel_alongy(field, data):
    return data[('boxlib', 'radial_velocity')]*(data[('index', 'y')] - cy)/data[('index', 'radius')]
def _radvel_alongz(field, data):
    return data[('boxlib', 'radial_velocity')]*(data[('index', 'z')] - cz)/data[('index', 'radius')]

ds.add_field(name=("gas", "radvel_alongx"),
            function=_radvel_alongx,
            units = "cm/s",
            take_log=False,
            sampling_type="local")

ds.add_field(name=("gas", "radvel_alongy"),
            function=_radvel_alongy,
            units = "cm/s",
            take_log=False,
            sampling_type="local")

ds.add_field(name=("gas", "radvel_alongz"),
            function=_radvel_alongz,
            units = "cm/s",
            take_log=False,
            sampling_type="local")

# grab convection zone
if len(argv) == 3:
    conv_sph = ds.sphere(ds.domain_center, (float(argv[2]), 'km'))
else:
    conv_sph = ds.sphere(ds.domain_center, (600, 'km'))

# density weighted means
density_radvel_base = conv_sph.mean(('boxlib', 'radial_velocity'), weight=('boxlib', 'rho'))
density_radvel_x = conv_sph.mean(('gas', 'radvel_alongx'), weight=('boxlib', 'rho'))
density_radvel_y = conv_sph.mean(('gas', 'radvel_alongy'), weight=('boxlib', 'rho'))
density_radvel_z = conv_sph.mean(('gas', 'radvel_alongz'), weight=('boxlib', 'rho'))
density_phi = np.arctan(density_radvel_y/density_radvel_x)
density_theta = np.arctan(np.sqrt(density_radvel_x**2 + density_radvel_y**2)/density_radvel_z)

print("Density weighted radvel")
print('radvel', density_radvel_base.in_units('km/s'))
print('radvel along x', density_radvel_x.in_units('km/s'))
print('radvel along y', density_radvel_y.in_units('km/s'))
print('radvel along z', density_radvel_z.in_units('km/s'))
print('angle phi', density_phi)
print('angle theta', density_theta)


# need to adjust values since arctan only output -pi/2 to pi/2
# move from IV to II and I to III.
if density_radvel_x < 0 and density_radvel_y > 0:
    density_phi += np.pi
if density_radvel_x < 0 and density_radvel_y < 0:
    density_phi += np.pi
if density_radvel_x > 0 and density_radvel_y < 0:
    density_phi += 2*np.pi

#if z is negative we need to adjust so 0 = horizontal and pi= down z
if density_radvel_z < 0:
    density_theta+= np.pi

# Unweighted
radvel_base = conv_sph.mean(('boxlib', 'radial_velocity'))
radvel_x = conv_sph.mean(('gas', 'radvel_alongx'))
radvel_y = conv_sph.mean(('gas', 'radvel_alongy'))
radvel_z = conv_sph.mean(('gas', 'radvel_alongz'))
phi = np.arctan(radvel_y/radvel_x)
theta = np.arctan(np.sqrt(radvel_x**2 + radvel_y**2)/radvel_z)

# need to adjust values since arctan only output -pi/2 to pi/2
# move from IV to II and I to III.
if radvel_x < 0 and radvel_y > 0:
    phi += np.pi
if radvel_x < 0 and radvel_y < 0:
    phi += np.pi
if radvel_x > 0 and radvel_y < 0:
    phi += 2*np.pi

#if z is negative we need to adjust so 0 = horizontal and pi= down z
if radvel_z < 0:
    theta+= np.pi

print()
print()
print("Unweighted Mean")
print('radvel', radvel_base.in_units('km/s'))
print('radvel along x', radvel_x.in_units('km/s'))
print('radvel along y', radvel_y.in_units('km/s'))
print('radvel along z', radvel_z.in_units('km/s'))
print('angle phi', phi)
print('angle theta', theta)
