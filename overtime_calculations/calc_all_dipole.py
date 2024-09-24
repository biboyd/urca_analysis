import yt
import numpy as np
from sys import argv
from os import listdir
import pandas as pd

yt.set_log_level(40)


def _radvel_alongx(field, data):
    return data[('gas', 'radial_velocity')]*(data[('index', 'x')] - cx)/data[('index', 'radius')]
def _radvel_alongy(field, data):
    return data[('gas', 'radial_velocity')]*(data[('index', 'y')] - cy)/data[('index', 'radius')]
def _radvel_alongz(field, data):
    return data[('gas', 'radial_velocity')]*(data[('index', 'z')] - cz)/data[('index', 'radius')]

all_files = listdir("plotfiles")
output = np.zeros((len(all_files), 8), dtype=np.float64) -1.

#previous data
try:
    old_dipole_info = np.load("dipole_info_overtime.npy")
    calc_all = False
    
except FileNotFoundError:
    calc_all = True
    print("Could not find old files, will be calculating for each file")
    
# grab convect zone array
try:
    conv_zone_arr = np.load("conv_zone_over_time.npy")
except FileNotFoundError:
    print("Could not find the conv zone file. Make sure to run convect_overtim.py first.")
    raise FileNotFoundError("conv_zone_over_time.npy")

for j, f in enumerate(all_files):
    if f[:3] == "plt" and len(f) == 10:
        #import dataset
        ds = yt.load(f"plotfiles/{f}")
        print(ds.basename)
        curr_time = ds.current_time 
        timestep = float(ds.basename.removeprefix("plt"))
        cx = ds.domain_center[0]
        cy = ds.domain_center[1]
        cz = ds.domain_center[2]
        
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
        
        #see if calced before
        if not calc_all and timestep in old_dipole_info[:, -1]:
            idx = np.argmin(np.abs(timestep -  old_dipole_info[:, -1]))
            output[j, :] = old_dipole_info[idx, :]
            print('already have saved')
            continue
            
        # calc convect zone
        # grab convection zone
        
        if len(argv) == 2:
            # take user input
            rad = (float(argv[1]), 'km')
            
        else:
            #read from file
            try:
                assert timestep in conv_zone_arr[:, 2], f"{timestep} not in conv zone array"
                idx = np.argmin(np.abs(timestep - conv_zone_arr[:, 2]))
                rad = ds.quan(conv_zone_arr[idx, 1], 'km')
            except AssertionError:
                print('no profile found for:', ds.basename)
                #fall back to a default
                rad = (600 , 'km')
                
        # make sure sphere is a real size lol.
        try:
            conv_sph = ds.sphere(ds.domain_center, rad)
        except yt.utilities.exceptions.YTSphereTooSmall:
            print('sphere too small ', rad, 'for file ', ds.basename)
            continue
        
        weighted=True
        # density weighted means
        if weighted:
            radvel_base = conv_sph.mean(('gas', 'radial_velocity'), weight=('boxlib', 'rho'))
            radvel_x = conv_sph.mean(('gas', 'radvel_alongx'), weight=('boxlib', 'rho'))
            radvel_y = conv_sph.mean(('gas', 'radvel_alongy'), weight=('boxlib', 'rho'))
            radvel_z = conv_sph.mean(('gas', 'radvel_alongz'), weight=('boxlib', 'rho'))
            phi = np.arctan(radvel_y/radvel_x)
            theta = np.arctan(np.sqrt(radvel_x**2 + radvel_y**2)/radvel_z)
            print("Density weighted Mean")
            
        else:
            # Unweighted
            radvel_base = conv_sph.mean(('gas', 'radial_velocity'))
            radvel_x = conv_sph.mean(('gas', 'radvel_alongx'))
            radvel_y = conv_sph.mean(('gas', 'radvel_alongy'))
            radvel_z = conv_sph.mean(('gas', 'radvel_alongz'))
            phi = np.arctan(radvel_y/radvel_x)
            theta = np.arctan(np.sqrt(radvel_x**2 + radvel_y**2)/radvel_z)
            print("Unweighted Mean")
            
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
            
        #print()
        #print()
        #print('radvel', radvel_base.in_units('km/s'))
        #print('radvel along x', radvel_x.in_units('km/s'))
        #print('radvel along y', radvel_y.in_units('km/s'))
        #print('radvel along z', radvel_z.in_units('km/s'))
        #print('angle phi', phi)
        #print('angle theta', theta)

        output[j, :] = [curr_time, radvel_base, radvel_x, radvel_y, radvel_z, phi, theta, timestep]


#sort array by sim time
sorted_output = output[ np.argsort(output[:, 0]) ]
sorted_output = sorted_output[ np.where(sorted_output[:, 0] > 0) ]
np.save("dipole_info_overtime.npy", sorted_output)
