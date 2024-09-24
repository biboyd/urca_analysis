import yt
from yt import derived_field
from yt.fields.api import ValidateParameter
import numpy as np
import unyt as u
from sys import argv

yt.enable_parallelism()

# define radvel direction fields
@derived_field(name=("gas", "urca_ratio"), sampling_type="local", 
               units="dimensionless", validators=[ValidateParameter(["center"])])
def _urca_ratio(field, data):
    return data[('boxlib', 'X(ne23)')]/data[('boxlib', 'X(na23)')]


# define radvel direction fields
@derived_field(name=("gas", "radvel_alongx"), sampling_type="local", 
               units="cm/s", validators=[ValidateParameter(["center"])])
def _radvel_alongx(field, data):
    c = data.get_field_parameter("center")[0]
    return data[('gas', 'radial_velocity')]*(data[('index', 'x')] - c)/data[('index', 'radius')]

@derived_field(name=("gas", "radvel_alongy"), sampling_type="local", 
               units="cm/s", validators=[ValidateParameter(["center"])])
def _radvel_alongy(field, data):
    c = data.get_field_parameter("center")[1]
    return data[('gas', 'radial_velocity')]*(data[('index', 'y')] - c)/data[('index', 'radius')]

@derived_field(name=("gas", "radvel_alongz"), sampling_type="local", 
               units="cm/s", validators=[ValidateParameter(["center"])])
def _radvel_alongz(field, data):
    c = data.get_field_parameter("center")[2]
    return data[('gas', 'radial_velocity')]*(data[('index', 'z')] - c)/data[('index', 'radius')]



def main(fname):
    # load dataset
    ds = yt.load(fname)

    # define convection zone size
    if len(argv) == 3:
        Rconv = float(argv[2]) * u.km
    else:
        Rconv = 600 * u.km

    conv_sph = ds.sphere(ds.domain_center, Rconv)

    # find dipole direction
    radvel_x = conv_sph.mean(('gas', 'radvel_alongx'), weight=('boxlib', 'rho'))
    radvel_y = conv_sph.mean(('gas', 'radvel_alongy'), weight=('boxlib', 'rho'))
    radvel_z = conv_sph.mean(('gas', 'radvel_alongz'), weight=('boxlib', 'rho'))

    # normalize dipole to unit vector
    dipole_direction = ds.arr([radvel_x, radvel_y, radvel_z])
    dipole = dipole_direction/np.linalg.norm(dipole_direction)

    # construct orthogonal vector
    ortho_to_dipole_direction = np.cross(dipole, np.array([1., 0., 0.]))
    ortho_to_dipole = ortho_to_dipole_direction / np.linalg.norm(ortho_to_dipole_direction)

    s = yt.SlicePlot(ds, ortho_to_dipole, 'vort', width=(1.2e3, 'km'), 
                     center=ds.domain_center, data_source=conv_sph, north_vector=dipole)
    
    s.set_zlim('all', 0., 18.)
    s.set_log('all', False)
    s.set_cmap('all', "viridis")
    s.save("vort_dipole/")

if __name__ == "__main__":
    main(argv[1])
