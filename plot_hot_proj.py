import yt
from sys import argv

ds = yt.load(argv[1], hint='amrex')
# sphere
sph = ds.sphere(ds.domain_center, (6e7, 'cm'))

hot = ds.cut_region(sph, "obj[('boxlib', 'tfromp')] > 6.e8")

proj = yt.AxisAlignedProjectionPlot(ds, normal=0, fields='tfromp', width=(1.2e8, 'cm'), 
                                     weight_field=('boxlib', 'rho'), data_source=hot )
proj.set_cmap('all', 'magma')
proj.annotate_timestamp()
proj.set_background_color('all')
proj.set_zlim('all', 7e7, 7e8)
    
proj.save()
