import numpy as np
import yt
import pynucastro as pyna
import unyt
from os import listdir

m_p = 1.672621777e-24 #g
m_n = 1.674927351e-24 #g
m_e = 9.10938291e-28 #g

ds = yt.load("plotfiles/plt0008400/")

comp_list=[]
bion_list=[]
for fld_type, fld in ds.field_list:
    if 'X(' == fld[:2]:
        comp_list.append(fld[2:-1])
        nuc = pyna.Nucleus(fld[2:-1])
        delta = (nuc.Z/nuc.A - 0.5)*(m_p + m_e - m_n) *unyt.g 
        bind = nuc.nucbind*unyt.MeV - delta*unyt.c**2
        bion_list.append(bind)

@yt.derived_field(name=('boxlib', 'specific_nuc_energy'), units='erg/g', sampling_type='local')
def _spec_nuc_energy(field, data):

    nuc_tot = 0.
    for c, b in zip(comp_list, bion_list):
        nuc_tot += data[('boxlib', f"X({c})")] * b
    return -nuc_tot * unyt.avogadros_number_mks.value / unyt.g

@yt.derived_field(name=('boxlib', 'nuc_energy'), units='erg', sampling_type='local')
def _spec_nuc_energy(field, data):
    return data[('boxlib', 'specific_nuc_energy')] * data[('gas' ,'cell_mass')]
    

rconv=7e7

t_list = []
nuc_list = []
mass_list = []

for f in np.sort(listdir("plotfiles/")):
    if '000' in f:
        ds = yt.load(f"plotfiles/{f}")
        sph = ds.sphere(ds.domain_center, (rconv, 'cm'))
        tot_nuc_energy = sph.sum('nuc_energy')
        nuc_list.append(tot_nuc_energy.value)    
        mass_list.append(sph.sum('mass').value)    
        t_list.append(ds.current_time)

save_arr = np.array((t_list, nuc_list, mass_list))
np.save("nuc_conv_energy_overtime.npy", save_arr)
