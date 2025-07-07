import yt
import unyt
import numpy as np
from sys import argv

def _mass(field, data):
    return data[("boxlib", "rho")]*data[("gas", "volume")] * unyt.g / unyt.cm**3

def _beta_rate(field, data):
    return data[("boxlib", "A23_beta_decay_rate")]/data[("boxlib", "X(Ne23)")] / unyt.s

def _ecap_rate(field, data):
    return data[("boxlib", "A23_electron_capture_rate")]/data[("boxlib", "X(Na23)")] / unyt.s

def _mass_beta(field, data):
    return data[("boxlib", "raw beta rate")] * data[("gas", "mass")]

def _mass_ecap(field, data):
    return data[("boxlib", "raw ecap rate")] * data[("gas", "mass")]

def _energy_rate(field, data):
    return data[('boxlib', 'Hnuc')]*unyt.erg/unyt.s/unyt.g * data[('gas', 'mass')]

def _tot_nu_loss(field, data):
    return  (data['thermal_nu_loss'] -1.0 * data[ 'A23_electron_capture_nu_loss']  -1.0 * data[ 'A23_beta_decay_nu_loss'] )* unyt.erg/unyt.g/unyt.s


def _tot_nu_energy_loss(field, data):
    return data[('boxlib', 'tot_nu_loss')]*data[('gas', 'mass')]

fname = argv[1]

# load plotfiles
ds = yt.load(f"{fname}")
ds_nu = yt.load(f"nu_loss.{fname}")

ds.add_field(
    name=("gas", "mass"),
    function=_mass,
    take_log=True,
    units='g',
    sampling_type="local")


ds.add_field(name=('boxlib', 'energy_rate'), 
            function=_energy_rate, 
            take_log=True, 
            units='erg/s', 
            sampling_type='local')

ds_nu.add_field(
    name=("gas", "mass"),
    function=_mass,
    take_log=True,
    units='g',
    sampling_type="local")


ds_nu.add_field(
    name=("boxlib", "raw beta rate"),
    function=_beta_rate,
    take_log=True,
    units='1/s',
    sampling_type="local")

ds_nu.add_field(
    name=("boxlib", "raw ecap rate"),
    function=_ecap_rate,
    take_log=True,
    units='1/s',
    sampling_type="local")

ds_nu.add_field(
    name=("boxlib", "mass beta"),
    function=_mass_beta,
    take_log=True,
    units='g/s',
    sampling_type="local")

ds_nu.add_field(
    name=("boxlib", "mass ecap"),
    function=_mass_ecap,
    take_log=True,
    units='g/s',
    sampling_type="local")

ds_nu.add_field(
    name=("boxlib", "tot_nu_loss"),
    function=_tot_nu_loss,
    take_log=True,
    units='erg/g/s',
    display_name="Energy Loss to Neutrino Emissions",
    sampling_type="local")


ds_nu.add_field(
    name=("boxlib", "tot_nu_energy_loss"),
    function=_tot_nu_energy_loss,
    take_log=True,
    units='erg/s',
    display_name="Total Power Loss to Neutrino Emissions",
    sampling_type="local")

tot_nu_loss = ds_nu.all_data().sum("tot_nu_energy_loss")
tot_energy_gen = ds.all_data().sum("energy_rate")

print("Total energy generated", tot_energy_gen)
print("Total energy lost to neutrinos", tot_nu_loss)
