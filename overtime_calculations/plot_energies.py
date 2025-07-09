import numpy as np
import matplotlib.pyplot as plt
import yt
import unyt
import pynucastro as pyna


def plot_energy(time_urca, stellar_urca, nuc_urca, nu_loss_urca,
                time_no, stellar_no, nuc_no, nu_loss_no):
    """
    inputs:
    ---------

    time_urca [arr]:    array holding the time for each snapshot in the Urca
                        simulation
    stellar_urca [arr]: sum of internal and gravitational energy for
                        each snapshot
    nuc_urca [arr]: The nuclear binding energy adjusted for mass loss
                    for each snaphshot
    nu_loss_urca [arr]: estimated total energy lost to neutrinos (erg)


    time_no [arr]:      array holding the time for each snapshot in the no
                        Urca simulation
                        simulation
    stellar_no [arr]: sum of internal and gravitational energy for
                        each snapshot
    nuc_no [arr]: The nuclear binding energy adjusted for mass loss
                    for each snaphshot
    nu_loss_no [arr]: estimated total energy lost to neutrinos (erg)

    outputs:
    ---------

    fig [plt.Figure] figure object of the plot
    ax [plt.axes] axes of the plot
    """

    fig, ax = plt.subplots(1, 1)

    # plot Urca portion
    l_urca, = ax.plot(time_urca, stellar_urca-stellar_urca[0], '-', label="Urca: Int+Grav")
    ax.plot(time_urca, np.abs(nuc_urca-nuc_urca[0]), '--', color=l_urca.get_color(), label="Urca: -Nuc")
    ax.plot(time_urca, nu_loss_urca, '.-', color=l_urca.get_color(), label="Urca: Nu loss")

    # plot No Urca portion
    l_no, = ax.plot(time_no, stellar_no-stellar_no[0], '-', label="No Urca: Int+Grav")
    ax.plot(time_no, np.abs(nuc_no-nuc_no[0]), '--', color=l_no.get_color(), label="No Urca: -Nuc")
    ax.plot(time_no, nu_loss_no, '.-', color=l_no.get_color(), label="No Urca: Nu loss")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Energy Difference (erg)")
    ax.grid()
    ax.legend()

    return fig, ax


def sum_energies(data_arr):
    """
    inputs:
    --------
    data_arr [arr]: output of load_data, includes time, energy and mass loss
                    data overtime for given simulation.

    outputs:
    ---------
    stellar_arr [arr]: sum of internal and gravitational energy for
                        each snapshot
    nuc_arr [arr]: The nuclear binding energy adjusted for mass loss
                    for each snaphshot
    nu_loss_arr [arr]: estimated total energy lost to neutrinos (erg)
    """

    time_arr, energy_arr, mass_loss_arr = data_arr

    # combine int and grav energies
    stellar_arr = energy_arr[1, :] + energy_arr[2, :]

    # adjust nuc for mass loss
    nuc_arr = energy_arr[0, :] - calc_nuc_mass_loss(mass_loss_arr)

    # grab nu losses
    nu_loss_arr = integrate_nu_loss(time_arr, energy_arr[3, :])
    # just add vals to match time array length
    nu_loss_arr = np.concat(([0.], nu_loss_arr, [nu_loss_arr[-1]]))

    return stellar_arr, nuc_arr, nu_loss_arr


def plot_energy_diff(time_urca, stellar_urca, nuc_urca, nu_loss_urca,
                     time_no, stellar_no, nuc_no, nu_loss_no):
    """
    inputs:
    ---------

    time_urca [arr]:    array holding the time for each snapshot in the Urca
                        simulation
    stellar_urca [arr]: sum of internal and gravitational energy for
                        each snapshot
    nuc_urca [arr]: The nuclear binding energy adjusted for mass loss
                    for each snaphshot
    nu_loss_urca [arr]: estimated total energy lost to neutrinos (erg)


    time_no [arr]:      array holding the time for each snapshot in the no
                        Urca simulation
                        simulation
    stellar_no [arr]: sum of internal and gravitational energy for
                        each snapshot
    nuc_no [arr]: The nuclear binding energy adjusted for mass loss
                    for each snaphshot
    nu_loss_no [arr]: estimated total energy lost to neutrinos (erg)

    outputs:
    ---------

    fig [plt.Figure] figure object of the plot
    ax [plt.axes] axes of the plot
    """


    fig, ax = plt.subplots(1, 1)

    # sum Energies
    energy_diff_urca = stellar_urca + nuc_urca + nu_loss_urca
    energy_diff_no = stellar_no + nuc_no + nu_loss_no

    # plot Urca portion
    ax.plot(time_urca, energy_diff_urca, '-', label="Urca")
    ax.plot(time_no, energy_diff_no, '-', label="No Urca")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Net Energy Difference (erg)")
    ax.grid()
    ax.legend()

    return fig, ax


def integrate_nu_loss(time_arr, nu_loss_arr):
    """
    inputs:
    --------

    time_arr [arr]: array holding the time for each snapshot
    nu_loss_arr [arr]: nuetrino loss rate for each snapshot (erg/s)

    outputs:
    ---------
    integrated_nu_loss [arr]: estimated total energy lost to neutrinos (erg)
    """

    # calc delta t midpoints (t_i+1 - t_i-1) / 2.
    delta_t_arr = (time_arr[2:] - time_arr[:-2])/2.

    # multiply nu loss rates by delta t
    integrated_nu_losses = np.cumsum(delta_t_arr * nu_loss_arr[1:-1])

    return integrated_nu_losses


def calc_nuc_mass_loss(mass_loss_arr, iso_name=None, composition=None):
    """
    inputs:
    ----------

    mass_loss_arr [arr] array of mass expanding out and being removed from our
                        summations
    iso_name [arr] array of mass expanding out and being lost

    outputs:
    ---------

    nuc_loss [arr]  the nuclear 'binding' energy of the mass which is lost off
                    the grid
    """

    # some nuclear data
    m_p = 1.672621777e-24*unyt.g  # g
    m_n = 1.674927351e-24*unyt.g  # g
    m_e = 9.10938291e-28*unyt.g  # g
    N_A = unyt.avogadros_number_mks.value / unyt.g  # nucleon/g

    if iso_name is None:

        iso_name = ['n', 'H1', 'He4', 'C12', 'C13', 'C14', 'N13', 'N14', 'O16',
                    'O17', 'O18', 'F18', 'F21', 'Ne20', 'Ne21', 'Ne22', 'Ne23',
                    'Na23', 'Na25', 'Mg24', 'Mg25']
    if composition is None:
        composition = [0., 0., 0., 0.4099082, 4e-05, 0., 0., 0., 0.576, 0., 0.,
                       0., 0., 0.000134, 3.74e-05, 0.0137, 0., 0.000142, 0.,
                       0., 3.84e-05]

    # calc specific bind per nucleon of given composition
    spec_bind = 0.
    for iso, comp in zip(iso_name, composition):
        nuc = pyna.Nucleus(iso)
        delta = (nuc.Z/nuc.A - 0.5)*(m_p + m_e - m_n)
        bind = nuc.nucbind*unyt.MeV - delta*unyt.c**2

        spec_bind += comp*bind

    # calc the bind energy of mass lost off the grid
    nuc_loss = np.zeros_like(mass_loss_arr)
    nuc_loss = mass_loss_arr*unyt.g * spec_bind * N_A

    return nuc_loss.in_cgs().value


def load_data(topdir='./'):
    """
    inputs:
    --------
    topdir  [string]: top directory that includes urca and no urca sim data


    outputs:
    ---------
    time_urca [arr]:    array holding the time for each snapshot in the Urca
                        simulation
    energy_urca [2d arr]:   array holding the respective nuclear, internal,
                            graviational energy, and nu loss rate for the Urca
                            simulation at the time corresponding to time_urca
    mass_loss_urca [arr]:   array holding the lost mass for each snapshot in
                            the Urca

    time_no [arr]:      array holding the time for each snapshot in the no
                        Urca simulation
                        simulation
    energy_no [2d arr]:     array holding the respective nuclear, internal,
                            graviational energy, and nu loss rate for the
                            no Urca simulation at the time corresponding to
                            time_no
    mass_loss_no [arr]: array holding the lost mass for each snapshot in the no
                        Urca simulation
    """

    # directory names rel to topdir
    urca_dir = 'bstate_on_urca_0.9Mconv_large/'
    no_urca_dir = 'no_urca_0.9Mconv_large'

    # load nuc data
    nuc_urca = np.load(f"{topdir}/{urca_dir}/nuc_energy_overtime.npy")
    nuc_urca = nuc_urca[:, np.argsort(nuc_urca[0, :])]
    nuc_urca = nuc_urca[:, 1:]  # skip first index as its t=0
    nuc_no = np.load(f"{topdir}/{no_urca_dir}/nuc_energy_overtime.npy")
    nuc_no = nuc_no[:, np.argsort(nuc_no[0, :])]

    # load internal energy data
    int_urca = np.load(f"{topdir}/{urca_dir}/int_energy_overtime.npy")
    int_urca = int_urca[:, np.argsort(int_urca[0, :])]
    int_urca = int_urca[:, 1:] # skip first index as its t=0
    int_no = np.load(f"{topdir}/{no_urca_dir}/int_energy_overtime.npy")
    int_no = int_no[:, np.argsort(int_no[0, :])]

    # load gravitational energy data
    grav_urca = np.load(f"{topdir}/{urca_dir}/grav_energy_overtime.npy")
    grav_urca = grav_urca[:, np.argsort(grav_urca[0, :])]
    grav_urca = grav_urca[:, 1:]  # skip first index as its t=0

    grav_no = np.load(f"{topdir}/{no_urca_dir}/grav_energy_overtime.npy")
    grav_no = grav_no[:, np.argsort(grav_no[0, :])]

    # load mass data
    mass_urca = np.load(f"{topdir}/{urca_dir}/mass_overtime.npy")
    mass_urca = mass_urca[:, np.argsort(mass_urca[0,:])]
    mass_urca = mass_urca[:, 1:]  # skip first index as its t=0

    mass_no = np.load(f"{topdir}/{no_urca_dir}/mass_overtime.npy")
    mass_no = mass_no[:, np.argsort(mass_no[0, :])]

    # load neutrino loss data
    nu_loss_urca = np.load(f"{topdir}/{urca_dir}/nu_loss_overtime.npy")
    nu_loss_urca = nu_loss_urca[np.argsort(nu_loss_urca[:, 0]), :]

    nu_loss_no = np.load(f"{topdir}/{no_urca_dir}/nu_loss_overtime.npy")
    nu_loss_no = nu_loss_no[np.argsort(nu_loss_no[:, 0]), :]

    # combine energy data
    energy_urca = np.vstack((nuc_urca[1, :], int_urca[1, :], grav_urca[1, :], nu_loss_urca[:, 1]))
    energy_no = np.vstack((nuc_no[1, :], int_no[1, :], grav_no[1, :], nu_loss_no[:, 1]))

    # calc mass loss
    mass_loss_urca = mass_urca[1, :] - mass_urca[1, 0]
    mass_loss_no = mass_no[1, :] - mass_no[1, 0]

    time_urca = mass_urca[0, :]
    time_no = mass_no[0, :]

    return time_urca, energy_urca, mass_loss_urca, time_no, energy_no, mass_loss_no


def main(data_dir="./", outdir='./'):
    # load in data
    overtime_data = load_data(topdir=data_dir)
    urca_data = overtime_data[:3]
    no_data = overtime_data[3:]

    energies_urca = sum_energies(urca_data)
    energies_no = sum_energies(no_data)

    # plot all energies
    fig_energy, ax = plot_energy(urca_data[0], *energies_urca, no_data[0], *energies_no)
    fig_energy.savefig(f"{outdir}/energies_overtime.png", dpi=300)

    # plot energy differences
    fig_diff, ax = plot_energy_diff(urca_data[0], *energies_urca, no_data[0], *energies_no)
    fig_diff.savefig(f"{outdir}/net_energy_diff_overtime.png", dpi=300)


if __name__ == "__main__":
    main()
