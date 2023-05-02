import yt
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

def doit(ds):

    # a FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this.
    
    left_edge =  ds.arr([1.96e8, 1.96e8, 1.96e8], 'cm')
    right_edge =  ds.arr([3.16e8, 3.16e8, 3.16e8], 'cm')

    max_level = ds.index.max_level

    low = left_edge #ds.domain_left_edge
    dims = np.array([1024, 1024, 1024], dtype=np.int64)

    nx, ny, nz = dims

    nindex_rho = 1.0 / 3.0

    Kk = np.zeros((nx // 2 + 1, ny // 2 + 1, nz // 2 + 1))

    for vel in [("boxlib", "velx"), ("boxlib", "vely"), ("boxlib", "velz")]:

        Kk += 0.5 * fft_comp(
            ds, ("boxlib", "rho"), vel, nindex_rho, max_level, low, dims
        )

    # wavenumbers
    L = (right_edge - left_edge).d

    kx = np.fft.rfftfreq(nx) * nx / L[0]
    ky = np.fft.rfftfreq(ny) * ny / L[1]
    kz = np.fft.rfftfreq(nz) * nz / L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0 / L)
    kmax = np.min(0.5 * dims / L)

    kbins = np.arange(kmin, kmax, kmin)
    N = len(kbins)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

    whichbin = np.digitize(k.flat, kbins)
    ncount = np.bincount(whichbin)

    E_spectrum = np.zeros(len(ncount) - 1)

    for n in range(1, len(ncount)):
        E_spectrum[n - 1] = np.sum(Kk.flat[whichbin == n])

    k = 0.5 * (kbins[0 : N - 1] + kbins[1:N])
    E_spectrum = E_spectrum[1:N]

    index = np.argmax(E_spectrum)
    kmax = k[index]
    Emax = E_spectrum[index]

    plt.loglog(k, E_spectrum)
    plt.loglog(k, Emax * (k / kmax) ** (-5.0 / 3.0), ls=":", color="0.5")

    plt.xlabel(r"$k (\mathrm{cm^{-1}})$")
    plt.ylabel(r"$\mathrm{E}(k)\mathrm{dk}$")

    plt.savefig(f"{ds.basename}_spectrum.png")
    return E_spectrum, k


def fft_comp(ds, irho, iu, nindex_rho, level, low, delta):

    cube = ds.covering_grid(level, left_edge=low, dims=delta, fields=[irho, iu])

    rho = cube[irho].d
    u = cube[iu].d

    nx, ny, nz = rho.shape

    # do the FFTs -- note that since our data is real, there will be
    # too much information here.  fftn puts the positive freq terms in
    # the first half of the axes -- that's what we keep.  Our
    # normalization has an '8' to account for this clipping to one
    # octant.
    ru = np.fft.fftn(rho**nindex_rho * u)[
        0 : nx // 2 + 1, 0 : ny // 2 + 1, 0 : nz // 2 + 1
    ]
    ru = 8.0 * ru / (nx * ny * nz)

    return np.abs(ru) ** 2

if __name__ == '__main__':
    dsname = argv[1]

    ds = yt.load(dsname, hint='amrex')

    E_arr, k_arr = doit(ds)

    np.save(f"{ds.basename}_Espec.npy", E_arr)
    np.save(f"{ds.basename}_k.npy", k_arr)
