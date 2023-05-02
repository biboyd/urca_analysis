import yt
from sys import argv


def urca(field, data):

    return ((data["boxlib", "X(ne23)"] + data["boxlib", "X(na23)"])-0.0005)/0.0005 


def main(infile):
    ds = yt.load(infile, hint='amrex')
    ds.add_field(
        name=("boxlib", "urca"),
        function=urca,
        sampling_type="local")

    s = yt.SlicePlot(ds, 'x', ('boxlib', 'urca'), width=(1.2e8, 'cm'))

    s.set_zlim(('boxlib', 'urca'), -1e-4, 1e-4)
    s.set_log(('boxlib', 'urca'), True, linthresh=1e-8)
    s.set_cmap(('boxlib', 'urca'), 'PiYG')
    s.save("plots_urca_excess/")


if __name__ == '__main__':
    f = argv[1]

    main(f)
