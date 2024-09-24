import yt
from sys import argv


def comp(field, data):

    return ((data["boxlib", "X(ne23)"] + data["boxlib", "X(na23)"] + data["boxlib", "X(c12)"] + data["boxlib", "X(o16)"])-1)


def main(infile):
    ds = yt.load(infile, hint='amrex')
    ds.add_field(
        name=("boxlib", "comp"),
        function=comp,
        sampling_type="local")

    s = yt.SlicePlot(ds, 'x', ('boxlib', 'comp'), width=(1.2e8, 'cm'))

    s.set_zlim(('boxlib', 'comp'), -1e-4, 1e-4)
    s.set_log(('boxlib', 'comp'), True, linthresh=1e-8)
    s.set_cmap(('boxlib', 'comp'), 'PiYG')
    s.save("plots_urca_excess/")


if __name__ == '__main__':
    f = argv[1]

    main(f)
