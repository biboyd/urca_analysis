"""
Want to make a quick script that will make slice plots and save them so I can look at them goodlike
"""

import numpy as np
import yt
from sys import argv

def plot_slice(ds, field):
    slice = yt.SlicePlot(ds, 'x', field, width=12e7)

    slice.save()


if __name__ == '__main__':
    plotfile = argv[1]

    field_str = argv[2]

   
    ds = yt.load(plotfile, hint='amrex')

    plot_slice(ds, ('boxlib', field_str))
