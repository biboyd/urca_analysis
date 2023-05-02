import yt
import matplotlib.pyplot as plt
from sys import argv

fname = argv[1]
if (fname[-1] == '/'):
    fname = fname[:-1]

ds = yt.load(fname, hint='amrex')

# get frbs
frbs = yt.SlicePlot(ds, 'z', [('boxlib', 'conv_actual'), ('boxlib', 'conv_adiabatic'), ('boxlib', 'conv_ledoux')], width=(1.2e3, 'km')).frb
act = frbs.data[('boxlib', 'conv_actual')]
ad = frbs.data[('boxlib', 'conv_adiabatic')]
led = frbs.data[('boxlib', 'conv_ledoux')]

s = act > ad
l = act > led

conv_crit = (s.astype(int) + l.astype(int))/2.

# plot some stuff
fig, ax = plt.subplots(1, 2, figsize=(10, 5 ))
ax[0].imshow(s.astype(int), cmap='RdBu')
ax[0].set_title("Schwarzchild Criterion")

ax[1].imshow(l.astype(int), cmap='RdBu')
ax[1].set_title("Ledoux Criterion")

fig.suptitle('Blue --> "Convectively Unstable"')
#fig.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes([0.85, 0.1, 0.075, 0.8])
#plt.colorbar(cax=cax)

fig.savefig(f"{fname}_convect.png")
