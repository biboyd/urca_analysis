import yt
import pandas as pd
import matplotlib.pyplot
import numpy as np
from sys import argv

@yt.derived_field(name=("gas", "ratio_ad"), sampling_type='local')
def _ratio_ad(field, data):
  return data['boxlib', 'del']/data['boxlib', 'del_ad']

@yt.derived_field(name=("gas", "ratio_ledoux"), sampling_type='local')
def _ratio_led(field, data):
  return data['boxlib', 'del']/data['boxlib', 'del_ledoux']


def plot_conv_grad(df):
  fig, ax = plt.subplots(1, 1)

  ax.plot(df0['radius']/1e5, df0['ratio_ledoux'],  
        color='tab:blue', label='Ledoux')

  ax.plot(df0['radius']/1e5, df0['ratio_ad'], 
        color='tab:blue', alpha=0.7, ls='--', 
        lw=4, label='Adiabat')
  ax.set_xlim(0., 650.)
  ax.set_ylim(-2.75, 1.2)

  ax.set_xlabel("Radius (km)")
  ax.set_ylabel("$\\nabla / \\nabla_{\\mathrm{conv}}$")
  ax.legend(framealpha=1)
  ax.grid()
  ax.minorticks_on()
  #add surrounding ticks
  ax.tick_params(axis='both', direction='in', which='both', top=True, right=True)

  b,t = ax.get_ylim()
  ax.vlines(conv_size_arr[-2, 1], b,t, color='black', lw=2, alpha=0.5, zorder=1)
  ax.vlines(415, b, t, colors='black', zorder=-1, alpha=0.5, lw=2, linestyle='-.')

  ax.set_ylim(b, t)

  fig.tight_layout()
  fig.savefig("convective_gradient_ratio.png", dpi=300)
  
  return fig

def main(ds):
  
  prof_fields = [('boxlib', 'ratio_ad'), ('boxlib', 'ratio_ledoux'), ('boxlib', 'del_ad'), ('boxlib', 'del')]
  prof = yt.create_profile(ds.all_data(), 'radius', prof_fields,weight_field='volume', n_bins=200, logs={'radius':False}, extrema={'radius':[0., 1e8]} )

  df = prof.to_dataframe(include_std=True).dropna()

  # plot profiles
  fig = plot_conv_grad(df)

  #save profile
  df.to_csv(f"{ds.basename}_conv_grad.csv")


if __name__ == "__main__":
  ds = yt.load(argv[0])
  main(ds)
