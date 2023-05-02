import yt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sys import argv

#file
plt_name=argv[1]
proj_name=argv[2]
yt.set_log_level(40)

# load data
ds = yt.load(plt_name, hint="maestro")
df_conv = pd.read_csv(f"{proj_name}", skiprows=2, header=None, delim_whitespace=True)

df_conv.columns = ['radius', 'actual', 'adiabatic', 'ledoux' ]
df_conv.head()

plt.plot(df_conv['radius']/1e5, df_conv['actual']-df_conv['adiabatic'], 'k', label="Schwarzchild")
plt.plot(df_conv['radius']/1e5, df_conv['actual']-df_conv['ledoux'], 'r--', label="Ledoux")
plt.ylim(-0.4, 0.2)
plt.legend()
left, right = plt.xlim(300, 600)
plt.hlines(0., left, right, color='k', linestyle='--', alpha=0.2)

plt.title(f"time: {ds.current_time:0.1f}")
plt.xlabel("Radius (km)")
plt.ylabel("Excess")
plt.savefig(f"plots_conv_projs/{ds.basename}_criteria.png")
