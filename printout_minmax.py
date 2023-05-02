import yt
from sys import argv
from os import listdir

indir=argv[1]
plt_list = sorted(listdir(indir))

yt.enable_parallelism()
ds_list=[]
for p in plt_list:
    ds_list.append(yt.load(f"{indir}/{p}", hint='amrex'))

for curr_ds in ds_list:
    fld=('boxlib', 'tfromp')
    
    fld_max=curr_ds.find_max(fld)[0].value
    fld_min=0#curr_ds.find_min(fld)[0].value
    
    if yt.is_root():
        print(curr_ds.filename)
        print(f"time: {curr_ds.current_time.value: 0.2f} s.")
        print(f"max: {fld_max: 0.2e} cm/s.")
        print(f"min: {fld_min: 0.2e} cm/s.\n\n")
