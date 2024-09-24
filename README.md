# URCA Analysis scripts
these are pretty basic plotting scripts.

The `make_slice.py` and `make_slice_w_grid.py` scripts are very basic and should just be used as a reference for how to make simple plots

The `plot_figures.py` is meant to be more useable with some flexibility built in. The default values are meant to reflect the `myfix_urca1024` and `myfix_urca2048` sims. But this is the most general purpose script for making slices. Run `python plot_figures.py --help` for details on what parameters you can edit. This plot will take in a `plotfile` and a `field` to plot. Then create a `SlicePlot` with the set parameters chosen (ie cmap, width etc.) and save a `.png` file.

## 3D Volume Rendering
These scripts under the `volume_rendering` directory make 3D volume renderings of various fields for the Urca problem. Most notably is the `radvel_volume_render.py` which produces a render of the Radial Velocity, showing a clear description of the velocity structure.

These scripts are quite computational expensive so most have a corresponding slurm script in `perlmutter_scripts/` or `sewulf_scripts/` 

## Calculating Quantities Over Time
These scripts take in a directory full of `plotfiles` and calculates some broad property (Radius of the Convection zone, Average entropy etc.) saves an array of the calculated quantity in some `.npy` file. This can then be quickly loaded using `numpy.load("file_overtime.npy")`

The most importatn one here is the `convect_overtime.py`. This calculates the radius of the convection zone and then a series of parameters in said convection zone:

* $R_{\mathrm{conv}}$: Radius of convection zone in km (needs 1D profiles of snapshots see `save_profiles.py`)
* $M_{\mathrm{conv}}$: Mass in the convection zone in solar mass 
* $V_{\mathrm{rms}}$:  The rms velocity in the convection zone in km/s
* $<V>$: The average velocity in the convection zone in km/s. Note: should use $V_{\mathrm{rms}}$ in most cases
* $\tau_{\mathrm{conv}}$:  A measure of the convective turnover time in seconds. Simply $2 R_{\mathrm{conv}} / V_{\mathrm{rms}}$

This information (especially the $R_{\mathrm{conv}}$) is used in many other analyses.

## Single snapshot calculations
These scripts are somewhat similar to the `overtime` scripts, but they only calculate a quantity on one plotfile. These are generally more complex quanttites ie the kinetic energy power spectrum, and so that is why they are more isolated scripts.

## 1D Profiles
In `save_plot_profiles` we there are scripts for both calculating 1D profiles (this is done by averaging values into radial bins) and plotting the profiles. This is often done to see how profiles vary with time.


## 2D Slices
These scripts in `slice_plotting` take a `plotfile` and makes a slice plot of a certain quantity/field. 

The `dipole_*.py` files aligns the plot in the direction of the velocity dipole. So the general flow of these plots should be "up". 

## Inputs Files to C scripts
There are a few analysis scripts in C that access the data more directly or use libraries that we need to compile and aren't available via python.

The sample inputs files for these scripts are placed in `c_scripts_inputs`

These relate to calculating the convective criteria, calculating the seperate neutrino losses, and shrinking plotfiles to remove unneeded/redundant fields.

## Others
Other scripts which don't fall into these categories are here. These are more helpful in cerrtain situations and are not doing any hard analysis.