#!/usr/bin/python3
"""
Plotting script to illustrate convergence of geometric properties with theta.

To use this script, edit the `local.in` input file so that ntheta is 6, 8, 16, 32, 64,
or 128, and run the main gs2 executable on that file. Copy the resulting `local.out.nc`
file to `new_<ntheta>.out.nc`. Once all these data files exist, run this script:
    $ ./plot_geometry.py
The resulting figure, geometry_vs_theta.png, shows that geometric properties converge
with increasing resolution. 

This test has not been automated, as the theta grids produced are non-uniform, so that
comparisons would require interpolation.

This script was written to test updates to the geometry calculations after the release
of 8.0.2. To see that the new approach is an improvement, repeat the above procedure
using GS2 <= v8.0.2, copy the data files to `old_<ntheta>.out.nc`, and uncomment the
`make_plot` lines below. The resulting plot shows that:
    1. the old and new methods converge to the same functions in the high ntheta limit
    2. the new method gives a better approximation to the converged function for any
       fixed ntheta.
    3. the new method with ntheta points and the old with with ntheta + 2 points give
       the same quality of approximation.

This script can be used to test any of the geometries in the directories `asym_geo_*`.
"""

from netCDF4 import Dataset
from matplotlib.pylab import *
import numpy as np

def make_plot(iname,col,mark,lw=0):
    """
    Read data files and plot geometry properties.
        iname: data file name without .out.nc suffix
        col: marker/line colour
        mark: marker style
        lw: linewidth, default 0 (no line)
    """

    fname = iname+".out.nc"
    d = Dataset(fname)
    theta = d.variables['theta']
    bmag = d.variables['bmag']
    bpol = d.variables['bpol']
    gbdrift = d.variables['gbdrift']
    jacob = d.variables['jacob']

    ax[0,0].plot(theta,bmag,marker=mark,color=col,lw=lw)
    ax[0,1].plot(theta,bpol,marker=mark,color=col,lw=lw)
    ax[1,0].plot(theta,gbdrift,marker=mark,color=col,lw=lw)
    ax[1,1].plot(theta,jacob,marker=mark,color=col,lw=lw)

fig, ax = plt.subplots(ncols=2,nrows=2, sharex=False, sharey=False)

make_plot("new_128",'b','',lw=0.3)
make_plot("new_64",'b','+')
make_plot("new_32",'b','x')
make_plot("new_16",'b','.')
make_plot("new_8",'c','.')
make_plot("new_6",'c','+')

# Uncomment the following lines to plot "old" data
#make_plot("old_64",'r','+')
#make_plot("old_32",'r','x')
#make_plot("old_16",'r','.')
#make_plot("old_8",'m','.')
#make_plot("old_6",'m','+')

# Set axis labels
ax[0,0].set_xlabel("theta")
ax[0,1].set_xlabel("theta")
ax[1,0].set_xlabel("theta")
ax[1,1].set_xlabel("theta")
ax[0,0].set_ylabel("bmag")
ax[0,1].set_ylabel("bpol")
ax[1,0].set_ylabel("gbdrift")
ax[1,1].set_ylabel("jacob")

plt.savefig("geometry_vs_theta.png",dpi=600)
