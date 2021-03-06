#!/usr/bin/env python

# generate figures in Getting Started section of User's Manual

# usage:
#   $ python basemapfigs.py FILEROOT [FIELD] [DPI]
# where
#   FILEROOT   root of NetCDF filename and output .png figures
#   FIELD      optional: one of {cbase, [csurf], mask, usurf}  (edit script to add more)
#   DPI        optional: resolution in dots per inch [200]
#
# equivalent usages:
#   $ python basemapfigs.py g20km_10ka_hy csurf 200
#   $ python basemapfigs.py g20km_10ka_hy csurf
#   $ python basemapfigs.py g20km_10ka_hy
#
# generate figs like those in Getting Started section of User's Manual:
#   $ for FLD in csurf usurf cbase mask; do ./basemapfigs.py g20km_10ka_hy ${FLD}; done
#
# crop out western Greenland with command like this (uses ImageMagick):
#   $ ./basemapfigs.py g20km_10ka_hy csurf 500
#   $ convert -crop 600x800+400+800 +repage g20km_10ka_hy-csurf.png g20km-detail.png
#
# batch generate figures from a parameter study like this:
#   $ for QQ in 0.1 0.5 1.0; do for EE in 1 3 6; do ../basemapfigs.py p10km_q${QQ}_e${EE} csurf 100; done; done
#   $ for QQ in 0.1 0.5 1.0; do for EE in 1 3 6; do convert -crop 274x486+50+6 +repage p10km_q${QQ}_e${EE}-csurf.png p10km-${QQ}-${EE}-csurf.png; done; done

from mpl_toolkits.basemap import Basemap

try:
    from netCDF4 import Dataset as NC
except:
    print "netCDF4 is not installed!"
    sys.exit(1)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import sys

if len(sys.argv) < 2:
    print "ERROR: first argument must be root of filename ..."
    sys.exit(1)
rootname = sys.argv[1]
try:
    nc = NC(rootname+'.nc', 'r')
except:
    print "ERROR: can't read from file %s.nc ..." % rootname
    sys.exit(2)

if len(sys.argv) >= 3:
  field = sys.argv[2]
else:
  field = 'csurf'

if len(sys.argv) >= 4:
  mydpi = float(sys.argv[3])
else:
  mydpi = 200

bluemarble = False  # FIXME: m.basemap() inverts earth!

if (field == 'csurf') | (field == 'cbase'):
  fill       = nc.variables[field]._FillValue
  logscale   = True
  contour100 = True
  myvmin     = 1.0
  myvmax     = 6.0e3
  ticklist   = [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
elif field == 'surfvelmag':
  fill       = 0.0
  logscale   = True
  contour100 = True
  myvmin     = 1.0
  myvmax     = 6.0e3
  ticklist   = [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
elif field == 'usurf':
  fill       = 0.0
  logscale   = False
  contour100 = False
  myvmin     = 1.0
  myvmax     = 3500.0
  ticklist   = [100, 500, 1000, 1500, 2000, 2500, 3000, 3500]
elif field == 'mask':
  fill       = -1.0
  logscale   = False
  contour100 = False
  myvmin     = 0.0
  myvmax     = 4.0
  ticklist   = [0, 1, 2, 3, 4]
else:
  print 'invalid choice for FIELD option'
  sys.exit(3)

# we need to know longitudes and latitudes corresponding to grid
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
if field == 'surfvelmag':
    lon = np.squeeze(lon).transpose()
    lat = np.squeeze(lat).transpose()

# x and y *in the dataset* are only used to determine plotting domain
# dimensions
if field == 'surfvelmag':
    x = nc.variables['x1'][:]
    y = nc.variables['y1'][:]
else:
    x = nc.variables['x'][:]
    y = nc.variables['y'][:]
width  = x.max() - x.min()
height = y.max() - y.min()

# load data and mask out ice-free areas
myvar = np.squeeze(nc.variables[field][:])
if field == 'surfvelmag':
    myvar = myvar.transpose()
myvar = np.ma.array(myvar, mask=(myvar == fill))

m = Basemap(width=1.1*width,    # width in projection coordinates, in meters
            height=1.05*height,      # height
            resolution='l',     # coastline resolution, can be 'l' (low), 'h'
                                # (high) and 'f' (full)
            projection='stere', # stereographic projection
            lat_ts=71,          # latitude of true scale
            lon_0=-41,          # longitude of the plotting domain center
            lat_0=72)           # latitude of the plotting domain center

#m.drawcoastlines()

# draw the Blue Marble background (requires PIL, the Python Imaging Library)
if bluemarble:  # seems to reverse N and S
    m.bluemarble()

# convert longitudes and latitudes to x and y:
xx,yy = m(lon, lat)

if contour100:
  # mark 100 m/a contour in black:
  m.contour(xx, yy, myvar, [100], colors="black")

# plot log color scale or not
if logscale:
  m.pcolormesh(xx,yy,myvar,
               norm=colors.LogNorm(vmin=myvmin, vmax=myvmax))
else:
  m.pcolormesh(xx,yy,myvar,vmin=myvmin,vmax=myvmax)

# add a colorbar:
plt.colorbar(extend='both',
             ticks=ticklist,
             format="%d")

# draw parallels and meridians
# labels kwarg is where to draw ticks: [left, right, top, bottom]
m.drawparallels(np.arange(-55.,90.,5.), labels = [1, 0, 0, 0])
m.drawmeridians(np.arange(-120.,30.,10.), labels = [0, 0, 0, 1])

outname = rootname+'-'+field+'.png'
print "saving image to file %s ..." % outname
plt.savefig(outname, dpi=mydpi, bbox_inches='tight')
