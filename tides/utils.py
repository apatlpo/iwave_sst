''' Utils for tidal processing
'''

import numpy as np
import xarray as xr
import dask.array as da
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
#from cmocean import cm
import time

g=9.81
cpd=2.*np.pi/86400. # radian/s

#------------------------------ plots ---------------------------------------

def grad(d,lon=None,lat=None,mask=None):
    """ Compute the gradient of data on lon/lat grid
    Should take mask as input
    """
    if lon is None:
        lon = 'longitude'
    if lat is None:
        lat = 'latitude'
    dx = ndiff(d,lon)/np.cos(np.pi/180.*d[lat])
    dy = ndiff(d,lat)
    return dx, dy

def ndiff(d,c):
    #c0 = d[c].values
    di = d.diff(c,label='lower')/d[c].diff(c,label='lower')/(111.e3)
    di = (di + di.roll(**{c:1}))*.5
    return di

def mom_inv(dvdx,dvdy,f,o):
    i = np.complex(0.,1.)
    uc = -g * (i*o*dvdx-f*dvdy)/(o**2-f**2)
    vc = -g * (f*dvdx+i*o*dvdy)/(o**2-f**2)
    return uc, vc

def ri2c(r,i):
    return r+np.complex(0,1.)*i

def c2ri(c):
    return np.real(c), np.imag(c)
