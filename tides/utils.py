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

#------------------------------ gradients ---------------------------------------

def grad(d, lon='longitude', lat='latitude'):
    """ Compute the gradient of data on lon/lat grid
    """
    dx = diff_lonlat(d,lon)
    dy = diff_lonlat(d,lat)
    return dx, dy

def diff_lonlat(d, c):
    """ Compute the gradient of a variable laid out on a regular lon/lat grid
    """
    if c is 'longitude':
        # treats dateline data correctly
        di = lon_extension(d)
    else:
        di = d   
    #
    di = di.diff(c,label='lower')/di[c].diff(c,label='lower')/(111.e3)
    di = (di + di.roll(**{c:1}))*.5
    #
    if c is 'longitude':
        di = di/np.cos(np.pi/180.*di['latitude'])
    return di

def lon_extension(v):
    """ Extends data array longitudinally in order to include dateline
    points in a computation of gradients
    """
    v0 = v.sel(longitude=0.)
    v0['longitude']=360.
    return xr.concat([v,v0],dim='longitude')

#
def mom_inv(dvdx, dvdy, f, o, r=1./20./86400.):
    """ Inverse a linearized momentum equation
    """
    j = np.complex(0.,1.)
    _o = o + j*r
    uc = -g * (j*_o*dvdx - f*dvdy)/(_o**2-f**2)
    vc = -g * (f*dvdx + j*_o*dvdy)/(_o**2-f**2)
    return uc, vc

def ri2c(r,i):
    return r+np.complex(0,1.)*i

def c2ri(c):
    return np.real(c), np.imag(c)

#------------------------------ bathymetry ---------------------------------------

def load_bathy(lon=None, lat=None, b='etopo2'):
    ''' Load bathymetry
    '''
    #
    if b is 'etopo2':
        #bfile = '/home2/pharos/othr/aponte/bathy/ETOPO2v2c_f4.nc'
        bfile = './ETOPO2v2c_f4.nc'
        hb = -xr.open_dataset(bfile)['z']
    #
    if lon is not None and lat is not None:
        # should use xESMF
        from scipy.interpolate import RectBivariateSpline
        lonb = hb['x'].values
        latb = hb['y'].values
        #
        iroll = np.where(lonb<0)[0][-1]+1
        lonb = np.roll(lonb,-iroll)
        hb = np.roll(hb.values,-iroll,axis=1)
        lonb[lonb<0] = lonb[lonb<0] + 360.
        # should add a test for lon type
        hi = RectBivariateSpline(lonb, latb, hb.T, kx=1, ky=1)(lon,lat).T
        return xr.DataArray(hi, coords={'latitude': lat, 'longitude': lon}, dims=('latitude', 'longitude'))
    else:
        return hb
        
    
        