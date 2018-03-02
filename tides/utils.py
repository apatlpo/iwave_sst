''' Utils for tidal processing
'''

import numpy as np
import xarray as xr
import dask.array as da
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
#from cmocean import cm
import time
from utide._ut_constants import ut_constants as utide

g=9.81
cpd=2.*np.pi/86400. # radian/s

#------------------------------ HRET reader -------------------------------------

def get_hret_ssh(constituents=['M2','N2','S2','K1','O1','P1'], lonb=None, latb=None, 
            hret='./Carrere_HRET_testing.nc', bathy=False):
    ''' Load HRET ssh
    
    Parameters
    ----------
    consituents: list of str
        Consistuents considered
    lonb: None, list, tuple
        Bounds for longitude, e.g.: lonb=(10.,100.)
    latb: None, list, tuple
        Bounds for latitude, e.g.: latb=(-10.,10.)
    hret: str
        Path to HRET netcdf file
    '''
    hret = xr.open_dataset(hret, chunks={'longitude': 500, 'latitude': 500})
    hret_constituents = ['M2','N2','S2','K1','O1','P1']
    for c in hret_constituents:
        if c not in constituents:
            del hret[c+'re']
            del hret[c+'im']        
    #
    omega = dict()
    for cst,o in zip(utide['const']['name'], utide['const']['freq']):
        if cst in constituents:
            omega[cst] = o*24. # cpd, input frequencies are cph
            print(cst+' omega=%e rad/s, %.3f cpd'%(o*2.*np.pi/3600., o*24.))
    #
    if lonb is not None:
        # should handle 360 wrapping
        hret = hret.where(hret['longitude']>=lonb[0], drop=True)
        hret = hret.where(hret['longitude']<=lonb[1], drop=True)
    if latb is not None:
        # should check conventions for longitude
        hret = hret.where(hret['latitude']>=latb[0], drop=True)
        hret = hret.where(hret['latitude']<=latb[1], drop=True)
    #
    if bathy:
        h = load_bathy(lon=hret.longitude, lat=hret.latitude)
        hret = hret.assign(h=h)
    #
    return hret, constituents, omega


def get_hret_uv(**kwargs):
    ''' Load HRET currents
    
    Parameters
    ----------
    consituents: list of str
        Consistuents considered
    lonb: None, list, tuple
        Bounds for longitude, e.g.: lonb=(10.,100.)
    latb: None, list, tuple
        Bounds for latitude, e.g.: latb=(-10.,10.)
    hret: str
        Path to HRET netcdf file
    '''
    hret, constituents, omega = get_hret_ssh(**kwargs)
    #
    omega_K1 = 15.04107*np.pi/180./3600. # deg/h -> rad/s
    print(omega_K1)
    f = 2*omega_K1*cpd*np.sin(np.pi/180.*hret['latitude'])
    #
    U = xr.Dataset()
    V = xr.Dataset()
    for cst in constituents:
        eta = ri2c(hret[cst+'re'], hret[cst+'im'])
        detadx_c, detady_c = grad(eta)
        U[cst], V[cst] = mom_inv(detadx_c, detady_c, f, omega[cst]*cpd)
    return U, V, constituents, omega

#------------------------------ gradients ---------------------------------------

#
def grad(d, lon='longitude', lat='latitude'):
    """ Compute the gradient of data on lon/lat grid
    """
    dx = di_lonlat(d,lon)
    dy = di_lonlat(d,lat)
    return dx, dy

def di_lonlat(d, c):
    """ Compute the gradient of a variable laid out on a regular lon/lat grid
    """
    if c is 'longitude':
        # treats dateline data correctly
        di = lon_extension(d)
    else:
        di = d   
    #
    dx = di[c].diff(c,label='lower')*111.e3
    di = di.diff(c,label='lower')/dx
    di = (di + di.shift(**{c:1}))*.5
    #
    if c is 'longitude':
        di = di/np.cos(np.pi/180.*di['latitude'])
    return di

#
def lap(d, lon='longitude', lat='latitude'):
    """ Compute the laplacian of data on lon/lat grid
    """
    d2x = di2_lonlat(d,lon)
    d2y = di2_lonlat(d,lat)
    d2x, d2y = xr.align(d2x, d2y, join='outer')
    lap = d2x + d2y
    return lap

def di2_lonlat(d, c):
    """ Compute the second derivative of a variable laid out on a regular lon/lat grid
    """
    if c is 'longitude':
        # treats dateline data correctly
        di = lon_extension(d)
    else:
        di = d   
    #
    dx = di[c].diff(c,label='lower')*111.e3    
    di = di.diff(c,label='lower')/dx
    di = (di - di.shift(**{c:1}))/dx
    #
    if c is 'longitude':
        di = di/np.cos(np.pi/180.*di['latitude'])**2
    return di

#
def lon_extension(v):
    """ Extends data array longitudinally in order to include dateline
    points in a computation of gradients
    """
    if v['longitude'].min()==0. and v['longitude'].max()==359.95:
        v0 = v.sel(longitude=0.)
        v0['longitude']=360.
        return xr.concat([v,v0],dim='longitude')
    else:
        return v
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
        
    
        