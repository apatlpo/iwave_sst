''' Utils for HW mask and sst processing
'''

import numpy as np
import xarray as xr
import dask.array as da
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
#from cmocean import cm
#import time
#import os
#from datetime import datetime
import ephem
import threading


r2d = 180./np.pi

#------------------------------ MASK ---------------------------------------

rint = xr.ufuncs.rint
fmod = xr.ufuncs.fmod

def process_raw_mask(mask):
    #QA:description = "
    #(2,1,0) Cloud Retrieval Algorithm Flag: 000=Outside of Scan, 001=No Cloud Mask,
    #        010=Clear, 011=Failed, 100=Successful: Low Confidence, 101=Successful: High Confidence, 110=TBD, 111=TBD; 
    #(4,3) Cloud Mask Confidence Level Flag: 00=Clear, 01=Probably Clear, 10=Probably Cloudy, 11=Cloudy; 
    #(6,5) Cloud Retrieval Phase Flag: 00=Clear, 01=Liquid Water, 10=Mixed or Uncertain, 11=Ice; 
    #(7) Spare: 0=TBD, 1=TBD;
    #(8) Sunglint Flag: 0=Yes, 1=No; 
    #(9) Snow Ice Background Possibility Flag: 0=Yes, 1=No; 
    #(11,10) Land/Water Flag: 00=Water, 01=Coastal, 10=TBD, 11=Land; 
    #(12) SOZ>80 or SAZ>70: 0=Yes, 1=No; 
    #(13) Subpixel Inhomogeneity Flag: 0=Yes, 1=No; 
    #(14) Multilayer Cloud Flag: 0=Yes, 1=No; 
    #(15) Inversion Layer Flag: 0=Yes, 1=No;" ;
    #
    fmask=xr.ones_like(mask)
    ### (2,1,0) Cloud Retrieval Algorithm Flag
    code=[]
    for i in range(3):
        code.append(rint(fmod(mask,2)))
        mask=mask//2
    fmask = fmask.where( (code[0]==0) & (code[1]==1) & (code[2]==0))
    ### (4,3) Cloud Mask Confidence Level Flag: 00=Clear, 01=Probably Clear, 10=Probably Cloudy, 11=Cloudy;
    code=[]
    for i in range(2):
        code.append(rint(fmod(mask,2)))
        mask=mask//2
    #fmask = fmask.where( (code[0]==0) & (code[1]==0) )
    #fmask = fmask.where( code[1]==1 )
    ### skips next
    for i in range(5):
        mask=mask//2
    ### (11,10) Land/Water Flag: 00=Water, 01=Coastal, 10=TBD, 11=Land; 
    code=[]
    for i in range(2):
        code.append(rint(fmod(mask,2)))
    #    #code.append(np.rint(mask%2))
        mask=mask//2
    fmask = fmask.where( (code[0]==0) & (code[1]==0) )
    #
    fmask = fmask.fillna(0.)
    return fmask

#
def plot_mask(mask, colorbar=False, title=None, vmin=0., vmax=1., savefig=None, offline=False, angles=None):
    MPL_LOCK = threading.Lock()
    with MPL_LOCK:
        if offline:
            plt.switch_backend('agg')
        fig = plt.figure(figsize=(10,10))
        #ax = plt.axes(projection=ccrs.Geostationary(central_longitude=140.0)) # may cause crash when distributed
        ax = fig.add_subplot(111, projection=ccrs.Geostationary(central_longitude=140.0))
        mask.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax,
                             x='longitude', y='latitude', add_colorbar=colorbar);
        #
        if angles is not None:
            im = (angles['angle2spec']*r2d).plot.contour(levels=[15.,30.,45.], colors=['orange'], 
                                                         ax=ax, transform=ccrs.PlateCarree(),
                                                         x='longitude', y='latitude', add_labels=True)
            im.clabel()
        #
        ax.coastlines(color='w')
        #
        if title is None:
            ax.set_title('HW cloud mask')
        else:
            ax.set_title(title)
        #
        if savefig is not None:
            fig.savefig(savefig, dpi=300)
            if offline:
                plt.close(fig)
        #
        if not offline:
            plt.show()
        #
        return fig, ax

#
def coarsen(fmask, dl, chunks=()):
    
    lon_bins = np.arange(fmask['longitude'].min().values, fmask['longitude'].max().values, dl)
    lat_bins = np.arange(fmask['latitude'].min().values, fmask['latitude'].max().values, dl)
    #
    lon_center = lon_bins[:-1]+dl
    lat_center = lat_bins[:-1]+dl

    # rename
    v1min, v1max, dv1= lon_bins[0], lon_bins[-1], dl
    v2min, v2max, dv2= lat_bins[0], lat_bins[-1], dl
    i1max = int(np.rint((v1max-v1min)/dv1))+1
    i2max = int(np.rint((v2max-v2min)/dv2))+1

    # meshgrid lon/lat, note: need transposing
    fmask = fmask.to_dataset()
    fmask['lon'] = (1.*fmask['longitude'] + 0.*fmask['latitude']).transpose()
    fmask['lat'] = (0.*fmask['longitude'] + 1.*fmask['latitude']).transpose()
    # need rechunking
    fmask = fmask.chunk(chunks)

    def get_index(v1,v2):
            ''' This function provides the index of (v1,v2) coupled value position
            in the 2D histogram array
            '''
            i1 = np.maximum(np.floor((v1-v1min)/dv1)+1,0)
            i1 = np.minimum(i1,i1max)
            i2 = np.maximum(np.floor((v2-v2min)/dv2)+1,0)
            i2 = np.minimum(i2,i2max)
            return i1+i2*(i1max+1)

    # sum QA over coarse grid cells
    v12 = da.map_blocks(get_index, fmask['lon'].data, fmask['lat'].data, dtype='float')
    h, lbins = da.histogram(v12, bins=np.arange(-.5,(i1max+1)*(i2max+1)+0.5,1.), weights=fmask['QA'].data)
    H = h.compute()
    # compute the number of points per grid cells
    hnorm, lbins = da.histogram(v12, bins=np.arange(-.5,(i1max+1)*(i2max+1)+0.5,1.))
    Hnorm = 1.*hnorm.compute()
    Hnorm[np.where(Hnorm==0)]=np.NaN
    # average the mask over coarse grid cells
    H = (H/Hnorm).reshape((i1max+1,i2max+1), order='F')

    cmask = xr.Dataset()
    #cmask['QA'] = (('longitude', 'latitude'), H[1:-1,1:-1].transpose())
    cmask['QA'] = (('longitude', 'latitude'), H[1:-1,1:-1])
    cmask.coords['longitude'] = (('longitude'),lon_center)
    cmask.coords['latitude'] = (('latitude'),lat_center)

    return cmask


#
def write_log(slog, clobber=False):
    if clobber:
        flog = open('hw_mask.log','w')
    else:
        flog = open('hw_mask.log','a')            
    flog.write(slog+'\n')
    flog.close()
    return

#
class twindow_manager():
    def __init__(self,threshold,Tmin,dl,t0):
        self.open = np.empty((0,4))
        self.closed = np.empty((0,4))
        # store other useful variables
        self.threshold = threshold
        self.Tmin = Tmin
        self.dl = dl
        self.t0 = t0
        #
        write_log('dl=%.2f deg, threshold=%.2f, Tmin=%.1f h'%(self.dl, self.threshold,self.Tmin*24.), clobber=True)
    def update(self,lon,lat,delt,time):
        if len(lon)>0:
            for llon, llat, ldelt in zip(lon,lat,delt):
                # test if open is empty
                if self.open.shape[0]>1:
                    # open has 2 elements at least
                    # test if window is already open
                    ij = np.where( (self.open[:,0]==llon) & (self.open[:,1]==llat) )
                    if len(ij[0])==0:
                        # items needs to be added to open list
                        self.open = np.concatenate((self.open,np.array([llon,llat,ldelt,time],ndmin=2)),axis=0)
                    else:
                        # items needs to be updated
                        self.open[ij[0],:] = np.array([llon,llat,ldelt,time])
                elif self.open.shape[0]==1:
                    # open has 1 element
                    if (self.open[0,0]==llon) & (self.open[0,1]==llat):
                        # items needs to be updated
                        self.open[0,:] = np.array([llon,llat,ldelt,time])
                    else:
                        # items needs to be added to open list
                        self.open = np.concatenate((self.open,np.array([llon,llat,ldelt,time],ndmin=2)),axis=0)
                else:
                    # open has 0 elements
                    self.open = np.array([llon,llat,ldelt,time],ndmin=2)
        # need to move inactive windows to closed
        if self.open.shape[0]>0:
            idel=[]
            for i in range(self.open[:,0].size):
                # search for matches in lon/lat
                ij = np.where( (lon==self.open[i,0]) & (lat==self.open[i,1]) )
                if len(ij[0])==0:
                    # the window needs to be closed
                    idel.append(i)
                    write_log('Stores: %.2f, %.2f, %.2f, %s ' \
                            %(self.open[i,0],self.open[i,1],self.open[i,2], \
                              str(t0+timedelta(seconds=self.open[i,3])) ) )
            self.open = np.delete(self.open,idel,axis=0)


# put everything in a function to see if it solves the memory issue
def process_mask_time(i, t, f, tagg, chunks, s):
    log = str(t)
    # load data
    mask = xr.open_dataset(f)['QA']
    # process
    fmask = process_raw_mask(mask)
    # coarsen
    cmask = coarsen(fmask, s.dl, chunks)
    # decimate
    mask = xr.ones_like(cmask['QA'])
    mask = mask.where(cmask['QA']>s.threshold) # keep only values above the threshold
    mask = mask.fillna(0.)
    #
    if i>0:
        delt = (t-s.tm1).total_seconds()/86400.
        if delt>0.25:
            # reset mask to 0 if time inverval between files exceeds 6h
            mask[:]=0.
        tagg += delt
        tagg *= mask
        tagg.compute()
        # store large values of tagg
        ij = np.where(tagg.values>=s.Tmin)
        s.update(tagg['longitude'].values[ij[0]], tagg['latitude'].values[ij[1]], \
                 tagg.values[ij], (t-s.t0).total_seconds())
    #
    #im1=i
    #s.im1 = i
    s.tm1 = t
    log += '  open: %d    closed: %d' %(s.open.shape[0],s.closed.shape[0])
    print('  open: %d    closed: %d' %(s.open.shape[0],s.closed.shape[0]))
    write_log(log)
    return tagg


#------------------------------ SST ---------------------------------------

#
def plot_sst(sst, colorbar=False, title=None, vmin=10., vmax=35., savefig=None, offline=False):
    if offline:
        plt.switch_backend('agg')
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    im = sst.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax,
                        x='lon', y='lat', add_colorbar=colorbar, cmap=cm.thermal)
    fig.colorbar(im)
    gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='k', alpha=0.5, linestyle='--')
    gl.xlabels_top=False
    ax.coastlines(color='k')
    #
    if title is None:
        ax.set_title('HW sst')
    else:
        ax.set_title(title)
    #
    if savefig is not None:
        fig.savefig(savefig, dpi=150)
        time.sleep(.1)
        plt.close(fig)
    #
    if not offline:
        plt.show()


        
#------------------------------ compute angle to specular reflection ---------------------------------------

# read himawari tle data
def read_hw_tle(time, tlepath):
    """ read Himawari8 TLE text files in order to build a pyephem body object
    Depending on the precision you want on the satellite position, this may 
    not be accurate enough
    
    TLE format: https://en.wikipedia.org/wiki/Two-line_element_set
    TLE databases:
        http://www.data.jma.go.jp/mscweb/en/operation8/orb_info/index.html
        http://celestrak.com/
    """
    file = tlepath+'tle_h8_%d.txt'%(time.year)
    f = open(file,'r')
    for line in f:
        if 'HIMAWARI-8' in line:
            tle=[line]
        elif line[0] == '1':
            tle.append(line)
            # store tle time
            lsplt = line.split()[3]
            tletime = ephem.Date('20'+lsplt[:2]+'/00')
            tletime+=float(lsplt[2:])
        elif line[0] == '2':
            tle.append(line)
            if np.abs(ephem.date(time)-tletime)<7.1:
                hw = ephem.readtle(tle[0], tle[1], tle[2])
                hw.compute(time)
                #print('hw: sublong=%f sublat=%f' % (hw.sublong*r2d, hw.sublat*r2d))
    f.close()
    return hw

##
cos = lambda theta: np.cos(theta)
sin = lambda theta: np.sin(theta)

def get_vector(body, az=None, alt=None):
    """ Get the local vector pointing toward a body given its position computed 
    given an observer (via pyephem) 
    The vector is in earth coordinates system: in e_r, e_east, e_north basis
    """
    if body is not None:
        return [sin(body.alt), cos(body.alt)*sin(body.az), cos(body.alt)*cos(body.az)]
    elif alt is not None and az is not None:
        return [sin(alt), cos(alt)*sin(az), cos(alt)*cos(az)]

def get_angle(v1, v2, acute=True):
    """ Compute the angle between two vectors
    """
    v1v2dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    nv1 = np.sqrt(v1[0]**2+v1[1]**2+v1[2]**2)
    nv2 = np.sqrt(v2[0]**2+v2[1]**2+v2[2]**2)
    angle = np.arccos( v1v2dot / ( nv1 * nv2 ) )
    if (acute == True):
        return angle
    else:
        return 2*np.pi - angle

def sun_zenith(sun, time):
    o = ephem.Observer()
    o.lon = 0.
    o.lat = 0.
    o.date = time
    o.elevation = 0.
    #
    sun.compute(o)
    v = get_vector(sun)
    sunlat = np.arcsin(v[2])
    sunlon = np.arctan2(v[1],v[0]) + o.lon
    return sunlon, sunlat

def spherical2cartesian(v, lon, lat):
    """ e_r, e_east, e_north to ex, ey, ez
    """
    M = [[cos(lat)*cos(lon), -sin(lon), -sin(lat)*cos(lon)],
         [cos(lat)*sin(lon) , cos(lon), -sin(lat)*sin(lon)],
         [sin(lat) , 0., cos(lat)]]
    vout = []
    for i in range(3):
        vout.append(lon*lat*0.) # enforces right coordinates
        for j in range(3):
            vout[i] = vout[i] + M[i][j]*v[j]
    return vout

def cartesian2spherical(v, lon, lat):
    """ ex, ey, ez to e_r, e_east, e_north
    """
    M = [[cos(lat)*cos(lon), -sin(lon), -sin(lat)*cos(lon)],
         [cos(lat)*sin(lon) , cos(lon), -sin(lat)*sin(lon)],
         [sin(lat) , 0., cos(lat)]]
    #for i in range(3): print(v[i]) # lon/lat, lon/lat, lat
    vout = []
    for i in range(3):
        vout.append(lon*lat*0.) # enforces right coordinates
        for j in range(3):
            vout[i] = vout[i] + M[j][i]*v[j]
    return vout

def get_azalt(lon, lat, r, olon, olat):
    """ compute azimuth and altitude from lon/lat/range
    """
    # compute vectors in a cartesian frame of reference
    v = [r*cos(lon)*cos(lat), r*sin(lon)*cos(lat), r*sin(lat)]
    #
    ro = ephem.earth_radius
    vo = [ro*cos(olon)*cos(olat), ro*sin(olon)*cos(olat), ro*sin(olat)]
    #
    for i in range(3): v[i] = v[i] - vo[i]
    #
    v = cartesian2spherical(v, olon, olat) # transform to spherical (er,eE,eN)
    nv = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    for i in range(3): v[i] = v[i]/nv
    az = np.arctan2(v[1],v[2])
    alt = np.arcsin(v[0])
    return az, alt

def angle2specular(sun_az, sun_alt, body_az, body_alt):
    """ Computes the angle to specular reflection
    """
    v_sun = get_vector(None, az=sun_az, alt=sun_alt)
    v1 = get_vector(None, az=body_az, alt=body_alt)    
    for i in range(3): v1[i] = v1[i] + v_sun[i]
    nv1 = np.sqrt(v1[0]**2+v1[1]**2+v1[2]**2)
    for i in range(3): v1[i] = v1[i]/nv1  
    # sun specular reflection:
    #v2 = [v_sun[0],-v_sun[1],-v_sun[2]]
    # normal
    v2 = [1.,0.,0.]
    return get_angle(v2,v1)

def get_reflection_angles(lon, lat, time, tlepath='./tle/'):
    """ Compute three relevant angles for Himawari sunglint analysis:
    angle to specular reflection, sun azimuth, sun altitude
    
    Parameters
    ----------
    lon, lat: xarray
        lon, lat in degrees
    time: datenum time
    """
    t = time.strftime('%Y/%m/%d %H:%M') #'2017/m/20 h:00'
    #
    sun = ephem.Sun()
    sun.compute(t)
    sunlon, sunlat = sun_zenith(sun, t)
    print('Sun zenith location: lon=%.1f, lat=%.1f' %(sunlon*r2d, sunlat*r2d))
    #
    hw = read_hw_tle(time, tlepath)
    print('Himawari position: sublong=%f sublat=%f' % (hw.sublong*r2d, hw.sublat*r2d))
    #
    ob = ephem.Observer()
    ob.lon = str(0.)
    ob.lat = str(0.)
    ob.elevation = 0.
    ob.date = t
    #
    sun.compute(ob)
    hw.compute(ob)
    # sun and himawari azimuth and altitude at observer positions
    sun_az, sun_alt = get_azalt(sunlon, sunlat, sun.earth_distance*ephem.meters_per_au, lon/r2d, lat/r2d)
    hw_az, hw_alt = get_azalt(hw.sublong, hw.sublat, hw.range, lon/r2d, lat/r2d)
    #
    angle2spec = angle2specular(sun_az, sun_alt, hw_az, hw_alt)
    #
    return xr.Dataset({'sun_az': sun_az, 'sun_alt': sun_alt,'angle2spec': angle2spec}), sun, hw
