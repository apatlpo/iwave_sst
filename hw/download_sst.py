
# coding: utf-8

# 
# # Play with JAXA Himawari 
# 
# In order to search for clear skies
# https://worldview.earthdata.nasa.gov/
# 
# ```
# ### /jma/
# ## Himawari Standard Data
# 
# ### /pub/
# ## Himawari Geophysical Parameter Data
# 
# ### /pub/himawari : lower resolution (2km), full disk
# #
# #Sea Surface Temperature (SST) Level 2 (near-real-time):
# #/pub/himawari
# #       +---/L2
# #             +---/SST
# #                   +---/[VER]_nc4_normal_nrt
# #                          +---/[YYYYMM]
# #                                 +---/[DD]
# 
# ### Sea Surface Temperature
#  YYYYMMDDhhmmss-JAXA-L2P_GHRSST-SSTskin-H08_AHI-vVER-v02.0-fvFVER.nc
# 
#  where YYYY: 4-digit year of observation start time (timeline);
#        MM: 2-digit month of timeline;
#        DD: 2-digit day of timeline;
#        hh: 2-digit hour of timeline;
#        mm: 2-gidit minutes of timeline;
#        ss: 2-digit seconds (fixed to "00");
#        VER: algorithm version; and
#        FVER: file version;.
# 
#  Example: 
#  20150728081000-JAXA-L2P_GHRSST-SSTskin-H08_AHI-v1.0-v02.0-fv01.0.nc
#  
# ```
# 


# Import libraries
import os, sys
import urllib
#from urllib.error import URLError, HTTPError
#from urllib import error
#import urllib2
import ftplib as ftplib
import numpy as np
import matplotlib.pyplot as plt
from  matplotlib.dates import date2num, datetime, num2date
from netCDF4 import Dataset


#dpath = '/Users/aponte/Current_projects/anr_jcjc/work/data/'
dpath = '/home/slyne/aponte/sst_fast/data_nwa_06/'; bbox_name = 'NWA_06'
#dpath = '/home/slyne/aponte/sst_fast/data_nwa_08/'; bbox_name = 'NWA_08'
#dpath = '/home/slyne/aponte/sst_fast/data_nwa_09/'; bbox_name = 'NWA_09'

# stdbuf -oL python himawari.py >& hw_nwa_06.log
# stdbuf -oL python himawari.py >& hw_nwa_08.log
# stdbuf -oL python himawari.py >& hw_nwa_09.log

# reference time for netcdf file
t0 = date2num(datetime.datetime(2015,7,1,0,0,0))

#
# download data
#
def readwrite_hw(time, lon, lat, bbox=None, first=False):

    # format date and time
    YYYYMM = time.strftime('%Y%m')
    DD = time.strftime('%d')
    YYYYMMDDhhmm = time.strftime('%Y%m%d%H%M')

    # SST data, nc format
    #VER='v100'; VERd='v1.0'
    VER='v101'; VERd='v1.1'
    FVER='01.0'

    # filename
    fdir='/pub/himawari/L2/SST/'+VER+'_nc4_normal_nrt/'+YYYYMM+'/'+DD+'/'
    fname=YYYYMMDDhhmm+'00-JAXA-L2P_GHRSST-SSTskin-H08_AHI_NRT-'+VERd+'-v02.0-fv'+FVER+'.nc'
    fname_bbox = YYYYMMDDhhmm+'_'+bbox_name+'.nc'

    # open ftp
    host, login, passwd = ftp_login()
    ftp = ftplib.FTP(host)
    ftp.login(login, passwd)

    if not os.path.isfile(dpath+fname_bbox):
        print 'No data file found, download from jaxa ftp: '+fname
        file = open(dpath+fname,'wb')
        #
        ftp = ftplib.FTP(host)
        ftp.login(login, passwd)
        ftp.cwd(fdir)
        try:
            ftp.retrbinary('RETR %s' % fname, file.write)
            file.close()
            ftp.quit()
        except:
            print "Error"
        #attempts=0
        #while attempts < 3:
        #   try:
        #        print 'test0'
        #        #urllib2.Request(ftp+fdir+fname, dpath+fname)
        #        urllib.urlretrieve(ftp+fdir+fname, dpath+fname)
        #        break
        #    except error.URLError as e:
        #    #except
        #        print 'test1'
        #        attempts += 1
        #        print ' - Failed attempt '+str(attempts)
        #        continue
        #    except HTTPError as e:
        #        print 'test2'
    else:
        print fname_bbox+' is already on disk'

    # very first file or not
    if first:
        fname_bbox = fname
    else:
        # create bbox file
        nc_create(dpath+fname_bbox, lon, lat)
        # load sst data
        if os.path.getsize(dpath+fname)>1:
            nc = Dataset(dpath+fname, 'r')
            sst = nc.variables['sea_surface_temperature'][0, bbox[2]:bbox[3], bbox[0]:bbox[1]]
            nc.close()
        else:
            flog = open(dpath+fname+'.log','wb')
            flog.write('Issue with this file')
            flog.close()
            sst = None
        # append
        nc_append(dpath+fname_bbox, sst, time)
        # delete full data file
        if sst is not None:
            os.remove(dpath+fname)

    return dpath+fname_bbox


def ftp_login():
    ''' Read ftp address, login and password
    '''
    f = open('ftp.login','r')
    host = f.readline().strip('\n')
    login = f.readline().strip('\n')
    passwd = f.readline().strip('\n')
    f.close()
    return host, login, passwd


#
# create netcdf file where data is appended
#
def nc_create(fname, lon,lat):

    nc = Dataset(fname, 'w', format='NETCDF4_CLASSIC', clobber=True)

    # create dimensions
    nc.createDimension('lon', lon.size)
    nc.createDimension('lat', lat.size)
    nc.createDimension('t', None)

    # create variables
    dtype = 'f8'
    nc_x = nc.createVariable('lon', dtype, ('lon'))
    nc_y = nc.createVariable('lat', dtype, ('lat'))
    nc_x[:] = lon
    nc_y[:] = lat
    #
    nc.createVariable('t', dtype, ('t'))
    #
    nc.createVariable('sst', dtype, ('t', 'lat', 'lon',))

    nc.close()

    return

#
# append to netcdf file where data is append
#
def nc_append(fname, sst, time):

    # open netcdf file
    nc = Dataset(fname, 'a', format='NETCDF4_CLASSIC')

    it = nc.variables['sst'].shape[0]
    nc.variables['t'][it] = date2num(time)-t0
    if sst is None:
        nc.variables['sst'][it,...] = 0.
    else:
        nc.variables['sst'][it,...] = sst

    nc.close()

    return


#
# Main part of the script
#
if __name__ == "__main__":

    # bounding box
    lonb = [105., 125.]
    latb = [-23., -10]

    # time line
    t1=datetime.datetime(2015,6,15,0,0,0); delt = 30. # in days
    # first file is in fact: 201507/07/20150707015000-JAX...
    #t1=datetime.datetime(2015,7,26,0,0,0); delt = 30. # in days
    #t1=datetime.datetime(2015,9,1,0,0,0); delt = 30.

    time = [t1+datetime.timedelta(minutes=10*n) for n in xrange(int(delt*24*60/10))]

    # download a first file in order to get coordinates
    fname = readwrite_hw(time[0], None, None, first=True)
    nc = Dataset(fname,'r')
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    nc.close()

    # find bounding box limits
    ilon1 = np.where(lon>lonb[0])[0][0]
    ilon2 = np.where( (lon<lonb[1]) & (lon>0.) )[0][-1]
    ilat1 = np.where(lat<latb[1])[0][0]
    ilat2 = np.where(lat>latb[0])[0][-1]
    bbox = [ilon1, ilon2, ilat1, ilat2]

    # download other files and store in one single file
    fnames = [readwrite_hw(t, lon[ilon1:ilon2], lat[ilat1:ilat2], bbox=bbox) for t in time]

#     for filen, t in zip(fname,time):
# 
#         nc = Dataset(filen, 'r')
#         sst = nc.variables['sea_surface_temperature'][0, ...]
#         nc.close()
# 
#         plt.figure(figsize=(10,10))
#         plt.pcolormesh(sst)
#         plt.savefig('figs/sst_'+t.strftime('%Y%m%d%H%M')+'.png',dpi=300)


