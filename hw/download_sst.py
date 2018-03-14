# coding: utf-8
# 
# Download JAXA Himawari SST data
#
# log:
#     python download_sst.py -c 120.,127.,-17,-10 -n NWAM -s 2017,4,1,0,0,0 -d 60 
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
from optparse import OptionParser
import ftplib as ftplib
import numpy as np
import datetime
#from netCDF4 import Dataset
import xarray as xr



#
# Parse inputs
#
def parse_input():
    """ Parse inputs
    """
    parser = OptionParser()
    parser.add_option("-c", "--coordinates", dest="c",
                      help="coordinate bounds, ex: -c lonmin,lonmax,latmin,latmax)", \
                      default='105.,125.,-40.,-10')
    parser.add_option("-n", "--name", dest="name",
                      help="box name", default='NWA')    
    parser.add_option("-s", "--timestart", dest="tstart",
                      help="starting time, ex: -s 2015,9,15,0,0,0", default='2015,9,15,0,0,0')
    parser.add_option("-d", "--delt", dest="delt",
                      help="time interval in days", default='10.')
    (options,args) = parser.parse_args()
    #
    c = [float(v) for v in options.c.split(',')]
    box = {'lon': slice(c[0],c[1]), 'lat': slice(np.max(c[2:]),np.min(c[2:])), 'name': options.name}
    #
    tstart = [int(v) for v in options.tstart.split(',')]
    tstart=datetime.datetime(*tstart)
    # 
    return box, tstart, float(options.delt)

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
# download data
#
def download_hw(times, rpath, dpath=None, box=None):
        
    # open ftp
    host, login, passwd = ftp_login()
    ftp = ftplib.FTP(host)
    ftp.login(login, passwd)
    print('Logged in JAXA ftp')

    # loop around times
    sstfiles=[]
    for time in times:
        
        # format date and time
        YYYYMM = time.strftime('%Y%m')
        DD = time.strftime('%d')
        YYYYMMDDhhmm = time.strftime('%Y%m%d%H%M')
        #
        #VER='v100'; VERd='v1.0'
        #VER='v101'; VERd='v1.1'
        VER = 'v102'; VERd = 'v1.2'
        FVER = '01.0'
        #
        #fdir='/pub/himawari/L2/SST/'+VER+'_nc4_normal_nrt/'+YYYYMM+'/'+DD+'/'
        fdir = '/pub/himawari/L2/SST/'+VER+'_nc4_normal_std/'+YYYYMM+'/'+DD+'/'
        #fname_pref = YYYYMMDDhhmm+'00-JAXA-L2P_GHRSST-SSTskin-H08_AHI_NRT-'+ \
        fname_pref = YYYYMMDDhhmm+'00-JAXA-L2P_GHRSST-SSTskin-H08_AHI-'+ \
                    VERd+'-v02.0-fv'+FVER+'.nc'
        if box is not None:
            fname_out = box['name']+'_'+YYYYMMDDhhmm+'.nc'
        try:
            # get ftp file names
            ftp.cwd(fdir)
            files_all = ftp.nlst()
            files = []
            for f in files_all:
                if fname_pref in f and f[-2:]=='nc':
                    files.append(f)
            #
            for fname in files:
                #
                fileraw = rpath+fname_pref
                sstfiles.append(fileraw)
                #
                if not os.path.isfile(fileraw):
                    print('No sst file found, download from jaxa ftp: '+fname)
                    lfile = open(rpath+fname,'wb')
                    try:
                        ftp.retrbinary('RETR %s' % fname, lfile.write)
                    except:
                        print('Error')
                    lfile.close()
                else:
                    print(fname_pref+' is already on disk')
                #
                if box is not None:
                    try:
                        process_raw_data(fileraw, box, fname_out)
                    except:
                        print('Cannot process '+fileraw)                
        except:
            print('Cannot access data in '+fdir)
    #
    ftp.quit()

    return sstfiles


def process_raw_data(fileraw,box,fname_out, \
                     variables=['sea_surface_temperature','wind_speed','solar_zenith_angle']):
    """ Open the raw file and zoom into box
    """
    ds = xr.open_dataset(fileraw).sel(lon=box['lon'],lat=box['lat'])
    #print(ds)
    for v in ds.data_vars:
        if v not in variables:
            ds = ds.drop(v)
    ds.to_netcdf(dpath+fname_out)
    print('  file reduced to box coordinates')


#
# Main part of the script
#
if __name__ == "__main__":

    box, tstart, delt = parse_input()
    print('box defined by following coordinates: ',box)
    print('start time: '+str(tstart))
    print('time interval considered: %.1f days'%delt)
    
    # output dir
    rpath = '/home/datawork-lops-osi/data/hw/sst/raw/'
    dpath = '/home/datawork-lops-osi/data/hw/sst/'+box['name']+'/'
    if not os.path.isdir(rpath): os.mkdir(rpath)
    if not os.path.isdir(dpath): os.mkdir(dpath)

    time = [tstart+datetime.timedelta(minutes=10*n) for n in range(int(delt*24*60/10))]

    # download files
    sstfiles = download_hw(time[:], rpath, dpath=dpath, box=box)



