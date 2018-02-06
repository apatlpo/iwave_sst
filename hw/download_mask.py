# coding: utf-8

# 
# # Play with JAXA Himawari 
#  
# In order to search for clear skies
# https://worldview.earthdata.nasa.gov/
# 
# ## ftp retrieval
# 
# 
# ## ftp structure at jaxa
# 
# ```
### /jma/hsd : higher resolution (500m), full disk+regions
#
# /jma/hsd
#       +---/[YYYYMM]
#              +---/[DD]
#                     +---/[hh]
# Full-disk
# HS_H08_YYYYMMDD_hhmm_Bbb_FLDK_Rjj_Skkll.DAT
# where YYYY: 4-digit year of observation start time (timeline);
#       MM: 2-digit month of timeline;
#       DD: 2-digit day of timeline;
#       hh: 2-digit hour of timeline;
#       mm: 2-gidit minutes of timeline;
#       bb: 2-digit band number (varies from "01" to "16");
#       jj: spatial resolution ("05": 0.5km, "10": 1.0km, "20": 2.0km);
#       kk: segment number (varies from "01" to "10"); and
#       ll: total number of segments (fixed to "10").
#
# examples:
# HS_H08_20150728_2200_B01_FLDK_R10_S0110.DAT
# HS_H08_20150728_1000_B03_FLDK_R05_S0510.DAT
#  
# ```
# 


# Import libraries
import os, sys
import ftplib as ftplib
import datetime


#
# download data
#
def download_hw(times):

    # open ftp
    host, login, passwd = ftp_login()
    ftp = ftplib.FTP(host)
    ftp.login(login, passwd)

    # loop around times
    maskfiles=[]
    for time in times:

        # format date and time
        YYYYMMDDhhmm = time.strftime('%Y%m%d%H%M')
        YYYYMM = time.strftime('%Y%m')
        DD = time.strftime('%d')
        hh = time.strftime('%H')
        mm = time.strftime('%M')

        # get mask data
        fdir='/pub/himawari/L2/CLP/bet/'+YYYYMM+'/'+DD+'/'+hh+'/'
        fname_pref='NC_H08_'+YYYYMM+DD+'_'+hh+mm+'_L2CLPbet_FLDK.02401_02401.nc'

        try:
            # get ftp file names
            ftp.cwd(fdir)
            files_all = ftp.nlst()
            files = []
            for f in files_all:
                if fname_pref in f:
                    files.append(f)
            #
            for fname in files:
                #
                fileraw = dpath+fname_pref
                print(fname)
                maskfiles.append(fileraw)
                #
                if not os.path.isfile(fileraw):
                    print('No cloud mask file found, download from jaxa ftp: '+fname)
                    lfile = open(dpath+fname,'wb')
                    try:
                        ftp.retrbinary('RETR %s' % fname, lfile.write)
                    except:
                        print('Error')
                    lfile.close()
                else:
                    print(fileraw+' is already on disk')
        except:
            print('Cannot access data in '+fdir)
            
            
    ftp.quit()

    return maskfiles


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
# Main part of the script
#



if __name__ == "__main__":

    # output dir
    #dpath = '/home/slyne/aponte/sst_fast/data_bin/';
    #dpath = '/home1/scratch/aponte/hw/mask/';
    dpath = '/home/datawork-lops-osi/data/hw/mask/';

    # time line
    # 201603 ->201706
    t1=datetime.datetime(2016,3,1,0,0,0);
    t1=datetime.datetime(2017,8,4,15,0,0);
    delt = 30.; # in days
    dt=30 # in minutes
    #
    #time = [t1+datetime.timedelta(minutes=dt*n) for n in xrange(int(delt*24*60/dt))]
    #time = [t1+datetime.timedelta(minutes=dt*n) for n in xrange(1)]
    #
    t2=datetime.datetime.today()
    t = t1
    time = [t]
    while t<=t2:
        t += datetime.timedelta(minutes=dt)
        time.append(t)
    print(time[0])
    print(time[-1])
    
    maskfiles = download_hw(time)
    
    print('---\n All files downloaded \n---')

    sys.exit()
