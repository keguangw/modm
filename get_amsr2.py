#===============================================================================
#
# This script download and interpolate amsr2 6.25km product
# from University of Bremen
# The uncertainty is calculated following Spreen et al. (2008), JGR
#
# Author: Keguang Wang, keguang.wang@met.no
# 25/02/2020
#
#===============================================================================
import sys, os
import numpy as np
import subprocess as sp
from pyproj import Proj
from datetime import datetime, timedelta
from netCDF4 import Dataset as nc

app  = 'barents-2.5km'
wdir = os.environ.get('METROMS_TMPDIR') + '/' + app + '/data'
grd  = os.environ.get('METROMS_APPDIR') + '/' + app + '/grid'

# set missing_value for output
set_missing_value = True

nodata  = -999
snodata = str(nodata)
mon  = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

#===============================================================================
def download_amsr2(date):

    web0 = 'https://seaice.uni-bremen.de/data/amsr2/asi_daygrid_swath/n6250/'
    web  = web0 + 'netcdf/' + date[:4] + '/'
    fil  = 'asi-AMSR2-n6250-' + date + '-v5.4.nc'
    f1   =  wdir + '/amsr2/sic_' + date  + '.nc'

    if os.path.exists(f1) is False:
       try:
          print ('download AMSR2 ice concentration ...')
          output = sp.check_call(['wget', web + fil])
          if output == 0:
             sp.call(['mv', fil, f1])
       except:
           print ('AMSR2 SIC data not available on ' + date + ' !!!')
           sys.exit()

    return f1

#===========================================================================
def nearest(X0,Y0,var0,X,Y,mask=None):

    print ('  interpolating data with nearest-neighbour method ...')

    #define the origin in old coordinate, and steps
    x0, y0 = X0[0,0], Y0[0,0]
    dx = X0[0,1] - X0[0,0]
    dy = Y0[1,0] - Y0[0,0]
    dd = dx * dy

    # interpolation
    i1 = np.floor((X-x0) / dx + 0.5).astype(int)
    i2 = np.floor((Y-y0) / dy + 0.5).astype(int)

    i1[i1 >= X0.shape[1]] = X0.shape[1] - 1
    i2[i2 >= X0.shape[0]] = X0.shape[0] - 1

    var = var0[...,i2,i1]

    if mask is not None:
       var[...,mask==0] = np.nan
    
    return var     

#================================================================
def err_amsr(sic0,mask):

    print ('  calculating uncertainties for amsr2 data ...')

    sic = sic0 + 0.0
    sic[sic == nodata] = 0
    sic = sic * 0.01

    # estimate uncertainty (see Spreen et al., 2008)
    Psw,   Epsw    = 82.0, 4.0
    Psi,   Epsi	   = 10.0, 4.0
    Tauw,  Etauw   = 0.27, 0.1
    Taui,  Etaui   = 0.14, 0.035
    #d3, d2, d1, d0 = 1.64e-5, -1.6e-3, 1.92e-2, 0.971
    d3, d2, d1, d0 = 5.587e-06, -5.218e-04, -1.226e-02, 1.116

    Ps   = sic * Psi   + (1-sic) * Psw
    Tau  = sic * Taui  + (1-sic) * Tauw
    Etau = sic * Etaui + (1-sic) * Etauw
    ac   = 1.1*np.exp(-2*Tau) - 0.11*np.exp(-Tau)
    P    = Ps * ac

    Ep2 = (Ps*Etau*(0.11*np.exp(-Tau)-2.2*np.exp(-2*Tau)))**2 + \
          (ac*(1-sic)*Epsw)**2 + (ac*sic*Epsi)**2
    err = np.abs(3*d3*P**2 + 2*d2*P + d1) * np.sqrt(Ep2) * 100
    err[mask == 0] = nodata
    err[np.isnan(sic0)] = nodata

    return err

#===========================================================================
def write_dim(f,Lon,Lat,dt,ensemble=[]):

    print ('  writing dimensions ...')

    f.createDimension('time',None)
    time = f.createVariable('time','f8',('time',))
    time.long_name = 'time'
    time.units = 'days since 1970-01-01 00:00:00'
    time[:] = dt.astype(np.float64)

    ne = len(ensemble)
    if ne > 0:
       f.createDimension('number',ne)
       ensid = f.createVariable('number','i4',('number',))
       ensid[:] = ensemble

    if len(Lon.shape) == 2:
       ny, nx = Lon.shape
       f.createDimension('lon',nx)
       f.createDimension('lat',ny)
       lonid = f.createVariable('longitude','f',('lat','lon'))
       latid = f.createVariable('latitude','f',('lat','lon'))
    elif len(Lon.shape) == 1:
       ny, = Lat.shape
       nx, = Lon.shape
       f.createDimension('lon',nx)
       f.createDimension('lat',ny)
       lonid = f.createVariable('longitude','f',('lon',))
       latid = f.createVariable('latitude','f',('lat',))
    else:
       print ('geographical dimenion definition no supported !!!')
       sys.exit()

    lonid.units = 'degree East'
    lonid[:] = Lon.astype(np.float32)

    latid.units = 'degree North'
    latid[:] = Lat.astype(np.float32)

#===========================================================================
def write_var(f,var,var_name,unit,long_name,ks=-1):

    print ('  saving ' + var_name + str(var.shape) + ' ...')

    try:
       varid = f.variables[var_name]
       print ('    shape of old data:', varid.shape)
       print ('    shape of new data:', var.shape)

       kn = var.shape[0]
       varid[ks:ks+kn,...] = var.astype(np.float32)
    except:
       if len(var.shape) <= 3:
          varid = f.createVariable(var_name,'f',('time','lat','lon'))
       elif len(var.shape) == 4:
          varid = f.createVariable(var_name,'f',('time','number','lat','lon'))

       if len(var.shape) == 2:
          ny, nx = var.shape
          var1 = np.zeros([1,ny,nx])
          var1[0,:,:] = var
          varid[:] = var1.astype(np.float32)
       else:
          varid[:] = var.astype(np.float32)
       varid.long_name = long_name
       varid.units = unit
       varid.time = 'time'
       varid.missing_value = snodata

    f.sync()

#===============================================================================
if len(sys.argv) <= 1:
   date = raw_input('Input date: e.g. 20100101\n')
elif len(sys.argv) == 2:
   date = sys.argv[1]

#===============================================================================
print ('Downloading AMSR2 SIC for ' + date + ' ...')
fname = download_amsr2(date)

#===============================================================================
print ('Reading basic information from original sic file ...')
f1 = nc(fname,'r')
x0 = f1.variables['x'][:]
y0 = f1.variables['y'][:]
X0, Y0 = np.meshgrid(x0,y0)
sic0 = f1.variables['z'][:]

#===============================================================================
print ('Projecting to original AMSR2 grid ...')
p = Proj(proj='stere',lon_0=-45,lat_0=90,lat_ts=70,ellps='WGS84')

fgrd = nc(grd + '/barents_grd.nc','r')
msk = fgrd.variables['mask_rho'][:]
Lon = fgrd.variables['lon_rho'][:]
Lat = fgrd.variables['lat_rho'][:]
X, Y = p(Lon,Lat,inverse=False)
ny, nx = msk.shape

#===============================================================================
print ('Preparing output file ...')
f = nc(wdir + '/amsr2_sic_' + date + '.nc','w')

t0 = datetime(1970,1,1)
t1 = datetime.strptime(date,'%Y%m%d')
dt = np.array([(t1-t0).days]) + 0.5

write_dim(f,Lon,Lat,dt)

#===============================================================================
print ('  interpolating sic ...')
sic = nearest(X0,Y0,sic0,X,Y,msk)
err = err_amsr(sic,msk)
write_var(f,sic,'obs_sic','%','sea ice concentration',set_missing_value)
write_var(f,err,'obs_err','%','error of SIC',set_missing_value)

f.close()


