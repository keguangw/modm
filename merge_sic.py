import sys, os
import numpy as np
from datetime import datetime, timedelta
from netCDF4 import Dataset as nc

#===============================================================================
app  = 'barents-2.5km'
wdir = os.environ.get('METROMS_TMPDIR') + '/' + app + '/data/'
grd  = os.environ.get('METROMS_APPDIR') + '/' + app + '/grid/'

tname  = 'time'
nodata = -999
snodata = str(nodata)

#===========================================================================
def read_data(sensor):

    fname0 = wdir + sensor + '_sic_' + date + '.nc'
    if os.path.exists(fname0):
       print('  Reading ' + sensor + ' SIC ...')
       f0 = nc(fname0)
       sic0 = f0.variables['obs_sic'][0,:,:]
       err0 = f0.variables['obs_err'][0,:,:]
    else:
       print('  ' + fname0 + ' not available ...')
       sic0, err0 = None, None

    return sic0, err0

#===========================================================================
def merge(sic1,err1,sic2,err2,mask):

#   calculating precision p = 1/err**2
    eps = 1.0e-16
    ny, nx = sic1.shape
    sic = np.zeros((ny,nx)) + nodata
    err = np.zeros((ny,nx)) + nodata

    sic1[~np.isfinite(sic1)] = nodata
    err1[~np.isfinite(err1)] = nodata
    sic2[~np.isfinite(sic2)] = nodata
    err2[~np.isfinite(err2)] = nodata

    for iy in range(ny):
        for ix in range(nx):
            if mask[iy,ix] != 0:
               if ((sic1[iy,ix] <= 100) & (sic1[iy,ix] >= 0)):
                  p1 = 1./(err1[iy,ix]**2 + eps)
               else:
                  p1 = 0.
                  sic1[iy,ix] = 0.

               if ((sic2[iy,ix] <= 100.) & (sic2[iy,ix] >= 0)):
                  p2 = 1./(err2[iy,ix]**2 + eps)
               else:
                  p2 = 0
                  sic2[iy,ix] = 0.

               p = p1 + p2
               if p > 0:
                  err[iy,ix] = np.sqrt(1.0/p)
                  sic[iy,ix] = (sic1[iy,ix]*p1 + sic2[iy,ix]*p2) / p

    return sic, err


#===========================================================================
def write_dim(f,Lon,Lat,dt):

    print ('  writing dimensions ...')

    f.createDimension('time',None)
    time = f.createVariable('time','f8',('time',))
    time.long_name = 'time'
    time.units = 'days since 1970-01-01 00:00:00'
    time[:] = dt.astype(np.float64)

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
def write_var(f,var,var_name,unit,long_name,lonname='lon', \
              latname='lat',ks=-1):

    print ('  saving ' + var_name + str(var.shape) + ' ...')

    try:
       varid = f.variables[var_name]
       print ('    shape of old data:', varid.shape)
       print ('    shape of new data:', var.shape)

       kn = var.shape[0]
       varid[ks:ks+kn,...] = var.astype(np.float32)
    except:
       varid = f.createVariable(var_name,'f',(tname,latname,lonname))

       if len(var.shape) == 2:
          ny, nx = var.shape
          var1 = np.zeros([1,ny,nx])
          var1[0,:,:] = var
          varid[:] = var1.astype(np.float32)
       else:
          varid[:] = var.astype(np.float32)
       varid.long_name = long_name
       varid.units = unit
       varid.time = tname
       varid.missing_value = snodata

    f.sync()

#===============================================================================
if len(sys.argv) <= 1:
   date = raw_input('Input date: e.g. 20100101\n')
elif len(sys.argv) == 2:
   date = sys.argv[1]

#===============================================================================
# read model grid file
fgrd = nc(grd + '/' + app + '_grd.nc','r')
msk = fgrd.variables['mask_rho'][:]
Lon = fgrd.variables['lon_rho'][:]
Lat = fgrd.variables['lat_rho'][:]
ny, nx = msk.shape

#===============================================================================
fout = wdir + 'merged_sic_' + date + '.nc'
f = nc(fout,'w')
t0 = datetime(1970,1,1)
t1 = datetime.strptime(date,'%Y%m%d')
dt = np.array([(t1-t0).days])
write_dim(f,Lon,Lat,dt)

#===============================================================================
print('Reading input data ...')
sic1, err1 = read_data('amsr2')
sic2, err2 = read_data('icechart')

if (sic1 is not None) and (sic2 is not None):
   sic, err = merge(sic1,err1,sic2,err2,msk)
elif sic1 is not None:
   sic, err = sic1, err1
elif sic2 is not None:
   sic, err = sic2, err2
else:
   print('Both data are missing, exit ...')
   sys.exis()

#===============================================================================
print ('Output results to ' + fout + '...')
write_var(f,sic,'obs_sic','%','merged SIC')
write_var(f,err,'obs_err','%','error of merged SIC')

f.close()


