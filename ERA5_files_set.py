import xarray as xr
import numpy as np
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
########################################################################################################################
dir_files = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/merged/'
out_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'
# Funciones ############################################################################################################
def Detrend(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
    dt = xrda - trend
    return dt

def Weights(data):
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180), (len(data.lon), 1)))
    data_w = data * weights
    return data_w

########################################################################################################################
variables =['hgt200', 'div200', 'slp', 'vp200', 'vp500', 'UV200', 'UV200']
name_variables = ['z', 'd', 'msl', 'w', 'w', 'u', 'v']

variables =['hgt750']
name_variables = ['z']

n = 0
for v in variables:
    n_v = name_variables[n]

    print(v)
    data = xr.open_dataset(dir_files + 'ERA5_' + v + '_50-20_mer.nc')

    if n_v == 'u':
        print('Drop v')
        data = data.drop('v')
    elif n_v == 'v':
        print('Drop u')
        data = data.drop('u')

    data = data.rename({n_v: 'var'})
    data = data.rename({'longitude': 'lon'})
    data = data.rename({'latitude': 'lat'})

    data = Weights(data)
    data = data.sel(lat=slice(20, -80))
    data = Detrend(data, 'time')

    print('to_netcdf...')
    if v == 'UV200':
        data.to_netcdf(out_dir + n_v + '_mer_d_w.nc')
    else:
        data.to_netcdf(out_dir + v + '_mer_d_w.nc')

    n += 1
########################################################################################################################
#probando hemisferio
########################################################################################################################
# variables =['hgt200', 'div200', 'slp', 'vp200', 'vp500', 'UV200', 'UV200']
# name_variables = ['z', 'd', 'msl', 'w', 'w', 'u', 'v']

n = 0
#for v in variables:
v = 'div200'
n_v='d'
    #n_v = name_variables[n]

print(v)
data = xr.open_dataset(dir_files + 'ERA5_' + v + '_50-20_mer.nc')
if n_v == 'u':
    print('Drop v')
    data = data.drop('v')
elif n_v == 'v':
    print('Drop u')
    data = data.drop('u')

data = data.rename({n_v: 'var'})
data = data.rename({'longitude': 'lon'})
data = data.rename({'latitude': 'lat'})

data = Weights(data)
data = data.sel(lat=slice(20, -90))
data = Detrend(data, 'time')

print('to_netcdf...')
if v == 'UV200':
    data.to_netcdf(out_dir + n_v + '_HS_mer_d_w.nc')
else:
    data.to_netcdf(out_dir + v + '_HS_mer_d_w.nc')





    n += 1