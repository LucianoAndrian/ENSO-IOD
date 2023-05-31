import xarray as xr
import numpy as np
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
###############################################################################
#dir_files = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/merged/'
#out_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'
dir_files = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/downloaded/'
out_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/1940_2020/'
# Funciones ###################################################################
def Detrend(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
    dt = xrda - trend
    return dt

def Weights(data):
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180),
                                   (len(data.lon), 1)))
    data_w = data * weights
    return data_w

###############################################################################
# U V va por RWS_makerterms.py
variables =['HGT200', 'HGT750']
name_variables = ['z', 'z']

for v, n_v in zip(variables, name_variables):

    print(v)
    data = xr.open_dataset(dir_files + 'ERA5_' + v + '_40-20.nc')

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
    data = data.rolling(time=3, center=True).mean()
    for mm, s_name in zip([10], ['SON']): # main month seasons
        aux = data.sel(time=data.time.dt.month.isin(mm))
        aux = Detrend(aux, 'time')

        print('to_netcdf...')
        if v == 'UV200':
            aux.to_netcdf(out_dir + n_v + '_' + s_name + '_mer_d_w.nc')
        else:
            aux.to_netcdf(out_dir + v + '_' + s_name + '_mer_d_w.nc')

###############################################################################