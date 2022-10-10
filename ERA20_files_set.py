import xarray as xr
import numpy as np
from ENSO_IOD_Funciones import xrFieldTimeDetrend

dir_files = '/datos/luciano.andrian/ncfiles/ERA20C/'
variables =['HGT200', 'psl', 'UV200']
name_variables = ['z', 'msl']

n = 0
for v in variables:
    n_v = name_variables[n]
    print(v)
    data = xr.open_dataset(dir_files + 'ERA-20C_' + v + '.nc')

    data = data.rename({n_v: 'var'})
    data = data.rename({'longitude': 'lon'})
    data = data.rename({'latitude': 'lat'})

    data20_49 = data.sel(time=slice('1920-01-01', '1949-12-31'))

    print('detrend...')
    data = xrFieldTimeDetrend(data, 'time')
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180), (len(data.lon), 1)))
    data = data * weights
    print('to_netcdf...')
    data.to_netcdf(dir_files + 'ERA20_' + v + '_full.nc')

    print('detrend...')
    data20_49 = xrFieldTimeDetrend(data20_49, 'time')
    weights = np.transpose(np.tile(np.cos(data20_49.lat * np.pi / 180), (len(data20_49.lon), 1)))
    data20_49 = data20_49 * weights
    print('to_netcdf...')
    data20_49.to_netcdf(dir_files + 'ERA20_' + v + '_20-49.nc')
    n += 1
