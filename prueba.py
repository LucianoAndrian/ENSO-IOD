import xarray as xr
import numpy as np

pwd = '/datos/luciano.andrian/ncfiles/'

seasons = ['DJF','JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA',
           'JAS', 'ASO', 'SON', 'OND', 'NDJ']

out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/'


variables = ['psl', 'hgt200','sf', 'div', 'vp']
var_name = ['psl', 'z', 'streamfunction', 'divergence','velocity_potential']
title_var = ['PSL', 'HGT', 'Psi', 'Divergence', 'Potential Velocity']


v = 2
var = xr.open_dataset(pwd + variables[v] + '.nc')
print(variables[v] + '.nc')

var = var.rename({'streamfunction': 'var'})
var = var.rename({'longitude': 'lon'})
var = var.rename({'latitude': 'lat'})

# Detrend
var = var - xr.polyval(var['time'], var.polyfit(dim='time', deg=1).var_polyfit_coefficients[0])
print('listo detrend')
weights = np.transpose(np.tile(np.cos(var.lat * np.pi / 180), (len(var.lon), 1)))
var = var * weights
print('writing nc')
var.to_netcdf(variables[v] + '.nc')

