import xarray as xr
import numpy as np
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
################################################################################
dir_files = '/pikachu/datos/luciano.andrian/observado/ncfiles/data_no_detrend/'
out_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/'
# Funciones ####################################################################
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

def ChangeLons(data, lon_name='lon'):
    data['_longitude_adjusted'] = xr.where(
        data[lon_name] < 0,
        data[lon_name] + 360,
        data[lon_name])

    data = (
        data
            .swap_dims({lon_name: '_longitude_adjusted'})
            .sel(**{'_longitude_adjusted': sorted(data._longitude_adjusted)})
            .drop(lon_name))

    data = data.rename({'_longitude_adjusted': 'lon'})

    return data
################################################################################

# T CRU
data = xr.open_dataset(dir_files + 't_cru0.25.nc')
data = data.drop('stn')
data = data.rename({'tmp':'var'})
data = ChangeLons(data)
data_40_20 = data.sel(time=slice('1940-01-16', '2020-12-16'))
del data

data_40_20 = Weights(data_40_20)
data_40_20 = data_40_20.sel(lat=slice(-80, 20)) # HS
data_40_20 = data_40_20.rolling(time=3, center=True).mean()
for mm, s_name in zip([7,10], ['JJA','SON']):  # en caso de sumar otras...
    aux = data_40_20.sel(time=data_40_20.time.dt.month.isin(mm))
    aux = Detrend(aux, 'time')

    aux.to_netcdf(out_dir + 'tcru_w_c_d_0.25_'+ s_name + '.nc')

del data_40_20, aux

# PP GPCC
# 0.25 hay que descargarlos en gz por separado...
#https://opendata.dwd.de/climate_environment/GPCC/full_data_monthly_v2022/025/
#data = xr.open_dataset(dir_files + 'pp_gpcc_v2022_0.25.nc')
# probando el monitoring de 2020 actualizado hasta 2023
data = xr.open_dataset(dir_files + 'pp_pgcc_v2020_1891-2023_1.nc')
data = data.rename({'precip':'var'})
data_40_20 = data.sel(time=slice('1940-01-16', '2020-12-16'))
del data

data_40_20 = Weights(data_40_20)
data_40_20 = data_40_20.sel(lat=slice(20, -80)) # HS
data_40_20 = data_40_20.rolling(time=3, center=True).mean()
for mm, s_name in zip([7, 10], ['JJA','SON']):  # en caso de sumar otras...
    aux = data_40_20.sel(time=data_40_20.time.dt.month.isin(mm))
    aux = Detrend(aux, 'time')
    aux.to_netcdf(out_dir + 'ppgpcc_w_c_d_1_'+ s_name + '.nc')

print('done')

#
# aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/
# data_obs_viejo/pp_20CR-V3.nc')
# pp_20cr = aux.sel(lon=slice(270, 330), lat=slice(-60, 15))
# data_20_20_int = data_20_20.interp(lon=pp_20cr.lon.values,
# lat=pp_20cr.lat.values)
# data_50_20_int = data_50_20.interp(lon=pp_20cr.lon.values,
# lat=pp_20cr.lat.values)
#
# data_20_20_int.to_netcdf(out_dir + 't_cru_d_w_c_1920-2020_1.nc')
# data_50_20_int.to_netcdf(out_dir + 't_cru_d_w_c_1950-2020_1.nc')
#
# # PP GPCC
# data = xr.open_dataset(dir_files + 'pp_gpcc_v2020_0.25.nc')
# data = data.rename({'precip':'var'})
#
# data_20_20 = data.sel(time=slice('1920-01-16', '2020-12-16'))
# data_50_20 = data.sel(time=slice('1950-01-16', '2020-12-16'))
#
# del data
#
# data_20_20 = Weights(data_20_20)
# data_50_20 = Weights(data_50_20)
#
# data_20_20 = data_20_20.sel(lon=slice(270, 330), lat=slice(15, -60))
# data_50_20 = data_50_20.sel(lon=slice(270, 330), lat=slice(15, -60))
#
# data_20_20 = Detrend(data_20_20, 'time')
# data_50_20 = Detrend(data_50_20, 'time')
#
# data_20_20.to_netcdf(out_dir + 'pp_gpcc_d_w_c_1920-2020_0.25.nc')
# data_50_20.to_netcdf(out_dir + 'pp_gpcc_d_w_c_1950-2020_0.25.nc')
#
# aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/
# data_obs_viejo/pp_20CR-V3.nc')
# pp_20cr = aux.sel(lon=slice(270, 330), lat=slice(-60, 15))
# data_20_20_int = data_20_20.interp(lon=pp_20cr.lon.values, l
# at=pp_20cr.lat.values)
# data_50_20_int = data_50_20.interp(lon=pp_20cr.lon.values,
# lat=pp_20cr.lat.values)
#
# data_20_20_int.to_netcdf(out_dir + 'pp_gpcc_d_w_c_1920-2020_1.nc')
# data_50_20_int.to_netcdf(out_dir + 'pp_gpcc_d_w_c_1950-2020_1.nc')
#
#
# #PP PREC:
# data = xr.open_dataset(dir_files + 'precip.mon.anom.nc')
# data = data.rename({'precip':'var'})
# data_50_20 = data.sel(time=slice('1950-01-01', '2020-12-01'))
# data_50_20 *= 30
# del data
#
# data_50_20 = Weights(data_50_20)
# data_50_20 = data_50_20.sel(lon=slice(270, 330), lat=slice(15, -60))
# data_50_20 = Detrend(data_50_20, 'time')
# data_50_20.to_netcdf(out_dir + 'pp_prec_d_w_c_1950-2020_2.5.nc')
#
#
#
# #PP CMAP:
# data = xr.open_dataset(dir_files + 'pp_cmap.nc')
# data = data.rename({'precip':'var'})
# data_50_20 = data.sel(time=slice('1950-01-01', '2020-12-01'))
# data_50_20 *= 30
# del data
#
# data_50_20 = Weights(data_50_20)
# data_50_20 = data_50_20.sel(lon=slice(270, 330), lat=slice(15, -60))
# data_50_20 = Detrend(data_50_20, 'time')
# data_50_20.to_netcdf(out_dir + 'pp_cmap_d_w_c_1979-2020_2.5.nc')
#
# ##############################################################################
# # Prueba hemisferio
# ##############################################################################
#
# # T CRU
# data = xr.open_dataset(dir_files + 't_cru0.25.nc')
# data = data.drop('stn')
# data = data.rename({'tmp':'var'})
# data = ChangeLons(data)
#
# #data_20_20 = data.sel(time=slice('1920-01-16', '2020-12-16'))
# data_50_20 = data.sel(time=slice('1950-01-16', '2020-12-16'))
#
# del data
#
# #data_20_20 = Weights(data_20_20)
# data_50_20 = Weights(data_50_20)
#
# #data_20_20 = data_20_20.sel(lat=slice(-90, 20))
# data_50_20 = data_50_20.sel(lat=slice(-90, 20))
#
# #data_20_20 = Detrend(data_20_20, 'time')
# data_50_20 = Detrend(data_50_20, 'time')
#
# #data_20_20.to_netcdf(out_dir + 't_cru_d_w_c_1920-2020_0.25.nc')
# data_50_20.to_netcdf(out_dir + 't_cru_HS_d_w_c_1950-2020_0.25.nc')
#
#
