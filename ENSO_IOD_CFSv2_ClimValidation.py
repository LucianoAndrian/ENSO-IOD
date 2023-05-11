"""
Validacion climatologica CFSv2 SON HGT200 contra ERA5
Periodo: 1981-2020 ?
Detrend ?
"""
################################################################################
import xarray as xr
import numpy as np
from ENSO_IOD_Funciones import SelectNMMEFiles
#------------------------------------------------------------------------------#
dir_hc = '/pikachu/datos/luciano.andrian/hindcast/'
dir_rt = '/pikachu/datos/luciano.andrian/real_time/'
out_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
v = 'hgt'
################################################################################
def fix_calendar(ds, timevar='time'):
    """
    agrega los dias a los archivos nc de NMME
    """
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds

def TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, main_month_season):

    for l in [0,1,2,3]:
        season_1982_1998 = data_1982_1998.sel(
            time=data_1982_1998.time.dt.month.isin(main_month_season - l), L=l)
        season_1999_2011 = data_1999_2011.sel(
            time=data_1999_2011.time.dt.month.isin(main_month_season - l), L=l)

        if l==0:
            season_clim_1982_1998 = season_1982_1998.mean(['r', 'time'])
            season_clim_1999_2011 = season_1999_2011.mean(['r', 'time'])

            season_anom_1982_1998 = season_1982_1998 - season_clim_1982_1998
            season_anom_1999_2011 = season_1999_2011 - season_clim_1999_2011
        else:
            season_clim_1982_1998 = xr.concat(
                [season_clim_1982_1998, season_1982_1998.mean(['r', 'time'])],
                dim='L')
            season_clim_1999_2011 = xr.concat(
                [season_clim_1999_2011, season_1999_2011.mean(['r', 'time'])],
                dim='L')

            aux_1982_1998 = season_1982_1998 - season_1982_1998.mean(
                ['r', 'time'])
            aux_1999_2011 = season_1999_2011 - season_1999_2011.mean(
                ['r', 'time'])

            season_anom_1982_1998 = xr.concat(
                [season_anom_1982_1998, aux_1982_1998], dim='time')
            season_anom_1999_2011 = xr.concat(
                [season_anom_1999_2011, aux_1999_2011], dim='time')

    return season_clim_1982_1998, season_clim_1999_2011, season_anom_1982_1998,\
           season_anom_1999_2011

def Detrend_Seasons(season_anom_1982_1998, season_anom_1999_2011,
                    main_month_season):

    for l in [0,1,2,3]:
        #1982-1998
        aux_season_anom_1982_1998 \
            = season_anom_1982_1998.sel(
            time=season_anom_1982_1998.time.dt.month.isin(main_month_season-l))

        aux = aux_season_anom_1982_1998.mean('r').polyfit(dim='time', deg=1)
        aux_trend = xr.polyval(aux_season_anom_1982_1998['time'],
                               aux.hgt_polyfit_coefficients)
        if l == 0:
            season_anom_1982_1998_detrened = \
                aux_season_anom_1982_1998 - aux_trend
        else:
            aux_detrend = aux_season_anom_1982_1998 - aux_trend
            season_anom_1982_1998_detrened = \
                xr.concat([season_anom_1982_1998_detrened, aux_detrend], dim='time')

    # 1999-2011
        aux_season_anom_1999_2011 \
            = season_anom_1999_2011.sel(
            time=season_anom_1999_2011.time.dt.month.isin(main_month_season - l))

        aux = aux_season_anom_1999_2011.mean('r').polyfit(dim='time', deg=1)
        aux_trend = xr.polyval(aux_season_anom_1999_2011['time'],
                               aux.hgt_polyfit_coefficients)
        if l==0:
            season_anom_1999_2011_detrend = aux_season_anom_1999_2011 - aux_trend
        else:
            aux_detrend = aux_season_anom_1999_2011 - aux_trend
            season_anom_1999_2011_detrend = \
                xr.concat([season_anom_1999_2011_detrend,
                           aux_detrend], dim='time')

    return season_anom_1982_1998_detrened, season_anom_1999_2011_detrend

def Anom_Detrend_SeasonRealTime(data_realtime, season_clim_1999_2011, main_month_season):

    for l in [0,1,2,3]:
        season_data = data_realtime.sel(
            time=data_realtime.time.dt.month.isin(main_month_season-l), L=l)
        aux_season_clim_1999_2011 = season_clim_1999_2011.sel(L=l)

        #Anomalia
        season_anom = season_data - aux_season_clim_1999_2011

        #Detrend
        aux = season_anom.mean('r').polyfit(dim='time', deg=1)
        aux_trend = xr.polyval(season_anom['time'],
                               aux.hgt_polyfit_coefficients)

        if l==0:
            season_anom_detrend = season_anom - aux_trend
        else:
            aux_detrend = season_anom - aux_trend
            season_anom_detrend = \
                xr.concat([season_anom_detrend, aux_detrend], dim='time')

    return season_anom_detrend

# CFSv2 ########################################################################
# Hindcast --------------------------------------------------------------------#
files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                        dir=dir_hc, All=True)
files = sorted(files, key=lambda x: x.split()[0])

# Abriendo todos los archivos
#xr no entiende la codificacion de Leads, r y las fechas
data = xr.open_mfdataset(files, decode_times=False)
data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
data = data.sel(L=[0.5, 1.5, 2.5, 3.5]) # Solo leads 0 1 2 3
data['L'] = [0,1,2,3]
data = xr.decode_cf(fix_calendar(data)) # corrigiendo fechas
data = data.sel(lat=slice(-80, 20))

#media movil de 3 meses para separar en estaciones
hindcast = data.rolling(time=3, center=True).mean()
hindcast = hindcast.load()
del data

# real - time -----------------------------------------------------------------#
files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                        dir=dir_rt, All=True)
files = sorted(files, key=lambda x: x.split()[0])

#abriendo todos los archivos
data = xr.open_mfdataset(files, decode_times=False)
data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
data = data.sel(L=[0.5, 1.5, 2.5, 3.5]) # Solo leads 0 1 2 3
data['L'] = [0,1,2,3]
data = xr.decode_cf(fix_calendar(data)) # corrigiendo fechas
data = data.sel(lat=slice(-80, 20), P=200)
data = data.drop('P')
data = data.drop('Z')

#media movil de 3 meses para separar en estaciones
realtime = data.rolling(time=3, center=True).mean()
realtime = realtime.load()
del data

# NO detrend ------------------------------------------------------------------#

def MakeSonClim(da):
    for l in [0,1,2,3]:
        da_son = da.sel(time=da.time.dt.month.isin(10 - l), L=l)
        if l == 0:
            da_son_clim  = da_son .mean(['r', 'time'])
        else:
            da_son_clim = xr.concat(
                [da_son_clim, da_son.mean(['r', 'time'])], dim='L')

    return da_son_clim.mean(['L'])

hind_son_clim = MakeSonClim(hindcast)
real_son_clim = MakeSonClim(realtime)

full = xr.concat([hindcast, realtime], dim='time')
full = MakeSonClim(full)


# detrend ---------------------------------------------------------------------#
# La tendencia se remueve igual que en "preSelect_HGT.py"

# 1982-1998, 1999-2011
data_1982_1998 = hindcast.sel(
    time=hindcast.time.dt.year.isin(np.linspace(1982,1998,17)))
data_1999_2011 = hindcast.sel(
    time=hindcast.time.dt.year.isin(np.linspace(1999,2011,13)))
# hindcast --------------------------------------------------------------------#
son_clim_82_98, son_clim_99_11, son_anom_82_98, son_anom_99_11 = \
    TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 10)
son_anom_82_98_detrend, son_anom_99_11_detrend = \
    Detrend_Seasons(son_anom_82_98, son_anom_99_11, 10)

son_hindcast_detrend = xr.concat(
    [son_anom_82_98_detrend, son_anom_99_11_detrend], dim='time')

# realtime --------------------------------------------------------------------#
son_realtime_detrend = Anom_Detrend_SeasonRealTime(realtime, son_clim_99_11, 10)

son_detrend = xr.concat([son_hindcast_detrend, son_realtime_detrend], dim='time')
# -----------------------------------------------------------------------------#
hind_son_clim_detrend = son_hindcast_detrend.mean(['time', 'r'])
real_son_clim_detrend = son_realtime_detrend.mean(['time', 'r'])
son_detrend_clim = son_detrend.mean(['time', 'r'])

# -----------------------------------------------------------------------------#
# ERA5 ------------------------------------------------------------------------#
# era5 NO
data_dir_era5 = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/' \
                'downloaded/'
data = xr.open_dataset(data_dir_era5 + 'ERA5_HGT200_40-20.nc')
data = data.rolling(time=3, center=True).mean()
data = data.sel(time=data.time.dt.month.isin(10))
data= data.rename({'z':'hgt', 'longitude':'lon', 'latitude':'lat'})

def Weights(data):
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180),
                                   (len(data.lon), 1)))
    data_w = data * weights
    return data_w

data = Weights(data)
data = data.sel(lat=slice(20, -80))

era = data.interp(lon=hindcast.lon.values, lat=hindcast.lat.values)
era_hind = era.sel(time=era.time.dt.year.isin(range(1981, 2012)))
era_real = era.sel(time=era.time.dt.year.isin(range(2011, 2020)))
#------------------------------------------------------------------------------#
# Detrend y Weights
data_dir_era5 = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/' \
                '1940_2020/'
data = xr.open_dataset(data_dir_era5 + 'HGT200_mer_d_w.nc')
data = data.sel(lat=slice(20, -80))
data = data.interp(lon=hindcast.lon.values, lat=hindcast.lat.values)
data = data.rolling(time=3, center=True).mean()
data = data.sel(time=data.time.dt.month.isin(10))

era_detrend = data.rename({'var':'hgt'})
era_hind_detrend = era.sel(time=era.time.dt.year.isin(range(1981, 2012)))
era_real_detrend = era.sel(time=era.time.dt.year.isin(range(2011, 2020)))
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# validacion
# NO detrend
def Weights2(data):
    weights = np.transpose(
        np.tile(np.cos(np.linspace(-80, 20, 101) * np.pi / 180),
                (len(data.lon), 1)))
    data_w = data * weights
    return data_w

# Hindcast
aux = era_hind.mean('time') - \
      Weights2(hind_son_clim.drop(['P']).__mul__(9.80665))

import matplotlib.pyplot as plt
plt.imshow(aux.hgt[:,:,0], vmin=-200, vmax=200);plt.colorbar();plt.show()


aux = era_real.mean('time') - \
      Weights2(real_son_clim.__mul__(9.80665))

import matplotlib.pyplot as plt
plt.imshow(aux.hgt, vmin=-200, vmax=200);plt.colorbar();plt.show()

