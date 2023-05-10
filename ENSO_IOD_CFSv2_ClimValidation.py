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

    return season_clim_1982_1998, season_clim_1999_2011, season_anom_1982_1998, season_anom_1999_2011

def Detrend_Seasons(season_anom_1982_1998, season_anom_1999_2011, main_month_season):

    for l in [0,1,2,3]:
        #1982-1998
        aux_season_anom_1982_1998 \
            = season_anom_1982_1998.sel(time=season_anom_1982_1998.time.dt.month.isin(main_month_season-l))

        aux = aux_season_anom_1982_1998.mean('r').polyfit(dim='time', deg=1)
        aux_trend = xr.polyval(aux_season_anom_1982_1998['time'], aux.hgt_polyfit_coefficients)
        if l == 0:
            season_anom_1982_1998_detrened = aux_season_anom_1982_1998 - aux_trend
        else:
            aux_detrend = aux_season_anom_1982_1998 - aux_trend
            season_anom_1982_1998_detrened = xr.concat([season_anom_1982_1998_detrened, aux_detrend], dim='time')

    # 1999-2011
        aux_season_anom_1999_2011 \
            = season_anom_1999_2011.sel(time=season_anom_1999_2011.time.dt.month.isin(main_month_season - l))

        aux = aux_season_anom_1999_2011.mean('r').polyfit(dim='time', deg=1)
        aux_trend = xr.polyval(aux_season_anom_1999_2011['time'], aux.hgt_polyfit_coefficients)
        if l==0:
            season_anom_1999_2011_detrend = aux_season_anom_1999_2011 - aux_trend
        else:
            aux_detrend = aux_season_anom_1999_2011 - aux_trend
            season_anom_1999_2011_detrend = xr.concat([season_anom_1999_2011_detrend, aux_detrend], dim='time')

    return season_anom_1982_1998_detrened, season_anom_1999_2011_detrend

def Anom_Detrend_SeasonRealTime(data_realtime, season_clim_1999_2011, main_month_season):

    for l in [0,1,2,3]:
        season_data = data_realtime.sel(time=data_realtime.time.dt.month.isin(main_month_season-l), L=l)
        aux_season_clim_1999_2011 = season_clim_1999_2011.sel(L=l)

        #Anomalia
        season_anom = season_data - aux_season_clim_1999_2011

        #Detrend
        aux = season_anom.mean('r').polyfit(dim='time', deg=1)
        aux_trend = xr.polyval(season_anom['time'], aux.hgt_polyfit_coefficients)

        if l==0:
            season_anom_detrend = season_anom - aux_trend
        else:
            aux_detrend = season_anom - aux_trend
            season_anom_detrend = xr.concat([season_anom_detrend, aux_detrend], dim='time')

    return season_anom_detrend
# CFSv2 -----------------------------------------------------------------------#
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
data = data.rolling(time=3, center=True).mean()

# 1982-1998, 1999-2011
data_1982_1998 = data.sel(time=data.time.dt.year.isin(np.linspace(1982,1998,17)))
data_1999_2011 = data.sel(time=data.time.dt.year.isin(np.linspace(1999,2011,13)))
#----------------------------------------------------------------------------------------------------------------------#
son_clim_82_98, son_clim_99_11, son_anom_82_98, son_anom_99_11 = \
    TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 10)
son_anom_82_98_detrend, son_anom_99_11_detrend = \
    Detrend_Seasons(son_anom_82_98, son_anom_99_11, 10)

son_hindcast_detrend = xr.concat([son_anom_82_98_detrend, son_anom_99_11_detrend], dim='time')






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
data = data.rolling(time=3, center=True).mean()

for l in [0,1,2,3]:
    son = data.sel(
        time=data.time.dt.month.isin(10 - l), L=l)
    if l == 0:
        son_clim = son.mean(['r', 'time'])
    else:
        son_clim = xr.concat([son_clim, son.mean(['r', 'time'])], dim='L')





# ERA5 ------------------------------------------------------------------------#
data_dir_era5 = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/' \
                '1940_2020/'

era = xr.open_dataset(data_dir_era5 + 'HGT200_mer_d_w.nc')
era = era.sel(lat=slice(20, -80))
era = era.interp(lon=data.lon.values, lat=data.lat.values)
era = era.rolling(time=3, center=True).mean()
era = era.sel(time=era.time.dt.month.isin(10))
era = era.sel(time=era.time.dt.year.isin(range(1981, 2012)))
era = era.rename({'var':'hgt'})


def Weights(data):
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180),
                                   (len(data.lon), 1)))
    data_w = data * weights
    return data_w

aux = Weights(era.mean('time')) - \
      Weights(son_hindcast_detrend.mean(['r', 'time']).drop(['P']).__mul__(9.80665))

import matplotlib.pyplot as plt
plt.imshow(aux.hgt[:,:,0], vmin=-200, vmax=200);plt.colorbar();plt.show()


# ok da un poco mejor...



