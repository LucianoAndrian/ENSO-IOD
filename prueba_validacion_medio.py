import xarray as xr
import numpy as np
from ENSO_IOD_Funciones import SelectNMMEFiles
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
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

# Funciones para CFSv2
def fix_calendar(ds, timevar='time'):
    """
    agrega los dias a los archivos nc de NMME
    """
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds

def TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, main_month_season):

    for l in [0,1,2,3]:
        season_1982_1998 = data_1982_1998.sel(time=data_1982_1998.time.dt.month.isin(main_month_season-l), L=l)
        season_1999_2011 = data_1999_2011.sel(time=data_1999_2011.time.dt.month.isin(main_month_season-l), L=l)

        if l==0:
            season_clim_1982_1998 = season_1982_1998.mean(['r', 'time'])
            season_clim_1999_2011 = season_1999_2011.mean(['r', 'time'])

            season_anom_1982_1998 = season_1982_1998 - season_clim_1982_1998
            season_anom_1999_2011 = season_1999_2011 - season_clim_1999_2011
        else:
            season_clim_1982_1998 = xr.concat([season_clim_1982_1998, season_1982_1998.mean(['r', 'time'])], dim='L')
            season_clim_1999_2011 = xr.concat([season_clim_1999_2011, season_1999_2011.mean(['r', 'time'])], dim='L')

            aux_1982_1998 = season_1982_1998 - season_1982_1998.mean(['r', 'time'])
            aux_1999_2011 = season_1999_2011 - season_1999_2011.mean(['r', 'time'])

            season_anom_1982_1998 = xr.concat([season_anom_1982_1998, aux_1982_1998], dim='time')
            season_anom_1999_2011 = xr.concat([season_anom_1999_2011, aux_1999_2011], dim='time')

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
########################################################################################################################
# ERA5
dir_files = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/merged/'
obs = xr.open_dataset(dir_files + 'ERA5_' + 'hgt200' + '_50-20_mer.nc')
obs = obs.rename({'z': 'var'})
obs = obs.rename({'longitude': 'lon'})
obs = obs.rename({'latitude': 'lat'})
#obs = obs.rolling(time=3, center=True).mean()
#----------------------------------------------------------------------------------------------------------------------#

aux = obs.sel(time=obs.time.dt.year.isin(np.arange(1982,2021)))

#obs_proc = obs_proc.sel(time=obs.time.dt.month.isin(10))
obs_proc = Detrend(aux, 'time')
obs_proc = obs_proc.sel(time=obs_proc.time.dt.month.isin(10))

obs_proc = Weights(obs_proc)
obs_proc = obs_proc.sel(lat=slice(20, -90))

obs_proc = obs_proc.mean('time')


# CFSv2
dir_hc = '/pikachu/datos/luciano.andrian/hindcast/'
dir_rt = '/pikachu/datos/luciano.andrian/real_time/'
out_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
v = 'hgt'
# HINDCAST ------------------------------------------------------------------------------------------------------------#
files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                        dir=dir_hc, All=True)
files = sorted(files, key=lambda x: x.split()[0])

#abriendo todos los archivos
data = xr.open_mfdataset(files, decode_times=False) #xr no entiende la codificacion de Leads, r y las fechas
data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
data = data.sel(L=[0.5, 1.5, 2.5, 3.5]) # Solo leads 0 1 2 3
data['L'] = [0,1,2,3]
data = xr.decode_cf(fix_calendar(data)) # corrigiendo fechas
data = data.sel(P=200)
data = data.drop('P')
#data = data.drop('Z')

#media movil de 3 meses para separar en estaciones
#data = data.rolling(time=3, center=True).mean()

# for l in [0, 1, 2, 3]:
#     season = data.sel(time=data.time.dt.month.isin(10 - l), L=l)
#     if l == 0:
#         season_clim = season.mean(['r', 'time'])
#     else:
#         season_clim = xr.concat([season_clim, season.mean(['r', 'time'])], dim='L')

#season_clim = season_clim.__mul__(9.80665)
# hindcast_raw = season_clim.load()
# hindcast_raw_w = Weights(hindcast_raw)

# procesamiento similar al preselect
# 1982-1998, 1999-2011
data_1982_1998 = data.sel(time=data.time.dt.year.isin(np.linspace(1982,1998,17)))
data_1999_2011 = data.sel(time=data.time.dt.year.isin(np.linspace(1999,2011,13)))
#----------------------------------------------------------------------------------------------------------------------#
son_clim_82_98, son_clim_99_11, son_anom_82_98, son_anom_99_11 = \
    TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 10)
son_anom_82_98_detrend, son_anom_99_11_detrend = \
    Detrend_Seasons(son_anom_82_98, son_anom_99_11, 10)

son_hindcast_detrend = xr.concat([son_anom_82_98_detrend, son_anom_99_11_detrend], dim='time')
# son_hindcast_detrend_m = son_hindcast_detrend.mean(['time', 'r'])
# son_hindcast_detrend_m_w = Weights(son_hindcast_detrend_m)
#
# son_hindcast_detrend_m = son_hindcast_detrend_m.sel(lat=slice(-90,20)).load()
# son_hindcast_detrend_m_w = son_hindcast_detrend_m_w.sel(lat=slice(-90,20)).load()

son_clim_99_11 = son_clim_99_11.load()


#agregando realtime
#----------------------------------------------------------------------------------------------------------------------#
# Real-time ------------------------------------------------------------------------------------------------------------#
files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                        dir=dir_rt, All=True)
files = sorted(files, key=lambda x: x.split()[0])

#abriendo todos los archivos
data = xr.open_mfdataset(files, decode_times=False) #xr no entiende la codificacion de Leads, r y las fechas
data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
data = data.sel(L=[0.5, 1.5, 2.5, 3.5]) # Solo leads 0 1 2 3
data['L'] = [0,1,2,3]
data = xr.decode_cf(fix_calendar(data)) # corrigiendo fechas
data = data.sel(P=200)
data = data.drop('P')
data = data.drop('Z')

realtime_raw = data.load()


son_realtime_detrend = Anom_Detrend_SeasonRealTime(data, son_clim_99_11, 10)

son_full = xr.concat([son_hindcast_detrend, son_realtime_detrend], dim='time')
son_full = son_full.load()
aux_son = Weights(son_full.mean(['r', 'time']))

#------#
import matplotlib.pyplot as plt
# aux_obs = obs_proc.reindex(lat=list(reversed(obs_proc.lat))).\
#     interp(lon=son_hindcast_detrend_m.lon.values, lat=son_hindcast_detrend_m.lat.values)



out_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
cfsv2 = xr.open_dataset(out_dir + 'hgt_son.nc')
aux = cfsv2.mean(['r', 'time']) - obs_proc.rename({'var':'hgt'})
plt.imshow(aux.hgt, cmap='RdBu_r', vmin=-500, vmax=500);plt.colorbar();plt.show()