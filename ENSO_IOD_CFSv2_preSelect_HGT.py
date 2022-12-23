"""
Pre-procesamiento HGT200
Anomalías respecto a la climatologia del hindcast y detrend de las anomalias
(similar a ENSO_IOD_fixCFSv2_DMI_N34.py)
"""
########################################################################################################################
import xarray as xr
import numpy as np
from ENSO_IOD_Funciones import SelectNMMEFiles
#######################################################################################################################
dir_hc = '/pikachu/datos/luciano.andrian/hindcast/'
dir_rt = '/pikachu/datos/luciano.andrian/real_time/'
out_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
v = 'hgt'
# Funciones ############################################################################################################
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

# usando SelectNMMEFiles con All=True,
# abre TODOS los archivos .nc de la ruta en dir

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
data = data.sel(lon=slice(30, 340), lat=slice(-80, 20), P=200)
data = data.drop('P')
#data = data.drop('Z')

#media movil de 3 meses para separar en estaciones
data = data.rolling(time=3, center=True).mean()

# 1982-1998, 1999-2011
data_1982_1998 = data.sel(time=data.time.dt.year.isin(np.linspace(1982,1998,17)))
data_1999_2011 = data.sel(time=data.time.dt.year.isin(np.linspace(1999,2011,13)))

#
# #- Climatologias y anomalias detrend por seasons ----------------------------------------------------------------------#
# #----------------------------------------------------------------------------------------------------------------------#
# jja_clim_82_98, jja_clim_99_11, jja_anom_82_98, jja_anom_99_11 = \
#     TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 7)
# jja_anom_82_98_detrend, jja_anom_99_11_detrend = \
#     Detrend_Seasons(jja_anom_82_98, jja_anom_99_11, 7)
#
# jja_hindcast_detrend = xr.concat([jja_anom_82_98_detrend, jja_anom_99_11_detrend], dim='time')
#
# jja_clim_99_11 = jja_clim_99_11.load()
#
# #----------------------------------------------------------------------------------------------------------------------#
# jas_clim_82_98, jas_clim_99_11, jas_anom_82_98, jas_anom_99_11 = \
#     TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 8)
# jas_anom_82_98_detrend, jas_anom_99_11_detrend = \
#     Detrend_Seasons(jas_anom_82_98, jas_anom_99_11, 8)
#
# jas_hindcast_detrend = xr.concat([jas_anom_82_98_detrend, jas_anom_99_11_detrend], dim='time')
#
# jas_clim_99_11 = jas_clim_99_11.load()
# #----------------------------------------------------------------------------------------------------------------------#
# aso_clim_82_98, aso_clim_99_11, aso_anom_82_98, aso_anom_99_11 = \
#     TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 9)
# aso_anom_82_98_detrend, aso_anom_99_11_detrend = \
#     Detrend_Seasons(aso_anom_82_98, aso_anom_99_11, 9)
#
# aso_hindcast_detrend = xr.concat([aso_anom_82_98_detrend, aso_anom_99_11_detrend], dim='time')
#
# aso_clim_99_11 = aso_clim_99_11.load()
#----------------------------------------------------------------------------------------------------------------------#
son_clim_82_98, son_clim_99_11, son_anom_82_98, son_anom_99_11 = \
    TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 10)
son_anom_82_98_detrend, son_anom_99_11_detrend = \
    Detrend_Seasons(son_anom_82_98, son_anom_99_11, 10)

son_hindcast_detrend = xr.concat([son_anom_82_98_detrend, son_anom_99_11_detrend], dim='time')

son_clim_99_11 = son_clim_99_11.load()
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
data = data.sel(lon=slice(30, 340), lat=slice(-80, 20), P=200)
data = data.drop('P')
data = data.drop('Z')

#media movil de 3 meses para separar en estaciones
data = data.rolling(time=3, center=True).mean()

#- Anomalias detrend por seasons --------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
# jja_realtime_detrend = Anom_Detrend_SeasonRealTime(data, jja_clim_99_11, 7)
# jas_realtime_detrend = Anom_Detrend_SeasonRealTime(data, jas_clim_99_11, 8)
# aso_realtime_detrend = Anom_Detrend_SeasonRealTime(data, aso_clim_99_11, 9)
son_realtime_detrend = Anom_Detrend_SeasonRealTime(data, son_clim_99_11, 10)

########################################################################################################################
# jja_f = xr.concat([jja_hindcast_detrend, jja_realtime_detrend], dim='time')
# jas_f = xr.concat([jas_hindcast_detrend, jas_realtime_detrend], dim='time')
# aso_f = xr.concat([aso_hindcast_detrend, aso_realtime_detrend], dim='time')
son_f = xr.concat([son_hindcast_detrend, son_realtime_detrend], dim='time')

# save ----------------------------------------------------------------------------------------------------------------#
# jja_f.to_netcdf(out_dir + 'hgt_jja.nc')
# jas_f.to_netcdf(out_dir + 'hgt_jas.nc')
# aso_f.to_netcdf(out_dir + 'hgt_aso.nc')
son_f.to_netcdf(out_dir + 'hgt_son.nc')
########################################################################################################################