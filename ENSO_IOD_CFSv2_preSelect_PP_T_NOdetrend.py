"""
Pre-procesamiento PP
Anomalías respecto a la climatologia del hindcast y detrend de las anomalias
(similar a ENSO_IOD_fixCFSv2_DMI_N34.py)
"""
########################################################################################################################
import xarray as xr
import numpy as np
import pymannkendall as mk
from ENSO_IOD_Funciones import SelectNMMEFiles
#######################################################################################################################
dir_hc = '/pikachu/datos/osman/nmme/monthly/hindcast/'
dir_rt = '/pikachu/datos/osman/nmme/monthly/real_time/'
out_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
v = 'tref'
save_nc = True
testtrend = False
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

# def Detrend_Seasons(season_anom_1982_1998, season_anom_1999_2011, main_month_season):
#
#     for l in [0,1,2,3]:
#         #1982-1998
#         aux_season_anom_1982_1998 \
#             = season_anom_1982_1998.sel(time=season_anom_1982_1998.time.dt.month.isin(main_month_season-l))
#
#         aux = aux_season_anom_1982_1998.mean('r').polyfit(dim='time', deg=1)
#         aux_trend = xr.polyval(aux_season_anom_1982_1998['time'], aux.prec_polyfit_coefficients)
#         if l == 0:
#             season_anom_1982_1998_detrened = aux_season_anom_1982_1998 - aux_trend
#         else:
#             aux_detrend = aux_season_anom_1982_1998 - aux_trend
#             season_anom_1982_1998_detrened = xr.concat([season_anom_1982_1998_detrened, aux_detrend], dim='time')
#
#     # 1999-2011
#         aux_season_anom_1999_2011 \
#             = season_anom_1999_2011.sel(time=season_anom_1999_2011.time.dt.month.isin(main_month_season - l))
#
#         aux = aux_season_anom_1999_2011.mean('r').polyfit(dim='time', deg=1)
#         aux_trend = xr.polyval(aux_season_anom_1999_2011['time'], aux.prec_polyfit_coefficients)
#         if l==0:
#             season_anom_1999_2011_detrend = aux_season_anom_1999_2011 - aux_trend
#         else:
#             aux_detrend = aux_season_anom_1999_2011 - aux_trend
#             season_anom_1999_2011_detrend = xr.concat([season_anom_1999_2011_detrend, aux_detrend], dim='time')
#
#     return season_anom_1982_1998_detrened, season_anom_1999_2011_detrend

def Anom_Detrend_SeasonRealTime(data_realtime, season_clim_1999_2011, main_month_season):

    for l in [0,1,2,3]:
        season_data = data_realtime.sel(time=data_realtime.time.dt.month.isin(main_month_season-l), L=l)
        aux_season_clim_1999_2011 = season_clim_1999_2011.sel(L=l)

        #Anomalia
        season_anom = season_data - aux_season_clim_1999_2011

        #Detrend
        # aux = season_anom.mean('r').polyfit(dim='time', deg=1)
        # aux_trend = xr.polyval(season_anom['time'], aux.prec_polyfit_coefficients)

        if l==0:
            #season_anom_detrend = season_anom - aux_trend
            season_anom_f = season_anom
        else:
            #aux_detrend = season_anom - aux_trend
            #season_anom_detrend = xr.concat([season_anom_detrend, aux_detrend], dim='time')
            season_anom_f = xr.concat([season_anom_f, season_anom], dim='time')

    #return season_anom_detrend
    return season_anom_f

def MKTest(x):

    return mk.original_test(x)[2]

def TrendTest(x, dim='time'):
    return xr.apply_ufunc(
        MKTest, x,
        input_core_dims=[[dim]],
        vectorize=True,
        output_dtypes=[float]
    )

def TestTrendMK(data_dask, main_month_season):

    print('loading data_dask...')
    data_ld = data_dask.load()
    print('done')
    for l in [0,1,2,3]:
        aux_sel = data_ld.sel(time=data_ld.time.dt.month.isin(main_month_season-l))
        ensamble = aux_sel.mean(['r'])
        print('testing ensemble at leadtime ' + str(l))
        aux = TrendTest(ensamble)

        if l==0:
            aux_f = aux
        else:
            aux_f = xr.concat([aux_f, aux], dim='L')

    return aux_f

def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, step=1,
         name_fig='fig', title='title'):

    import cartopy.feature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    comp_var = comp['var']
    fig = plt.figure(figsize=(5, 6), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([270,330, -60,20], crs_latlon)

    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                     levels=levels, transform=crs_latlon, colors=cmap)
    # cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    # cb.ax.tick_params(labelsize=8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.LAND, facecolor='white')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=7)

    plt.title(title, fontsize=10)
    plt.tight_layout()

    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()
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
data = data.sel(lon=slice(275, 330), lat=slice(-60, 15))
data = data.sel(r=slice(1,24))

#media movil de 3 meses para separar en estaciones
data = data.rolling(time=3, center=True).mean()

# 1982-1998, 1999-2011
data_1982_1998 = data.sel(time=data.time.dt.year.isin(np.linspace(1982,1998,17)))
data_1999_2011 = data.sel(time=data.time.dt.year.isin(np.linspace(1999,2011,13)))

#- Climatologias y anomalias detrend por seasons ----------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
# jja_clim_82_98, jja_clim_99_11, jja_anom_82_98, jja_anom_99_11 = \
#     TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 7)
# jja_anom_82_98_detrend, jja_anom_99_11_detrend = \
#     Detrend_Seasons(jja_anom_82_98, jja_anom_99_11, 7)

#jja_hindcast_detrend = xr.concat([jja_anom_82_98_detrend, jja_anom_99_11_detrend], dim='time')
# jja_hindcast = xr.concat([jja_anom_82_98, jja_anom_99_11], dim='time')
#
# jja_clim_99_11 = jja_clim_99_11.load()

# #----------------------------------------------------------------------------------------------------------------------#
# jas_clim_82_98, jas_clim_99_11, jas_anom_82_98, jas_anom_99_11 = \
#     TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 8)
# jas_anom_82_98_detrend, jas_anom_99_11_detrend = \
#     Detrend_Seasons(jas_anom_82_98, jas_anom_99_11, 8)

# jas_hindcast_detrend = xr.concat([jas_anom_82_98_detrend, jas_anom_99_11_detrend], dim='time')
# jas_hindcast = xr.concat([jas_anom_82_98, jas_anom_99_11], dim='time')
#
# jas_clim_99_11 = jas_clim_99_11.load()
#----------------------------------------------------------------------------------------------------------------------#
# aso_clim_82_98, aso_clim_99_11, aso_anom_82_98, aso_anom_99_11 = \
#     TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 9)
# aso_anom_82_98_detrend, aso_anom_99_11_detrend = \
#     Detrend_Seasons(aso_anom_82_98, aso_anom_99_11, 9)

#aso_hindcast_detrend = xr.concat([aso_anom_82_98_detrend, aso_anom_99_11_detrend], dim='time')
# aso_hindcast = xr.concat([aso_anom_82_98, aso_anom_99_11], dim='time')
#
# aso_clim_99_11 = aso_clim_99_11.load()
#----------------------------------------------------------------------------------------------------------------------#
son_clim_82_98, son_clim_99_11, son_anom_82_98, son_anom_99_11 = \
    TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 10)
# son_anom_82_98_detrend, son_anom_99_11_detrend = \
#     Detrend_Seasons(son_anom_82_98, son_anom_99_11, 10)

#son_hindcast_detrend = xr.concat([son_anom_82_98_detrend, son_anom_99_11_detrend], dim='time')
son_hindcast = xr.concat([son_anom_82_98, son_anom_99_11], dim='time')

son_clim_99_11 = son_clim_99_11.load()

#----------------------------------------------------------------------------------------------------------------------#
# Real-time ------------------------------------------------------------------------------------------------------------#
files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                        dir=dir_rt, All=True)
files = sorted(files, key=lambda x: x.split()[0])

if v=='prec':
    # abriendo todos los archivos
    data = xr.open_mfdataset(files, decode_times=False)  # xr no entiende la codificacion de Leads, r y las fechas
    data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data = data.sel(L=[0.5, 1.5, 2.5, 3.5])  # Solo leads 0 1 2 3
    data['L'] = [0, 1, 2, 3]
    data = xr.decode_cf(fix_calendar(data))  # corrigiendo fechas
    data = data.sel(time=data.time.dt.year.isin(np.linspace(2011, 2020, 10)))
    data = data.sel(lon=slice(275, 330), lat=slice(-60, 15))
    data = data.sel(r=slice(1, 24))
    # media movil de 3 meses para separar en estaciones
    data = data.rolling(time=3, center=True).mean()
else:
    print(v)
    print('(es una mierda...)')
    files = [x for x in files if "_2022" not in x and '_2021' not in x]

    # para evitar: ValueError: Resulting object does not have monotonic global indexes along dimension
    # en xr.open_mfdataset
    files0 = files[0:252]
    files1 = files[253:len(files)]

    data0 = xr.open_mfdataset(files0, decode_times=False).sel(
        L=[0.5, 1.5, 2.5, 3.5], M=slice(1, 24), X=slice(275, 330), Y=slice(-60, 15))

    data1 = xr.open_mfdataset(files1, decode_times=False).sel(
        L=[0.5, 1.5, 2.5, 3.5], M=slice(1, 24), X=slice(275, 330), Y=slice(-60, 15))

    data0 = data0.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data0['L'] = [0, 1, 2, 3]
    data0 = xr.decode_cf(fix_calendar(data0))  # corrigiendo fechas
    data0 = data0.sel(time=data0.time.dt.year.isin(2011))
    # media movil de 3 meses para separar en estaciones
    data0 = data0.rolling(time=3, center=True).mean()

    data1 = data1.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data1['L'] = [0, 1, 2, 3]
    data1 = xr.decode_cf(fix_calendar(data1))  # corrigiendo fechas
    data1 = data1.sel(time=data1.time.dt.year.isin(np.linspace(2012, 2020, 9)))
    # media movil de 3 meses para separar en estaciones
    data1 = data1.rolling(time=3, center=True).mean()

    data = xr.concat([data0, data1], dim='time')

#- Anomalias detrend por seasons --------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
# jja_realtime_detrend = Anom_Detrend_SeasonRealTime(data, jja_clim_99_11, 7)
# jas_realtime_detrend = Anom_Detrend_SeasonRealTime(data, jas_clim_99_11, 8)
# aso_realtime_detrend = Anom_Detrend_SeasonRealTime(data, aso_clim_99_11, 9)
# son_realtime_detrend = Anom_Detrend_SeasonRealTime(data, son_clim_99_11, 10)

# jja_realtime = Anom_Detrend_SeasonRealTime(data, jja_clim_99_11, 7)
# jas_realtime = Anom_Detrend_SeasonRealTime(data, jas_clim_99_11, 8)
# aso_realtime = Anom_Detrend_SeasonRealTime(data, aso_clim_99_11, 9)
son_realtime = Anom_Detrend_SeasonRealTime(data, son_clim_99_11, 10)

########################################################################################################################
# jja_f = xr.concat([jja_hindcast, jja_realtime], dim='time')
# jas_f = xr.concat([jas_hindcast, jas_realtime], dim='time')
# aso_f = xr.concat([aso_hindcast, aso_realtime], dim='time')
son_f = xr.concat([son_hindcast, son_realtime], dim='time')

# save ----------------------------------------------------------------------------------------------------------------#
if save_nc:
    # jja_f.to_netcdf(out_dir + v + '_jja_nodetrend.nc')
    # jas_f.to_netcdf(out_dir + v + '_jas_nodetrend.nc')
    # aso_f.to_netcdf(out_dir + v + '_aso_nodetrend.nc')
    son_f.to_netcdf(out_dir + v + '_son_nodetrend.nc')

########################################################################################################################

if testtrend:
    out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Composites/PP/TrendTest/'

    # jja_trend_test = TestTrendMK(data_dask=jja_f, main_month_season=7)
    # jas_trend_test = TestTrendMK(data_dask=jas_f, main_month_season=8)
    # aso_trend_test = TestTrendMK(data_dask=aso_f, main_month_season=9)
    son_trend_test = TestTrendMK(data_dask=son_f, main_month_season=10)

    # jja_trend_test = jja_trend_test.rename({v: 'var'})
    # jas_trend_test = jas_trend_test.rename({v: 'var'})
    # aso_trend_test = aso_trend_test.rename({v: 'var'})
    son_trend_test = son_trend_test.rename({v: 'var'})


    for l in [0, 1, 2, 3]:
        # Plot(comp=jja_trend_test.sel(L=l), levels=[0, 0.05], cmap='firebrick', dpi=200, save=True,
        #      title=v + ' JJA - Significant Trend - Leadtime ' + str(l), name_fig=v + 'jja_trend_l_' + str(l))
        #
        # Plot(comp=jas_trend_test.sel(L=l), levels=[0, 0.05], cmap='firebrick', dpi=200, save=True,
        #      title=v + ' JAS - Significant Trend - Leadtime ' + str(l), name_fig=v + 'jas_trend_l_' + str(l))
        #
        # Plot(comp=aso_trend_test.sel(L=l), levels=[0, 0.05], cmap='firebrick', dpi=200, save=True,
        #      title=v + ' ASO - Significant Trend - Leadtime ' + str(l), name_fig=v + 'aso_trend_l_' + str(l))

        Plot(comp=son_trend_test.sel(L=l), levels=[0, 0.05], cmap='firebrick', dpi=200, save=True,
             title=v + ' SON - Significant Trend - Leadtime ' + str(l), name_fig=v + 'son_trend_l_' + str(l))

