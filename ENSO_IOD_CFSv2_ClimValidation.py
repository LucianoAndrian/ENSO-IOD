"""
Validacion climatologica CFSv2 SON HGT200 contra ERA5
Periodo: 1981-2020
Con y sin tendencia
"""
################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from ENSO_IOD_Funciones import SelectNMMEFiles
import warnings
warnings.filterwarnings("ignore")
################################################################################
dir_hc = '/pikachu/datos/luciano.andrian/hindcast/'
dir_rt = '/pikachu/datos/luciano.andrian/real_time/'
out_dir = '/pikachu/datos/luciano.andrian/val_clim_cfsv2/'
v = 'hgt'
#------------------------------------------------------------------------------#
compute = True
save = True
dpi = 300
################################################################################
def fix_calendar(ds, timevar='time'):
    """
    agrega los dias a los archivos nc de NMME
    """
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds

def Plot(comp, comp_var, levels, save, dpi, title, name_fig, out_dir,
         color_map, cmap):

    import matplotlib.pyplot as plt
    import cartopy.feature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import cartopy.crs as ccrs
    fig_size = (9, 3.5)
    extent = [0, 359, -80, 20]
    xticks = np.arange(0, 360, 30)
    yticks = np.arange(-80, 20, 10)
    crs_latlon = ccrs.PlateCarree()

    levels_contour = levels.copy()
    if isinstance(levels, np.ndarray):
        levels_contour = levels[levels != 0]
    else:
        levels_contour.remove(0)


    fig = plt.figure(figsize=fig_size, dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_extent(extent, crs=crs_latlon)
    im = ax.contourf(comp.lon, comp.lat, comp_var, levels=levels,
                     transform=crs_latlon, cmap=cmap, extend='both')

    cb = plt.colorbar(im, fraction=0.042, pad=0.035, shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='white', edgecolor=color_map)
    # ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color=color_map, linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(xticks, crs=crs_latlon)
    ax.set_yticks(yticks, crs=crs_latlon)
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

def Weights(data):
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180),
                                   (len(data.lon), 1)))
    data_w = data * weights
    return data_w
################################################################################
if compute:
    print('compute = True')
    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    print('procesando hindcast...')

    files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                            dir=dir_hc, All=True)
    files = sorted(files, key=lambda x: x.split()[0])
    data = xr.open_mfdataset(files, decode_times=False)
    data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data = data.sel(L=[0.5, 1.5, 2.5, 3.5])  # Solo leads 0 1 2 3
    data['L'] = [0, 1, 2, 3]
    data = xr.decode_cf(fix_calendar(data))  # corrigiendo fechas
    data = data.sel(lat=slice(-80, 20))

    # media movil de 3 meses para separar en estaciones
    hindcast = data.rolling(time=3, center=True).mean()
    hindcast = hindcast.load()

    # media climatologica hindcast sin filtrar tendencia
    hindcast2 = hindcast.mean(['r', 'time', 'L']).sel(P=200)
    #hindcast2.to_netcdf(out_dir + 'hindcast_cfsv2_meanclim_son.nc')
    del data

    # Detrend -----------------------------------------------------------------#
    print('Detrend...')
    # dos climatologias
    data_1982_1998 = hindcast.sel(
        time=hindcast.time.dt.year.isin(np.linspace(1982, 1998, 17)))
    data_1999_2011 = hindcast.sel(
        time=hindcast.time.dt.year.isin(np.linspace(1999, 2011, 13)))

    for l in [0, 1, 2, 3]:
        # 1982_1998
        season_1982_1998 = data_1982_1998.sel(
            time=data_1982_1998.time.dt.month.isin(10 - l), L=l)
        #tendencia
        aux = season_1982_1998.mean('r').polyfit(dim='time', deg=1)
        aux_trend = xr.polyval(season_1982_1998['time'],
                              aux.hgt_polyfit_coefficients)

        if l == 0:
            season_1982_1998_detrened = season_1982_1998 - aux_trend
        else:
            aux_detrend = season_1982_1998 - aux_trend

            season_1982_1998_detrened = \
                xr.concat([season_1982_1998_detrened, aux_detrend],
                          dim='time')

        # 1999-2011
        season_1999_2011 = data_1999_2011.sel(
            time=data_1999_2011.time.dt.month.isin(10 - l), L=l)
        #tendencia
        aux = season_1999_2011.mean('r').polyfit(dim='time', deg=1)
        aux_trend = xr.polyval(season_1999_2011['time'],
                           aux.hgt_polyfit_coefficients)

        if l == 0:
            season_1999_2011_detrend = season_1999_2011 - aux_trend
            season_clim_1999_2011 = season_1999_2011.mean(['r', 'time'])
        else:
            aux_detrend = season_1999_2011 - aux_trend

            season_anom_1999_2011_detrend = \
                xr.concat([season_1999_2011_detrend,
                           aux_detrend], dim='time')

            season_clim_1999_2011 = xr.concat(
                [season_clim_1999_2011,
                 season_1999_2011.mean(['r', 'time'])], dim='L')

    son_hindcast_detrend = xr.concat(
        [season_1982_1998_detrened,
         season_1999_2011_detrend], dim='time')
    son_hindcast_detrend_mean = \
        son_hindcast_detrend.mean(['r', 'time']).sel(P=200)

    # son_hindcast_detrend.to_netcdf(
    #     out_dir + 'hindcast_cfsv2_meanclim_detrend_son.nc')

    # -------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    print('procesando realtime...')
    files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                            dir=dir_rt, All=True)
    files = sorted(files, key=lambda x: x.split()[0])

    data = xr.open_mfdataset(files,
                             decode_times=False)
    data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data = data.sel(L=[0.5, 1.5, 2.5, 3.5])  # Solo leads 0 1 2 3
    data['L'] = [0, 1, 2, 3]
    data = xr.decode_cf(fix_calendar(data))  # corrigiendo fechas
    data = data.sel(lat=slice(-80, 20), P=200)
    data = data.drop('P')
    data = data.drop('Z')

    # media movil de 3 meses para separar en estaciones
    real_time = data.rolling(time=3, center=True).mean()
    real_time = real_time.load()
    real_time2 = real_time.mean(['r', 'time', 'L'])
    #real_time2.to_netcdf(out_dir + 'real_time_cfsv2_meanclim_son.nc')
    del data

    # Detrend -----------------------------------------------------------------#
    print('Detrend...')
    # la climatologia usada es la de hindcast 1998-2011
    for l in [0,1,2,3]:
        season_data = real_time.sel(
            time=real_time.time.dt.month.isin(10-l), L=l)
        aux_clim_1999_2011 = season_clim_1999_2011.sel(L=l)
        #Anomalia
        season_anom = season_data - aux_clim_1999_2011
        #Detrend
        aux = season_anom.mean('r').polyfit(dim='time', deg=1)
        aux_trend = xr.polyval(season_anom['time'],
                               aux.hgt_polyfit_coefficients)

        if l==0:
            son_realtime_detrend = season_anom - aux_trend
        else:
            aux_detrend = season_anom - aux_trend
            son_realtime_detrend = \
                xr.concat([son_realtime_detrend, aux_detrend], dim='time')

    son_realtime_detrend_mean = \
        son_realtime_detrend.mean(['r', 'time']).sel(P=200)
    # son_realtime_detrend.to_netcdf(
    #     out_dir + 'realtime_cfsv2_meanclim_detrend_son.nc')

else:
    print('compute = False')
    print('Opening CFSv2 saved files')
    hindcast = xr.open_dataset(out_dir + 'hindcast_cfsv2_meanclim_son.nc')
    son_hindcast_detrend = \
            xr.open_dataset(out_dir + 'hindcast_cfsv2_meanclim_detrend_son.nc')

    real_time = xr.open_dataset(out_dir + 'real_time_cfsv2_meanclim_son.nc')
    son_realtime_detrend = \
            xr.open_dataset(out_dir + 'realtime_cfsv2_meanclim_detrend_son.nc')


###############################################################################
# ERA
print('ERA5...')
data_dir_era5 = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/' \
                        'downloaded/'
data = xr.open_dataset(data_dir_era5 + 'ERA5_HGT200_40-20.nc')
data = data.sel(time=data.time.dt.year.isin(range(1981, 2012)))
data = data.rolling(time=3, center=True).mean()
data = data.sel(time=data.time.dt.month.isin(10)) #son
data = data.rename({'z':'hgt', 'longitude':'lon', 'latitude':'lat'})
data = data.sel(time=data.time.dt.year.isin(range(1982, 2020)))
# era_hind = data.sel(time=data.time.dt.year.isin(range(1981, 2012)))
# era_real = data.sel(time=data.time.dt.year.isin(range(2011, 2020)))

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# era_hind = era_hind.interp(lat=np.arange(-90,91)[::-1], lon=np.arange(0,361))
# era_hind = era_hind.mean('time').__mul__(1/9.8)
# era_real = era_real.interp(lat=np.arange(-90,91)[::-1], lon=np.arange(0,361))
# era_real = era_real.mean('time').__mul__(1/9.8)
era_clim = data.interp(lat=np.arange(-80,21)[::-1], lon=np.arange(0,360))
era_clim = era_clim.mean('time').__mul__(1/9.8)


# sin tendencia
era_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/1940_2020/'
data = xr.open_dataset(era_dir + 'HGT200_SON_mer_d_w.nc')
data = data.rename({'var':'hgt'})
data = data.sel(time=data.time.dt.year.isin(range(1981, 2021)))

era_full = data.sel(time=data.time.dt.year.isin(range(1982, 2021)))
era_full = era_full.interp(
    lat=np.arange(-80,21)[::-1], lon=np.arange(0,360))
era_full = era_full.__mul__(1/9.8)
#
# era_hind_d = data.sel(time=data.time.dt.year.isin(range(1981, 2012)))
# era_hind_d = era_hind_d.interp(
#     lat=np.arange(-90,91)[::-1], lon=np.arange(0,361))
# era_hind_d = era_hind_d.mean('time').__mul__(1/9.8)
#
# era_real_d = data.sel(time=data.time.dt.year.isin(range(2011, 2020)))
# era_real_d = era_real_d.interp(
#     lat=np.arange(-90,91)[::-1], lon=np.arange(0,360))
# era_real_d = era_real_d.mean('time').__mul__(1/9.8)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
print('Plot...')
# total sin tendencia
# aux = xr.concat([era_hind_d, era_real_d], dim='time').mean('time') + \
#       Weights(xr.concat([era_hind, era_real], dim='time').mean('time'))

aux_hind = son_hindcast_detrend + hindcast2
aux_hind = Weights(aux_hind.sel(P=200))

aux_real = son_realtime_detrend.sel(P=200) + real_time2
aux_real = Weights(aux_real)

cfsv2 = xr.concat([aux_hind, aux_real], dim='time')
era = era_full + Weights(era_clim)

dif = cfsv2.mean(['time','r']) - era.mean('time')
Plot(dif, dif.hgt, np.arange(-100, 120, 20), save, dpi,
     'CFSv2 hindcast + real time - ERA5 (detrended)', 'dif_hind_real_d.jpg',
     out_dir, 'k', 'RdBu_r')


# Testeo ----------------------------------------------------------------------#
print('Testeo de diferencia de medias con t-student por leadtime')
#from scipy import stats
# # promediando todos los r (y leadtimes, sino hay condicion sobre el tiempo)
# aux = stats.ttest_ind(
#     cfsv2.mean('r').hgt,
#     era.isel(lat=slice(None, None, -1)).hgt,
#     equal_var=False)
# plt.imshow(aux[0]);plt.colorbar();plt.show()

# por leadtimes
from scipy.stats import ttest_ind
pvalue = []
for m in [7, 8, 9, 10]:
    cfsv2_monthly_mean = cfsv2.sel(
        time=cfsv2.time.dt.month.isin(m)).mean('r').hgt
    era_inverted = era.isel(lat=slice(None, None, -1)).hgt
    # test
    pvalue.append(
        ttest_ind(cfsv2_monthly_mean, era_inverted, equal_var=False)[1])

# promedio de pvalue por leadtime
pvalue = sum(pvalue) / len(pvalue)

Plot(dif, dif.where(pvalue<0.05).hgt, np.arange(-100, 120, 20), save, dpi,
     'CFSv2 hindcast + real time - ERA5 (detrended)',
     'dif_hind_real_d_tested.jpg',
     out_dir, 'k', 'RdBu_r')
#------------------------------------------------------------------------------#
print('done')
#------------------------------------------------------------------------------#

# print('plot...')
# # hindcast sin proc
# dif = Weights(hindcast) - Weights(era_hind)
# Plot(dif, dif.hgt, np.arange(-100, 120, 20), save, dpi,
#      'CFSv2 hindcast - ERA5', 'dif_hind.jpg', out_dir, 'k', 'RdBu_r')
#
# # Real time sin proc
# dif = Weights(real_time) - Weights(era_real)
# Plot(dif, dif.hgt, np.arange(-100, 120, 20), save, dpi,
#      'CFSv2 real time - ERA5', 'dif_real.jpg', out_dir, 'k', 'RdBu_r')
#
# # Sin tendencia ---------------------------------------------------------------#
# # hindcast
# dif1 = Weights(son_hindcast_detrend) + Weights(hindcast)
# dif2 = era_hind_d + Weights(era_hind) # era_hind_d ya fue pesado por lat.
# dif = dif1 - dif2
# Plot(dif, dif.hgt, np.arange(-100, 120, 20), save, dpi,
#      'CFSv2 hindcast - ERA5 (detrended)', 'dif_hind_d.jpg', out_dir,
#      'k', 'RdBu_r')
#
# # realtime
# dif1 = Weights(son_realtime_detrend) + Weights(real_time)
# dif2 = era_real_d + Weights(era_real)
# dif = dif1 - dif2
# Plot(dif, dif.hgt, np.arange(-100, 120, 20), save, dpi,
#      'CFSv2 real time - ERA5 (detrended)', 'dif_realtime_d.jpg', out_dir,
#      'k', 'RdBu_r')
#
# #------------------------------------------------------------------------------#
# # total -----------------------------------------------------------------------#
# # total sin proc
# era_full = Weights(
#     xr.concat([era_hind, era_real], dim='time').mean('time'))
# cfsv2_full = Weights(
#     xr.concat([hindcast, real_time], dim='time').mean('time'))
#
# dif = cfsv2_full - era_full
# Plot(dif, dif.hgt, np.arange(-100, 120, 20), save, dpi,
#      'CFSv2 hindcast + real time - ERA5', 'dif_hind_real.jpg', out_dir,
#      'k', 'RdBu_r')
