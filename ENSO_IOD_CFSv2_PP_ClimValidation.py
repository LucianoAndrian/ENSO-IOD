"""
Validacion climatologica CFSv2 SON HGT200 contra ERA5
Periodo: 1981-2020
Con y sin tendencia
"""
################################################################################
compute = True
save = True
################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from ENSO_IOD_Funciones import SelectNMMEFiles
import warnings
warnings.filterwarnings("ignore")
################################################################################
dir_hc = '/pikachu/datos/osman/nmme/monthly/hindcast/'
dir_rt = '/pikachu/datos/osman/nmme/monthly/real_time/'
out_dir = '/pikachu/datos/luciano.andrian/val_clim_cfsv2/'
v = 'prec'
#------------------------------------------------------------------------------#
if save:
    dpi = 300
else:
    dpi = 100
################################################################################
def fix_calendar(ds, timevar='time'):
    """
    agrega los dias a los archivos nc de NMME
    """
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds

def Plot(comp, comp_var, levels, save, dpi, title, name_fig, out_dir, cmap):

    import matplotlib.pyplot as plt
    import cartopy.feature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import cartopy.crs as ccrs

    fig_size = (5, 6)
    extent = [270, 330, -60, 20]
    xticks = np.arange(275, 330+1, 10)
    yticks = np.arange(-60, 15+1, 10)
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

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.BORDERS, facecolor='k')
    ax.add_feature(cartopy.feature.OCEAN, zorder=10, facecolor='white',
                   edgecolor='k')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-', zorder=12)
    ax.set_xticks(xticks, crs=crs_latlon)
    ax.set_yticks(yticks, crs=crs_latlon)
    # al usar ocean como mascara
    # y las girdlanes encima
    # los bordes del plot quedan tapadaos por las gridlines
    for k, spine in ax.spines.items():
        spine.set_zorder(13)

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
    data = data.sel(lat=slice(-70, 20))
    data = data.sel(lon=slice(275, 331))

    # media movil de 3 meses para separar en estaciones
    hindcast = data.rolling(time=3, center=True).mean()
    hindcast = hindcast.load()

    # media climatologica hindcast sin filtrar tendencia
    hindcast2 = (hindcast.mean(['r', 'time', 'L'])/
                 hindcast.std(['r', 'time', 'L']))
    hindcast3 = hindcast.mean(['r', 'time', 'L'])
    hindcast2.to_netcdf(f"{out_dir}hindcast_{v}_cfsv2_meanclim_son.nc")
    hindcast3.to_netcdf(f"{out_dir}hindcast_no-norm_{v}"
                        f"_cfsv2_meanclim_son.nc")
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
                              aux.prec_polyfit_coefficients)

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
                           aux.prec_polyfit_coefficients)

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
        son_hindcast_detrend.mean(['r', 'time'])

    son_hindcast_detrend.to_netcdf(f"{out_dir}hindcast_{v}_cfsv2_meanclim_"
                                   f"detrend_son.nc")

    # -------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    print('procesando realtime...')
    files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                            dir=dir_rt, All=True)
    files = sorted(files, key=lambda x: x.split()[0])
    # files = [x for x in files if "_2022" not in x and '_2021' not in x]
    #
    # # para evitar: ValueError:
    # # Resulting object does not have monotonic global indexes along dimension
    # # en xr.open_mfdataset
    # files0 = files[0:252]
    # files1 = files[253:len(files)]

    data = xr.open_mfdataset(files,
                             decode_times=False)
    data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data = data.sel(L=[0.5, 1.5, 2.5, 3.5])  # Solo leads 0 1 2 3
    data['L'] = [0, 1, 2, 3]
    data = xr.decode_cf(fix_calendar(data))  # corrigiendo fechas
    data = data.sel(lat=slice(-70, 20))
    data = data.sel(lon=slice(275, 331))

    # media movil de 3 meses para separar en estaciones
    real_time = data.rolling(time=3, center=True).mean()
    real_time = real_time.load()
    real_time2 = (real_time.mean(['r', 'time', 'L'])/
                  real_time.std(['r', 'time', 'L']))
    real_time3 = real_time.mean(['r', 'time', 'L'])
    real_time2.to_netcdf(f"{out_dir}real_time_{v}_cfsv2_meanclim_son.nc")
    real_time3.to_netcdf(
        f"{out_dir}real_time_no-norm_{v}_cfsv2_meanclim_son.nc")
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
                               aux.prec_polyfit_coefficients)

        if l==0:
            son_realtime_detrend = season_anom - aux_trend
        else:
            aux_detrend = season_anom - aux_trend
            son_realtime_detrend = \
                xr.concat([son_realtime_detrend, aux_detrend], dim='time')

    son_realtime_detrend_mean = \
        son_realtime_detrend.mean(['r', 'time'])
    son_realtime_detrend.to_netcdf(f"{out_dir}realtime_{v}_cfsv2_meanclim_"
                                   f"detrend_son.nc")


print('Opening CFSv2 saved files')
hindcast = xr.open_dataset(f"{out_dir}hindcast_{v}_cfsv2_meanclim_son.nc")
hindcast_no_norm = (
    xr.open_dataset(f"{out_dir}hindcast_no-norm_{v}_cfsv2_meanclim_son.nc"))
son_hindcast_detrend = (
    xr.open_dataset(f"{out_dir}hindcast_{v}_cfsv2_meanclim_detrend_son.nc"))

real_time = xr.open_dataset(f"{out_dir}real_time_{v}_cfsv2_meanclim_son.nc")
real_time_no_norm = (
    xr.open_dataset(f"{out_dir}real_time_no-norm_{v}_cfsv2_meanclim_son.nc"))
son_realtime_detrend = (
    xr.open_dataset(f"{out_dir}realtime_{v}_cfsv2_meanclim_detrend_son.nc"))


###############################################################################
print('GPCC...')
data_dir_pp_clim = ('/pikachu/datos/luciano.andrian/observado/ncfiles/'
                    'data_no_detrend/')
data = xr.open_dataset(data_dir_pp_clim + 'pp_pgcc_v2020_1891-2023_1.nc')
data = data.sel(lat=slice(20, -80), lon=slice(275, 331),
                time=data.time.dt.year.isin(range(1982, 2020)))
data = data.interp(lat=hindcast.lat.values, lon=hindcast.lon.values)
data = data.rolling(time=3, center=True).mean()
data = data.sel(time=data.time.dt.month.isin(10)) #son
data = data.rename({'precip':'prec'})
pp_clim = data.mean('time')
del data

# sin tendencia
data_dir_pp = ('/pikachu/datos/luciano.andrian/observado/ncfiles/'
               'data_obs_d_w_c/')
data = xr.open_dataset(data_dir_pp + 'ppgpcc_w_c_d_1_SON.nc')
data = data.sel(time=data.time.dt.year.isin(range(1981, 2021)),
                lon=slice(275,331))
data = data.rename({'var':'prec'})
data = data.interp(lat=hindcast.lat.values, lon=hindcast.lon.values)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
print('Plot...')
# total sin tendencia
aux_hind = son_hindcast_detrend + hindcast
aux_real = son_realtime_detrend + real_time
cfsv2 = xr.concat([aux_hind, aux_real], dim='time')

aux_hind_no_norm = son_hindcast_detrend + hindcast_no_norm
aux_real_no_norm = son_realtime_detrend + real_time_no_norm
cfsv2_no_norm = xr.concat([aux_hind_no_norm,
                           aux_real_no_norm], dim='time')*30

pp = data + pp_clim
#------------------------------------------------------------------------------#
# VEEER normalización
pp_norm = pp/pp.std('time')
dif = cfsv2.mean(['time','r']) - pp_norm.mean('time')
dif = Weights(dif)

dif_no_norm = cfsv2_no_norm.mean(['time','r']) - pp.mean('time')
dif_no_norm = Weights(dif_no_norm)

from matplotlib import colors
cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC',
                                 '#B4E2DB',
                                 'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07',
                                 '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')
scale_pp = [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2]
scale_pp_no_norm = np.linspace(-60,60,13)

# Plot(dif, dif.prec, scale_pp, save, dpi,
#      'CFSv2 hindcast + real time - GPCC (Norm y detr)',
#      'prec_dif_hind_real_d.jpg', out_dir, cbar_pp)


# Plot(dif_no_norm, dif_no_norm.prec, scale_pp_no_norm, save, dpi,
#      'CFSv2 hindcast + real time - GPCC (No-Norm y detr)',
#      'prec_no-norm_dif_hind_real_d.jpg', out_dir, cbar_pp)


# Testeo ----------------------------------------------------------------------#
print('Testeo de diferencia de medias con t-student por leadtime')

# pvalue = []
# for m in [7, 8, 9, 10]:
#     cfsv2_monthly_mean = cfsv2.sel(
#         time=cfsv2.time.dt.month.isin(m)).mean('r').prec
#     pp_prec = pp_norm.prec#.isel(lat=slice(None, None, -1)).prec
#     # test
#     pvalue.append(
#         ttest_ind(cfsv2_monthly_mean, pp_prec, equal_var=False)[1])
#
# # promedio de pvalue por leadtime
# pvalue = sum(pvalue) / len(pvalue)
#
# Plot(dif, dif.where(pvalue<0.05).prec, scale_pp, save, dpi,
#      'CFSv2 hindcast + real time - GPCC (Norm y detr)',
#      'prec_dif_hind_real_d_tested.jpg', out_dir, cbar_pp)


pvalue = []
for m in [7, 8, 9, 10]:
    cfsv2_monthly_mean = cfsv2_no_norm.sel(
        time=cfsv2_no_norm.time.dt.month.isin(m)).mean('r').prec
    pp_prec = pp.prec#.isel(lat=slice(None, None, -1)).prec
    # test
    pvalue.append(
        ttest_ind(cfsv2_monthly_mean, pp_prec, equal_var=False)[1])

# promedio de pvalue por leadtime
pvalue = sum(pvalue) / len(pvalue)

Plot(dif_no_norm, dif_no_norm.where(pvalue<0.05).prec, scale_pp_no_norm,
     save, dpi,
     'CFSv2 hindcast + real time - GPCC (No-Norm y detr)',
     'prec_no-norm_dif_hind_real_d_tested.jpg', out_dir, cbar_pp)

#------------------------------------------------------------------------------#
print('done')
#------------------------------------------------------------------------------#
