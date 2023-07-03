"""
ENSO vs IOD Regression
T, PP y HGT200
"""
################################################################################
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
import xarray as xr
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import os
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import colors
from ENSO_IOD_Funciones import Nino34CPC
from ENSO_IOD_Funciones import DMI
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

from ENSO_IOD_Funciones import ComputeWithEffect, ComputeWithoutEffect, WAF
################################################################################
#era5_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'
era5_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/1940_2020/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/paper1/1940_2020/' \
          'regression/'
################################################################################
save = True
dpi = 300
full_season = False
text = False
waf = True
# Functions ####################################################################


def PlotReg(data, data_cor, levels=np.linspace(-100,100,2), cmap='RdBu_r'
            , dpi=100, save=False, title='\m/', name_fig='fig_PlotReg', sig=True
            ,two_variables = False, data2=None, data_cor2=None,
            levels2 = np.linspace(-100,100,2), sig2=True, step=1,SA=False,
            color_map = '#d9d9d9', color_sig='magenta', sig_point=False,
            r_crit=1, waf=False, data_waf=None, px=None, py=False,
            waf_scale=1 / 1000, step_waf=10):


    from numpy import ma

    levels_contour = levels.copy()
    if isinstance(levels_contour, np.ndarray):
        levels_contour = levels_contour[levels_contour != 0]
    else:
        levels_contour.remove(0)

    crs_latlon = ccrs.PlateCarree()
    if SA:
        fig = plt.figure(figsize=(5, 6), dpi=dpi)
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
        ax.set_extent([270,330, -60,20], crs=crs_latlon)
    else:
        fig = plt.figure(figsize=(9, 3.5), dpi=dpi)
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
        ax.set_extent([0, 359, -80, 20], crs=crs_latlon)

    ax.contour(data.lon[::step], data.lat[::step], data[::step, ::step],
               linewidths=.5, alpha=0.5,
               levels=levels_contour, transform=crs_latlon, colors='black')

    im = ax.contourf(data.lon[::step], data.lat[::step], data[::step,::step],
                     levels=levels,
                     transform=crs_latlon, cmap=cmap, extend='both')
    if sig:
        if sig_point:
            colors_l = [color_sig, color_sig]
            cs = ax.contourf(data_cor.lon, data_cor.lat,
                             data_cor.where(np.abs(data_cor) > np.abs(r_crit)),
                             transform=crs_latlon, colors='none',
                             hatches=["...", "..."],
                             extend='lower')
            for i, collection in enumerate(cs.collections):
                collection.set_edgecolor(colors_l[i % len(colors_l)])

            for collection in cs.collections:
                collection.set_linewidth(0.)
            # para hgt200 queda mejor los dos juntos
            ax.contour(data_cor.lon[::step], data_cor.lat[::step],
                       data_cor[::step, ::step],
                       levels=np.linspace(-r_crit, r_crit, 2),
                       colors=color_sig, transform=crs_latlon, linewidths=1)

        else:
            ax.contour(data_cor.lon[::step], data_cor.lat[::step],
                       data_cor[::step, ::step],
                       levels=np.linspace(-r_crit, r_crit, 2),
                       colors=color_sig, transform=crs_latlon, linewidths=1)


    if two_variables:
        ax.contour(data2.lon, data2.lat, data2, levels=levels2,
                   colors='k', transform=crs_latlon, linewidths=1)
        if sig2:
            ax.contour(data_cor2.lon, data_cor2.lat, data_cor2,
                       levels=np.linspace(-r_crit, r_crit, 2),
                       colors='forestgreen', transform=crs_latlon, linewidths=1)

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='white', edgecolor=color_map)
    #ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    # ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=
    
    if waf:
        from numpy import ma
        Q60 = np.nanpercentile(np.sqrt(np.add(np.power(px, 2),
                                              np.power(py, 2))), 60)
        M = np.sqrt(np.add(np.power(px, 2), np.power(py, 2))) < Q60
        # mask array
        px_mask = ma.array(px, mask=M)
        py_mask = ma.array(py, mask=M)

        Q99 = np.nanpercentile(np.sqrt(np.add(np.power(px, 2),
                                              np.power(py, 2))), 99)
        M = np.sqrt(np.add(np.power(px, 2), np.power(py, 2))) > Q99
        # mask array
        px_mask = ma.array(px_mask, mask=M)
        py_mask = ma.array(py_mask, mask=M)

        # plot vectors
        lons, lats = np.meshgrid(data_waf.lon.values, data_waf.lat.values)
        ax.quiver(lons[::step_waf, ::step_waf], lats[::step_waf, ::step_waf],
                  px_mask[0, ::step_waf, ::step_waf],
                  py_mask[0, ::step_waf, ::step_waf],
                  transform=crs_latlon, pivot='tail', width=1.5e-3,
                  headwidth=3, alpha=1,
                  headlength=2.5, color='k', scale=waf_scale)
    
    if SA:
        ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5, zorder=17)
        #ax.add_feature(cartopy.feature.COASTLINE)
        ocean = cartopy.feature.NaturalEarthFeature('physical', 'ocean',
                                                    scale='50m',
                                                    facecolor='white', alpha=1)
        ax.add_feature(ocean, linewidth=0.2, zorder=15)
        ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
        ax.set_yticks(np.arange(-60, 20, 20), crs=crs_latlon)

        ax2 = ax.twinx()
        ax2.set_yticks([])
        #ax2.set_xticks([])

    else:
        ax.set_xticks(np.arange(0, 360, 30), crs=crs_latlon)
        ax.set_yticks(np.arange(-80, 20, 10), crs=crs_latlon)
        ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
        ax.coastlines(color=color_map, linestyle='-', alpha=1)

    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-', zorder=20)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #ax2.spines['left'].set_color('k')

    ax.tick_params(labelsize=7)
    plt.title(title, fontsize=10)
    plt.tight_layout()

    if save:
        print('save: ' + out_dir + name_fig + '.jpg')
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()

    else:
        plt.show()



def MakerMaskSig(data):
    mask_sig = data.where((data < -1 * r_crit) | (data > r_crit))
    mask_sig = mask_sig.where(np.isnan(mask_sig), 1)

    return mask_sig

def Detrend(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
    dt = xrda - trend
    return dt

def OrdenarNC_wTime_fromW(data):
    newdata = xr.Dataset(
        data_vars=dict(
            var=(['time', 'lat', 'lon'], data['var'][0, :, :, :].values)

        ),
        coords=dict(
            lon=(['lon'], data.lon),
            lat=(['lat'], data.lat),
            time=(['time'], data.time)
        )
    )
    return newdata


def MakeXr(data):
    return xr.Dataset(
        data_vars=dict(
            var=(['lat', 'lon'], data)

        ),
        coords=dict(
            lon=(['lon'], data.lon),
            lat=(['lat'], data.lat)
        )
    )


def ComputeWaf(reg_output, data_clim, hpalevel):
    reg_output = MakeXr(reg_output)
    px, py = WAF(data_clim, reg_output, data_clim.lon, data_clim.lat,
                 reshape=True, variable='var',
                 hpalevel=hpalevel)

    weights = np.transpose(np.tile(-2 * np.cos(data_clim.lat.values * 1 * np.pi / 180) + 2.1, (359, 1)))
    weights_arr = np.zeros_like(px)
    weights_arr[0, :, :] = weights
    px *= weights_arr
    py *= weights_arr

    return px, py
################################################################################
# variables = ['pp_gpcc_or', 't_cru', 'hgt200_HS_mer_d_w', 'pp_gpcc_0.25', 'pp_cmap', 'pp_prec']
# name_var = ['precip']
# title_var = ['PP GPCC', 'T Cru', 'HGT200 ERA5', 'PP GPCC 0.25 pr', 'PP CMAP', 'PP PREC']
# seasons = [7, 10] # main month
# seasons_name = ['JJA', 'SON']
# interp = False
# two_variables=False
# SA = [True, True, False, True, True, True]
sig = True

#variables = ['hgt200_HS_mer_d_w', 'hgt750_mer_d_w']
variables = ['HGT200_SON_mer_d_w', 'HGT750_SON_mer_d_w']
name_var = ['precip']
title_var = ['HGT200 ERA5', 'HGT750 ERA5']
seasons = [10] # main month
seasons_name = ['SON']
interp = False
two_variables=False
SA = [False, False]

scales = [[-150,-100,-75,-50,-25,-15,0,15,25,50,75,100,150],
          [-150,-100,-75,-50,-25,-15,0,15,25,50,75,100,150]]


cbar = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55',
                              '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3',
                              '#2064AF', '#014A9B'][::-1])
cbar.set_over('#641B00')
cbar.set_under('#012A52')
cbar.set_bad(color='white')


cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169',
                                 '#79C8BC', '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13',
                                 '#6A3D07', '#543005'][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')

cmap = [cbar, cbar]
plt.rcParams['hatch.linewidth'] = 0.5




dmi_or = DMI(filter_bwa=False, start_per='1920', end_per='2020')[2]
n34_or = Nino34CPC(xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc"),
                   start=1920, end=2020)[0]

periodos = [[1940,2020]]
t_critic = 1.66 # es MUY similar (2 digitos) para ambos períodos
y1 = 0
p = periodos[0]
r_crit = np.sqrt(1 / (((np.sqrt((p[1] - p[0] + y1) - 2) / t_critic) ** 2) + 1))

for v, v_count, hpalevel, waf_scale in zip(variables,[0,1], [200,750],
                                           [1/800, 1/500]):
    # indices: -----------------------------------------------------------------
    dmi = dmi_or.sel(time=slice(str(p[0] + y1) + '-01-01', str(p[1]) + '-12-01'))
    n34 = n34_or.sel(time=slice(str(p[0] + y1) + '-01-01', str(p[1]) + '-12-01'))

    data = xr.open_dataset(era5_dir + v + '.nc')
    data = data.sel(time=slice(str(p[0] + y1) + '-01-01', str(p[1]) + '-12-01'))
    time_original = data.time

    # Anomaly ------------------------------------------------------------------
    # data = data.groupby('time.month') - \
    #        data.groupby('time.month').mean('time',skipna=True)
    data = data - data.mean('time')

    # 3-month running mean -----------------------------------------------------
    #data = data.rolling(time=3, center=True).mean()

    # --- WAF --- #
    if waf:
        data_sf = xr.open_dataset(era5_dir + 'sf_from_UV' +
                                  str(hpalevel) + '_w.nc')
        data_sf = data_sf.sel(variable='var') # WTF
        #data_sf = Detrend(OrdenarNC_wTime_fromW(
        data_sf = data_sf.rename({'streamfunction': 'var'})

        data_sf = data_sf.sel(time=slice(str(p[0])+'-01-01',str(p[1])+'-12-31'))
        # Climatology SON ------------------------------------------------------
        #aux = data_sf.rolling(time=3, center=True).mean()
        data_clim = data_sf.sel(time=data_sf.time.dt.month.isin(10)).mean('time')
        data_sf = data_sf.sel(time=data_sf.time.dt.month.isin(10))
        # Anomaly --------------------------------------------------------------
        data_sf = data_sf - data_sf.mean('time')
        # data_sf = data_sf.groupby('time.month') -\
        #           data_sf.groupby('time.month').mean('time', skipna=True)
        # 3-month running mean -------------------------------------------------
        # data_sf = data_sf.rolling(time=3, center=True).mean()

    # Seasons ------------------------------------------------------------------
    for s, s_count in zip(seasons_name, np.arange(0,len(seasons_name))):
        aux_n34, aux_corr_n34, aux_dmi, \
        aux_corr_dmi, aux_n34_2, aux_corr_n34_2, \
        aux_dmi_2, aux_corr_dmi_2 = \
            ComputeWithEffect(data=data, data2=None,
                              n34=n34.sel(time=n34.time.dt.month.isin(10)),
                              dmi=dmi.sel(time=dmi.time.dt.month.isin(10)),
                              two_variables=two_variables, m=seasons[s_count],
                              full_season=False, time_original=time_original)

        if waf:
            print('#--- WAF form SF ---#')
            aux_waf_aux_n34, aux_waf_aux_corr_n34, aux_waf_aux_dmi, \
            aux_waf_aux_corr_dmi, aux_waf_aux_n34_2, aux_waf_aux_corr_n34_2, \
            aux_waf_aux_dmi_2, aux_waf_aux_corr_dmi_2 = \
                ComputeWithEffect(data=data_sf, data2=None,
                                  n34=n34.sel(time=n34.time.dt.month.isin(10)),
                                  dmi=dmi.sel(time=dmi.time.dt.month.isin(10)),
                                  two_variables=False,
                                  m=seasons[s_count], full_season=False,
                                  time_original=time_original)
            # N34 WAF
            px_n34, py_n34 = ComputeWaf(aux_waf_aux_n34, data_clim, hpalevel)

            #DMI WAF
            px_dmi, py_dmi = ComputeWaf(aux_waf_aux_dmi, data_clim, hpalevel)
        else:
            px_n34 = None
            py_n34 = None
            px_dmi = None
            py_dmi = None


        print('Plot')
        PlotReg(data=aux_n34 * MakerMaskSig(aux_corr_n34), data_cor=aux_corr_n34,
                levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
                title=title_var[v_count] + '_' + s +'_'+ str(p[0] + y1) + '_' +
                      str(p[1]) + '_Niño3.4',
                name_fig=v + '_' + s + '_' + str(p[0] + y1) + '_' + str(p[1]) +
                         '_N34',
                save=save, sig=False, two_variables=True, data2=aux_n34,
                sig2=False, levels2=scales[v_count], SA=SA[v_count], step=1,
                color_map='grey', color_sig='k', sig_point=True, r_crit=r_crit,
                waf=True, px=px_n34, py=py_n34, data_waf=data_sf,
                waf_scale = 1/1000, step_waf=4)

        PlotReg(data=aux_dmi * MakerMaskSig(aux_corr_dmi), data_cor=aux_corr_dmi,
                levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
                title=title_var[v_count] + '_' + s + '_' + str(p[0] + y1) + '_'+
                      str(p[1]) + '_DMI',
                name_fig=v + '_' + s + '_' + str(p[0] + y1) + '_' + str(p[1]) +
                         '_DMI',
                save=save, sig=False,
                two_variables=True, data2=aux_dmi, sig2=False,
                levels2=scales[v_count],
                SA=SA[v_count], step=1, color_map='grey', color_sig='k',
                sig_point=False, r_crit=r_crit,
                waf=True, px=px_dmi, py=py_dmi, data_waf=data_sf,
                waf_scale = 1/1000, step_waf=4)

        del aux_n34, aux_dmi, aux_n34_2, aux_dmi_2, aux_corr_dmi, aux_corr_n34, \
            aux_corr_dmi_2, aux_corr_n34_2
        # Compute Without Effect -----------------------------------------------
        aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi = \
            ComputeWithoutEffect(data, n34.sel(time=n34.time.dt.month.isin(10)),
                                 dmi.sel(time=dmi.time.dt.month.isin(10)),
                                 seasons[s_count], time_original)

        aux_n34_wodmi_2 = 0
        aux_corr_n34_2 = 0
        aux_dmi_won34_2 = 0
        aux_corr_dmi_2 = 0

        if waf:
            print('#--- WAF form SF ---#')
            aux_waf_aux_n34_wodmi, aux_waf_aux_corr_n34, aux_waf_aux_dmi_won34, \
            aux_waf_aux_corr_dmi = \
                ComputeWithoutEffect(data_sf,
                                     n34.sel(time=n34.time.dt.month.isin(10)),
                                     dmi.sel(time=dmi.time.dt.month.isin(10)),
                                     seasons[s_count],
                                     time_original)

            # N34_wo_DMI WAF
            px_n34_wodmi, py_n34_wodmi = ComputeWaf(aux_waf_aux_n34_wodmi,
                                                    data_clim, hpalevel)

            #DMI_wo_N34 WAF
            px_dmi_won34, py_dmi_won34 = ComputeWaf(aux_waf_aux_dmi_won34,
                                                    data_clim, hpalevel)

        print('Plot...')
        PlotReg(data=aux_n34_wodmi * MakerMaskSig(aux_corr_n34),
                data_cor=aux_corr_n34,
                levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
                title=title_var[v_count] + '_' + s +'_' + str(p[0] + y1) + '_' +
                      str(p[1]) + '_Niño3.4 -{DMI}',
                name_fig=v + '_' + s + str(p[0] + y1) + '_' + str(p[1]) +
                         '_N34_wodmi',
                save=save, sig=False, two_variables=True, data2=aux_n34_wodmi,
                sig2=False, levels2=scales[v_count], SA=SA[v_count], step=1,
                color_map='grey', color_sig='k', sig_point=True, r_crit=r_crit,
                waf=True, px=px_n34_wodmi, py=py_n34_wodmi, data_waf=data_sf,
                waf_scale = 1/1000, step_waf=4)

        PlotReg(data=aux_dmi_won34 * MakerMaskSig(aux_corr_dmi),
                data_cor=aux_corr_dmi,
                levels=scales[v_count], cmap=cmap[v_count], dpi=200,
                title=title_var[v_count] + '_' + s + '_' + str(p[0] + y1) + '_'+
                      str(p[1]) + '_DMI -{N34}',
                name_fig=v + '_' + s + str(p[0] + y1) + '_' + str(p[1]) +
                         '_DMI_woN34',
                save=save, sig=False, two_variables=True, data2=aux_dmi_won34,
                sig2=False, levels2=scales[v_count], SA=SA[v_count], step=1,
                color_map='grey', color_sig='k', sig_point=True, r_crit=r_crit,
                waf=True, px=px_dmi_won34, py=py_dmi_won34, data_waf=data_sf,
                waf_scale=1/2000, step_waf=4)

        del aux_n34_wodmi, aux_dmi_won34, aux_corr_dmi, aux_corr_n34, \
            aux_n34_wodmi_2, aux_dmi_won34_2, aux_corr_dmi_2, aux_corr_n34_2
################################################################################
# Temp y PP#####################################################################
################################################################################
#
def OpenObsDataSet(name, sa=True,
                   dir='/pikachu/datos/luciano.andrian/observado/ncfiles'
                        '/data_obs_d_w_c/'):

    aux = xr.open_dataset(dir + name + '.nc')
    if sa:
        aux2 = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if len(aux2.lat) > 0:
            return aux2
        else:
            aux2 = aux.sel(lon=slice(270, 330), lat=slice(-60, 15))
            return aux2
    else:
        return aux

################################################################################
variables_tpp = ['ppgpcc_w_c_d_1', 'tcru_w_c_d_0.25']
name_var = ['var']
title_var = ['PP GPCC', 'T Cru']
seasons = [7, 10] # main month
seasons_name = ['JJA', 'SON']
SA = [True, False]
aux_name = ['HS', 'SA'] # esto funciona con sa, [True] = 'SA'

scales = [np.linspace(-15, 15, 13),   #pp
          [-.6,-.4,-.2,-.1,-.05,0,0.05,0.1,0.2,0.4,0.6]] #t
cmap = [cbar_pp, cbar]

periodos = [[1940,2020]]
t_critic = 1.66 # es MUY similar (2 digitos) para ambos períodos
y1 = 0
p = periodos[0]
r_crit = np.sqrt(1 / (((np.sqrt((p[1] - p[0] + y1) - 2) / t_critic) ** 2) + 1))

for v, v_count in zip(variables_tpp, [0,1]):

    plt.rcParams['hatch.linewidth'] = 2
    # indices: ----------------------------------------------------------------#
    dmi = dmi_or.sel(
        time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))
    n34 = n34_or.sel(
        time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))

    for sa in SA:
        for s, mm in zip(seasons_name, seasons):
            data = OpenObsDataSet(name=v + '_' + s, sa=sa)
            data = data.sel(
                time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))
            time_original = data.time
            aux_n34, aux_corr_n34, aux_dmi, aux_corr_dmi, aux_n34_2,\
            aux_corr_n34_2, aux_dmi_2, aux_corr_dmi_2 = \
                ComputeWithEffect(data=data, data2=None,
                                  n34=n34.sel(time=n34.time.dt.month.isin(mm)),
                                  dmi=dmi.sel(time=dmi.time.dt.month.isin(mm)),
                                  two_variables=two_variables, m=mm,
                                  full_season=False,
                                  time_original=time_original)

            print('Plot')
            PlotReg(data=aux_n34, data_cor=aux_corr_n34,
                    levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
                    title=title_var[v_count] + '_' + s +
                          '_' + str(p[0] + y1) + '_' + str(p[1]) + '_Niño3.4',
                    name_fig=v + '_' + s + '_' + str(p[0] + y1) +
                             '_' + str(p[1]) + '_N34_' + aux_name[sa],
                    save=save, sig=True, sig_point=True,
                    two_variables=False,
                    SA=sa, step=1,
                    color_map='k', color_sig='k')

            PlotReg(data=aux_dmi, data_cor=aux_corr_dmi,
                    levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
                    title=title_var[v_count] + '_' + s +
                          '_' + str(p[0] + y1) + '_' + str(p[1]) + '_DMI',
                    name_fig=v + '_' + s + '_' + str(p[0] + y1) +
                             '_' + str(p[1]) + '_DMI_' + aux_name[sa],
                    save=save, sig=True, sig_point=True,
                    two_variables=False,
                    SA=sa, step=1, color_map='k', color_sig='k')

            del aux_n34, aux_dmi, aux_n34_2, aux_dmi_2, aux_corr_dmi, \
                aux_corr_n34, aux_corr_dmi_2, aux_corr_n34_2

            aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi = \
                ComputeWithoutEffect(data,
                                     n34.sel(time=n34.time.dt.month.isin(mm)),
                                     dmi.sel(time=dmi.time.dt.month.isin(mm)),
                                     mm, time_original)

            aux_n34_wodmi_2 = 0
            aux_corr_n34_2 = 0
            aux_dmi_won34_2 = 0
            aux_corr_dmi_2 = 0

            print('Plot...')
            PlotReg(data=aux_n34_wodmi, data_cor=aux_corr_n34,
                    levels=scales[v_count], cmap=cmap[v_count], dpi=200,
                    title=title_var[v_count] + '_' + s +
                          '_' + str(p[0] + y1) + '_' + str(
                        p[1]) + '_Niño3.4 -{DMI}',
                    name_fig=v + '_' + s + str(p[0] + y1) + '_' + str(
                        p[1]) + '_N34_wodmi_' + aux_name[sa],
                    save=save, sig=True, sig_point=True,
                    two_variables=False,
                    SA=sa, step=1, color_map='k', color_sig='k')

            PlotReg(data=aux_dmi_won34, data_cor=aux_corr_dmi,
                    levels=scales[v_count], cmap=cmap[v_count], dpi=200,
                    title=title_var[v_count] + '_' + s +
                          '_' + str(p[0] + y1) + '_' + str(
                        p[1]) + '_DMI -{N34}',
                    name_fig=v + '_' + s + str(p[0] + y1) + '_' + str(
                        p[1]) + '_DMI_woN34_' + aux_name[sa],
                    save=save, sig=True, sig_point=True,
                    two_variables=False,
                    SA=sa, step=1, color_map='k', color_sig='k')
################################################################################