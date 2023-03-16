"""
ENSO vs IOD Regression
T, PP y HGT200
"""
########################################################################################################################
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
########################################################################################################################
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/regression/fix/'
hgt_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'

########################################################################################################################
save = False
dpi = 100
full_season = False
text = False
waf = True
# Functions ############################################################################################################
def OpenDataSet(name, interp=False, lat_interp=None, lon_interp=None):

    if name == 'pp_gpcc':
        # GPCC2018
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/' +
                              'pp_gpcc_d_w_c_1950-2020_1.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        #pp_gpcc = pp_gpcc.rename({'precip': 'var'})

        return pp_gpcc
    elif name == 'pp_gpcc_0.25':
        # GPCC2018
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/' +
                              'pp_gpcc_d_w_c_1950-2020_0.25.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        #pp_gpcc = pp_gpcc.rename({'precip': 'var'})

        return pp_gpcc
    elif name == 'pp_gpcc_or':
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_no_detrend/' +
                              'pp_gpcc_v2020_0.25.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        return pp_gpcc
    elif name == 'pp_cmap':
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/' +
                              'pp_cmap_d_w_c_1979-2020_2.5.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        return pp_gpcc
    elif name == 't_cru':
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/' +
                              't_cru_d_w_c_1920-2020_0.25.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(-60, 15))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        return pp_gpcc
    elif name == 'pp_prec':
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/' +
                              'pp_prec_d_w_c_1950-2020_2.5.nc')
        aux = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            aux = aux.interp(lon=lon_interp, lat=lat_interp)
        return aux
#
# def LinearReg(xrda, dim, deg=1):
#     # liner reg along a single dimension
#     aux = xrda.polyfit(dim=dim, deg=deg, skipna=True)
#     return aux
#
# def LinearReg1_D(dmi, n34):
#     import statsmodels.formula.api as smf
#
#     df = pd.DataFrame({'dmi': dmi.values, 'n34': n34.values})
#
#     result = smf.ols(formula='n34~dmi', data=df).fit()
#     n34_pred_dmi = result.params[1] * dmi.values + result.params[0]
#
#     result = smf.ols(formula='dmi~n34', data=df).fit()
#     dmi_pred_n34 = result.params[1] * n34.values + result.params[0]
#
#     return n34 - n34_pred_dmi, dmi - dmi_pred_n34
#
# def is_months(month, mmin, mmax):
#     return (month >= mmin) & (month <= mmax)
#
# def RegWEffect(n34, dmi,data=None, data2=None, m=9,two_variables=False):
#     var_reg_n34_2=0
#     var_reg_dmi_2=1
#
#     data['time'] = n34
#      #print('Full Season')
#     aux = LinearReg(data.groupby('month')[m], 'time')
#     # aux = xr.polyval(data.groupby('month')[m].time, aux.var_polyfit_coefficients[0]) + \
#     #       aux.var_polyfit_coefficients[1]
#     var_reg_n34 = aux.var_polyfit_coefficients[0]
#
#     data['time'] = dmi
#     aux = LinearReg(data.groupby('month')[m], 'time')
#     var_reg_dmi = aux.var_polyfit_coefficients[0]
#
#     if two_variables:
#         print('Two Variables')
#
#         data2['time'] = n34
#         #print('Full Season data2, m ignored')
#         aux = LinearReg(data2.groupby('month')[m], 'time')
#         var_reg_n34_2 = aux.var_polyfit_coefficients[0]
#
#         data2['time'] = dmi
#         aux = LinearReg(data2.groupby('month')[m], 'time')
#         var_reg_dmi_2 = aux.var_polyfit_coefficients[0]
#
#     return var_reg_n34, var_reg_dmi, var_reg_n34_2, var_reg_dmi_2
#
# def RegWOEffect(n34, n34_wo_dmi, dmi, dmi_wo_n34, m=9, datos=None):
#
#     datos['time'] = n34
#
#     aux = LinearReg(datos.groupby('month')[m], 'time')
#     aux = xr.polyval(datos.groupby('month')[m].time, aux.var_polyfit_coefficients[0]) +\
#           aux.var_polyfit_coefficients[1]
#
#     #wo n34
#     var_regdmi_won34 = datos.groupby('month')[m]-aux
#
#     var_regdmi_won34['time'] = dmi_wo_n34.groupby('time.month')[m] #index wo influence
#     var_dmi_won34 = LinearReg(var_regdmi_won34,'time')
#
#     #-----------------------------------------#
#
#     datos['time'] = dmi
#     aux = LinearReg(datos.groupby('month')[m], 'time')
#     aux = xr.polyval(datos.groupby('month')[m].time, aux.var_polyfit_coefficients[0]) + \
#           aux.var_polyfit_coefficients[1]
#
#     #wo dmi
#     var_regn34_wodmi = datos.groupby('month')[m]-aux
#
#     var_regn34_wodmi['time'] = n34_wo_dmi.groupby('time.month')[m] #index wo influence
#     var_n34_wodmi = LinearReg(var_regn34_wodmi,'time')
#
#     return var_n34_wodmi.var_polyfit_coefficients[0],\
#            var_dmi_won34.var_polyfit_coefficients[0],\
#            var_regn34_wodmi,var_regdmi_won34
#
# def Corr(datos, index, time_original, m=9):
#     aux_corr1 = xr.DataArray(datos.groupby('month')[m]['var'],
#                              coords={'time': time_original.groupby('time.month')[m].values,
#                                      'lon': datos.lon.values, 'lat': datos.lat.values},
#                              dims=['time', 'lat', 'lon'])
#     aux_corr2 = xr.DataArray(index.groupby('time.month')[m],
#                              coords={'time': time_original.groupby('time.month')[m]},
#                              dims={'time'})
#
#     return xr.corr(aux_corr1, aux_corr2, 'time')

def PlotReg(data, data_cor, levels=np.linspace(-100,100,2), cmap='RdBu_r'
            , dpi=100, save=False, title='\m/', name_fig='fig_PlotReg', sig=True
            ,two_variables = False, data2=None, data_cor2=None, levels2 = np.linspace(-100,100,2)
            , sig2=True, step=1,SA=False, color_map = '#d9d9d9', color_sig='magenta', sig_point=False, r_crit=1,
            waf=False, px=None, py=None, aux_waf_scale=1/1000, step_waf=5, data_waf=None):
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

    ax.contour(data.lon[::step], data.lat[::step], data[::step, ::step], linewidths=.5, alpha=0.5,
               levels=levels_contour, transform=crs_latlon, colors='black')

    im = ax.contourf(data.lon[::step], data.lat[::step], data[::step,::step],levels=levels,
                     transform=crs_latlon, cmap=cmap, extend='both')
    if sig:
        if sig_point:
            colors_l = [color_sig, color_sig]
            cs = ax.contourf(data_cor.lon, data_cor.lat, data_cor.where(np.abs(data_cor) > np.abs(r_crit)),
                             transform=crs_latlon, colors='none', hatches=["...", "..."],
                             extend='lower')
            for i, collection in enumerate(cs.collections):
                collection.set_edgecolor(colors_l[i % len(colors_l)])

            for collection in cs.collections:
                collection.set_linewidth(0.)
            # para hgt200 queda mejor los dos juntos
            ax.contour(data_cor.lon[::step], data_cor.lat[::step], data_cor[::step, ::step],
                       levels=np.linspace(-r_crit, r_crit, 2),
                       colors=color_sig, transform=crs_latlon, linewidths=1)

        else:
            ax.contour(data_cor.lon[::step], data_cor.lat[::step], data_cor[::step, ::step],
                       levels=np.linspace(-r_crit, r_crit, 2),
                       colors=color_sig, transform=crs_latlon, linewidths=1)


    if two_variables:
        ax.contour(data2.lon, data2.lat, data2, levels=levels2,
                   colors='k', transform=crs_latlon, linewidths=1)
        if sig2:
            ax.contour(data_cor2.lon, data_cor2.lat, data_cor2, levels=np.linspace(-r_crit, r_crit, 2),
                       colors='forestgreen', transform=crs_latlon, linewidths=1)

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='white', edgecolor=color_map)
    #ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    # ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=
    
    if waf:
        Q60 = np.percentile(np.sqrt(np.add(np.power(px, 2), np.power(py, 2))), 0)
        M = np.sqrt(np.add(np.power(px, 2), np.power(py, 2))) < Q60
        # mask array
        px_mask = ma.array(px, mask=M)
        py_mask = ma.array(py, mask=M)
        # plot vectors
        lons, lats = np.meshgrid(data_waf.lon.values, data_waf.lat.values)
        ax.quiver(lons[::step_waf, ::step_waf],
                  lats[::step_waf, ::step_waf],
                  px_mask[0, ::step_waf, ::step_waf],
                  py_mask[0, ::step_waf, ::step_waf], transform=crs_latlon,pivot='tail',
                  width=0.0020,headwidth=4.1, alpha=1, color='k', scale=aux_waf_scale, scale_units=None)
                  #, scale=1/10)#, width=1.5e-3, headwidth=3.1,  # headwidht (default3)
                  #headlength=2.2)  # (default5))
    
    if SA:
        ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5, zorder=17)
        #ax.add_feature(cartopy.feature.COASTLINE)
        ocean = cartopy.feature.NaturalEarthFeature('physical', 'ocean',
                                                    scale='50m', facecolor='white', alpha=1)
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


# def ComputeWithEffect(data=None, data2=None, n34=None, dmi=None,
#                      two_variables=False, full_season=False,
#                      time_original=None,m=9):
#     print('Reg...')
#     print('#-- With influence --#')
#     aux_n34, aux_dmi, aux_n34_2, aux_dmi_2 = RegWEffect(data=data, data2=data2,
#                                                        n34=n34.__mul__(1 / n34.std('time')),
#                                                        dmi=dmi.__mul__(1 / dmi.std('time')),
#                                                        m=m, two_variables=two_variables)
#     if full_season:
#         print('Full Season')
#         n34 = n34.rolling(time=5, center=True).mean()
#         dmi = dmi.rolling(time=5, center=True).mean()
#
#     print('Corr...')
#     aux_corr_n34 = Corr(datos=data, index=n34, time_original=time_original, m=m)
#     aux_corr_dmi = Corr(datos=data, index=dmi, time_original=time_original, m=m)
#
#     aux_corr_dmi_2 = 0
#     aux_corr_n34_2 = 0
#     if two_variables:
#         print('Corr2..')
#         aux_corr_n34_2 = Corr(datos=data2, index=n34, time_original=time_original, m=m)
#         aux_corr_dmi_2 = Corr(datos=data2, index=dmi, time_original=time_original, m=m)
#
#     return aux_n34, aux_corr_n34, aux_dmi, aux_corr_dmi, aux_n34_2, aux_corr_n34_2, aux_dmi_2, aux_corr_dmi_2
#
# def ComputeWithoutEffect(data, n34, dmi, m):
#     # -- Without influence --#
#     print('# -- Without influence --#')
#     print('Reg...')
#     # dmi wo n34 influence and n34 wo dmi influence
#     dmi_wo_n34, n34_wo_dmi = LinearReg1_D(n34.__mul__(1 / n34.std('time')),
#                                           dmi.__mul__(1 / dmi.std('time')))
#
#     # Reg WO
#     aux_n34_wodmi, aux_dmi_won34, data_n34_wodmi, data_dmi_won34 = \
#         RegWOEffect(n34=n34.__mul__(1 / n34.std('time')),
#                    n34_wo_dmi=n34_wo_dmi,
#                    dmi=dmi.__mul__(1 / dmi.std('time')),
#                    dmi_wo_n34=dmi_wo_n34,
#                    m=m, datos=data)
#
#     print('Corr...')
#     aux_corr_n34 = Corr(datos=data_n34_wodmi, index=n34_wo_dmi, time_original=time_original,m=m)
#     aux_corr_dmi = Corr(datos=data_dmi_won34, index=dmi_wo_n34, time_original=time_original,m=m)
#
#     return aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi


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
########################################################################################################################
# variables = ['pp_gpcc_or', 't_cru', 'hgt200_HS_mer_d_w', 'pp_gpcc_0.25', 'pp_cmap', 'pp_prec']
# name_var = ['precip']
# title_var = ['PP GPCC', 'T Cru', 'HGT200 ERA5', 'PP GPCC 0.25 pr', 'PP CMAP', 'PP PREC']
# seasons = [7, 10] # main month
# seasons_name = ['JJA', 'SON']
# interp = False
# two_variables=False
# SA = [True, True, False, True, True, True]
sig = True

variables = ['hgt200_HS_mer_d_w', 'hgt750_mer_d_w']
name_var = ['precip']
title_var = ['HGT200 ERA5', 'HGT750 ERA5']
seasons = [7, 10] # main month
seasons_name = ['JJA', 'SON']
interp = False
two_variables=False
SA = [False, False]

scales = [[-150,-100,-75,-50,-25,-15,0,15,25,50,75,100,150], [-150,-100,-75,-50,-25,-15,0,15,25,50,75,100,150]]


cbar = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar.set_over('#641B00')
cbar.set_under('#012A52')
cbar.set_bad(color='white')


cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC', '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07', '#543005'][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')

cmap = [cbar, cbar]

periodos = [[1950,2020]]
t_critic = 1.66 # es MUY similar (2 digitos) para ambos períodos

dmi_or = DMI(filter_bwa=False, start_per='1920', end_per='2020')[2]
n34_or = Nino34CPC(xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc"),
                   start=1920, end=2020)[0]

v_count = 0
for v in variables:
    if v == 'pp_cmap':
        y1 = 29
    else:
        y1 = 0
    if v != 'hgt200_HS_mer_d_w' and v != 'hgt750_mer_d_w':
        plt.rcParams['hatch.linewidth'] = 2
        for p in periodos:
            r_crit = np.sqrt(1 / (((np.sqrt((p[1] - p[0]) - 2) / t_critic) ** 2) + 1))

            # indices: ------------------------------------------------------------------------------------------------#
            dmi = dmi_or.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))
            n34 = n34_or.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))

            data = OpenDataSet(name=v, interp=False)
            data = data.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))
            time_original = data.time

            if v == 'pp_gpcc_or':
                data = data.rename({name_var[v_count]: 'var'})
                dmi = dmi.sel(time=slice(str(p[0]) + '-' + '0'+str(data.time[0].dt.month.values) +'-01',
                                         str(time_original[-1].dt.year.values) + '-12-01'))
                n34 = n34.sel(time=slice(str(p[0]) + '-' + '0'+str(data.time[0].dt.month.values) +'-01',
                                         str(time_original[-1].dt.year.values) + '-12-01'))
            elif v == 'pp_gpcc_0.25':
                dmi = dmi.sel(time=slice(str(p[0]) + '-' + '0'+str(data.time[0].dt.month.values) +'-01',
                                         str(time_original[-1].dt.year.values) + '-12-01'))
                n34 = n34.sel(time=slice(str(p[0]) + '-' + '0'+str(data.time[0].dt.month.values) +'-01',
                                         str(time_original[-1].dt.year.values) + '-12-01'))
            elif v == 'pp_cmap':
                dmi = dmi.sel(time=slice('1979-' + '0'+str(data.time[0].dt.month.values) +'-01',
                                         str(time_original[-1].dt.year.values) + '-12-01'))
                n34 = n34.sel(time=slice('1979-'+ '0'+str(data.time[0].dt.month.values) +'-01',
                                         str(time_original[-1].dt.year.values) + '-12-01'))

            # Anomaly -------------------------------------------------------------------------------------------------#
            data = data.groupby('time.month') - data.groupby('time.month').mean('time', skipna=True)

            # 3-month running mean ------------------------------------------------------------------------------------#
            data = data.rolling(time=3, center=True).mean()

            # Seasons -------------------------------------------------------------------------------------------------#
            s_count = 0
            for s in seasons_name:
                aux_n34, aux_corr_n34, aux_dmi, \
                aux_corr_dmi, aux_n34_2, aux_corr_n34_2, \
                aux_dmi_2, aux_corr_dmi_2 = ComputeWithEffect(data=data, data2=None, n34=n34, dmi=dmi,
                                                              two_variables=two_variables, m=seasons[s_count],
                                                              full_season=False, time_original=time_original)

                print('Plot')
                PlotReg(data=aux_n34, data_cor=aux_corr_n34,
                        levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
                        title=title_var[v_count] + '_' + s +
                              '_' + str(p[0] + y1) + '_' + str(p[1]) + '_Niño3.4',
                        name_fig=v + '_' + s + '_' + str(p[0] + y1) +
                                 '_' + str(p[1]) + '_N34',
                        save=save, sig=True,sig_point=True,
                        two_variables=False,
                        SA=SA[v_count], step=1,
                        color_map='k', color_sig='k')

                PlotReg(data=aux_dmi, data_cor=aux_corr_dmi,
                        levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
                        title=title_var[v_count] + '_' + s +
                              '_' + str(p[0] + y1) + '_' + str(p[1]) + '_DMI',
                        name_fig=v + '_' + s + '_' + str(p[0] + y1) +
                                 '_' + str(p[1]) + '_DMI',
                        save=save, sig=True,sig_point=True,
                        two_variables=False,
                        SA=SA[v_count], step=1, color_map='k', color_sig='k')


                del aux_n34, aux_dmi, aux_n34_2, aux_dmi_2, aux_corr_dmi, aux_corr_n34, \
                    aux_corr_dmi_2, aux_corr_n34_2

                aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi = \
                    ComputeWithoutEffect(data, n34, dmi, seasons[s_count])

                aux_n34_wodmi_2 = 0
                aux_corr_n34_2 = 0
                aux_dmi_won34_2 = 0
                aux_corr_dmi_2 = 0

                print('Plot...')
                PlotReg(data=aux_n34_wodmi, data_cor=aux_corr_n34,
                        levels=scales[v_count], cmap=cmap[v_count], dpi=200,
                        title=title_var[v_count] + '_' + s +
                              '_' + str(p[0] + y1) + '_' + str(p[1]) + '_Niño3.4 -{DMI}',
                        name_fig=v + '_' + s + str(p[0] + y1) + '_' + str(p[1]) + '_N34_wodmi',
                        save=save, sig=True,sig_point=True,
                        two_variables=False,
                        SA=SA[v_count], step=1, color_map='k', color_sig='k')

                PlotReg(data=aux_dmi_won34, data_cor=aux_corr_dmi,
                        levels=scales[v_count], cmap=cmap[v_count], dpi=200,
                        title=title_var[v_count] + '_' + s +
                              '_' + str(p[0] + y1) + '_' + str(p[1]) + '_DMI -{N34}',
                        name_fig=v + '_' + s + str(p[0] + y1) + '_' + str(p[1]) + '_DMI_woN34',
                        save=save, sig=True,sig_point=True,
                        two_variables=False,
                        SA=SA[v_count], step=1,  color_map='k', color_sig='k')

                del aux_n34_wodmi, aux_dmi_won34, aux_corr_dmi, aux_corr_n34, \
                    aux_n34_wodmi_2, aux_dmi_won34_2, aux_corr_dmi_2, aux_corr_n34_2

                s_count +=1

    else:
        plt.rcParams['hatch.linewidth'] = 0.5
        p = periodos[0]
        r_crit = np.sqrt(1 / (((np.sqrt((p[1] - p[0] + y1) - 2) / t_critic) ** 2) + 1))
        print("r_crit = " + str(r_crit))
        # indices: ------------------------------------------------------------------------------------------------#
        dmi = dmi_or.sel(time=slice(str(p[0] + y1) + '-01-01', str(p[1]) + '-12-01'))
        n34 = n34_or.sel(time=slice(str(p[0] + y1) + '-01-01', str(p[1]) + '-12-01'))

        data = xr.open_dataset(hgt_dir + v + '.nc')
        data = data.sel(time=slice(str(p[0] + y1) + '-01-01', str(p[1]) + '-12-01'))
        time_original = data.time

        # Anomaly -------------------------------------------------------------------------------------------------#
        data = data.groupby('time.month') - data.groupby('time.month').mean('time', skipna=True)

        # 3-month running mean ------------------------------------------------------------------------------------#
        data = data.rolling(time=3, center=True).mean()

        # --- WAF --- #
        if waf:
            data_sf = xr.open_dataset(hgt_dir + 'sf.nc')
            data_sf = Detrend(OrdenarNC_wTime_fromW(data_sf.rename({'streamfunction': 'var'})), 'time')
            data_sf = data_sf.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))
            # Climatology SON -----------------------------------------------------------------------------------------#
            aux = data_sf.rolling(time=3, center=True).mean()
            data_clim = aux.sel(time=aux.time.dt.month.isin(10)).mean('time')
            # Anomaly -------------------------------------------------------------------------------------------------#
            data_sf = data_sf.groupby('time.month') - data_sf.groupby('time.month').mean('time', skipna=True)
            # 3-month running mean ------------------------------------------------------------------------------------#
            data_sf = data_sf.rolling(time=3, center=True).mean()



        # Seasons -------------------------------------------------------------------------------------------------#
        s_count = 0
        for s in seasons_name:
            aux_n34, aux_corr_n34, aux_dmi, \
            aux_corr_dmi, aux_n34_2, aux_corr_n34_2, \
            aux_dmi_2, aux_corr_dmi_2 = ComputeWithEffect(data=data, data2=None, n34=n34, dmi=dmi,
                                                          two_variables=two_variables, m=seasons[s_count],
                                                          full_season=False, time_original=time_original)

            if waf:
                print('#--- WAF form SF ---#')
                aux_waf_aux_n34, aux_waf_aux_corr_n34, aux_waf_aux_dmi, \
                aux_waf_aux_corr_dmi, aux_waf_aux_n34_2, aux_waf_aux_corr_n34_2, \
                aux_waf_aux_dmi_2, aux_waf_aux_corr_dmi_2 = ComputeWithEffect(data=data_sf, data2=None, n34=n34,
                                                                              dmi=dmi,
                                                                              two_variables=False, m=seasons[s_count],
                                                                              full_season=False,
                                                                              time_original=time_original)

                aux_waf_aux_n34 = MakeXr(aux_waf_aux_n34)
                px_aux_n34, py_aux_n34 = WAF(data_clim, aux_waf_aux_n34, data_clim.lon, data_clim.lat, reshape=True, variable='var')
                px_xr_aux_n34 = xr.DataArray(px_aux_n34[0, :, :])
                py_xr_aux_n34 = xr.DataArray(py_aux_n34[0, :, :])

                aux_waf_aux_dmi = MakeXr(aux_waf_aux_dmi)
                px_aux_dmi, py_aux_dmi = WAF(data_clim, aux_waf_aux_dmi, data_clim.lon, data_clim.lat, reshape=True, variable='var')
                px_xr_aux_dmi = xr.DataArray(px_aux_dmi[0, :, :])
                py_xr_aux_dmi = xr.DataArray(py_aux_dmi[0, :, :])

                weights = np.transpose(np.tile(-2 * np.cos(np.arange(-90, 89) * 1 * np.pi / 180) + 2.1, (359, 1)))
                weights_arr = np.zeros_like(px_aux_n34)
                weights_arr[0, :, :] = weights


            print('Plot')
            PlotReg(data=aux_n34*MakerMaskSig(aux_corr_n34), data_cor=aux_corr_n34,
                    levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
                    title=title_var[v_count] + '_' + s +
                          '_' + str(p[0] + y1) + '_' + str(p[1]) + '_Niño3.4',
                    name_fig=v + '_' + s + '_' + str(p[0] + y1) +
                             '_' + str(p[1]) + '_N34',
                    save=save, sig=False,
                    two_variables=True, data2=aux_n34,
                    sig2=False, levels2=scales[v_count],
                    SA=SA[v_count], step=1,
                    color_map='#grey',
                    color_sig='k', sig_point=True, r_crit=r_crit,
                    waf=True, px=px_aux_n34*weights_arr, py=py_aux_n34*weights_arr, data_waf=data_sf, aux_waf_scale=1/1000)

            PlotReg(data=aux_dmi*MakerMaskSig(aux_corr_dmi), data_cor=aux_corr_dmi,
                    levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
                    title=title_var[v_count] + '_' + s +
                          '_' + str(p[0] + y1) + '_' + str(p[1]) + '_DMI',
                    name_fig=v + '_' + s + '_' + str(p[0] + y1) +
                             '_' + str(p[1]) + '_DMI',
                    save=save, sig=False,
                    two_variables=True, data2=aux_dmi,
                    sig2=False, levels2=scales[v_count],
                    SA=SA[v_count], step=1,
                    color_map='#grey',
                    color_sig='k', sig_point=True, r_crit=r_crit,
                    waf=True, px=px_aux_dmi*weights_arr, py=py_aux_dmi*weights_arr, data_waf=data_sf, aux_waf_scale=1 / 1000)

            del aux_n34, aux_dmi, aux_n34_2, aux_dmi_2, aux_corr_dmi, aux_corr_n34, \
                aux_corr_dmi_2, aux_corr_n34_2

            aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi = \
                ComputeWithoutEffect(data, n34, dmi, seasons[s_count], time_original)

            aux_n34_wodmi_2 = 0
            aux_corr_n34_2 = 0
            aux_dmi_won34_2 = 0
            aux_corr_dmi_2 = 0

            if waf:
                print('#--- WAF form SF ---#')
                aux_waf_aux_n34_wodmi, aux_waf_aux_corr_n34, aux_waf_aux_dmi_won34, aux_waf_aux_corr_dmi = \
                    ComputeWithoutEffect(data_sf, n34, dmi, seasons[s_count], time_original)

                aux_waf_aux_n34 = MakeXr(aux_waf_aux_n34_wodmi)
                px_aux_n34, py_aux_n34 = WAF(data_clim, aux_waf_aux_n34, data_clim.lon, data_clim.lat, reshape=True, variable='var')
                px_xr_aux_n34 = xr.DataArray(px_aux_n34[0, :, :])
                py_xr_aux_n34 = xr.DataArray(py_aux_n34[0, :, :])

                aux_waf_aux_dmi = MakeXr(aux_waf_aux_dmi_won34)
                px_aux_dmi, py_aux_dmi = WAF(data_clim, aux_waf_aux_dmi, data_clim.lon, data_clim.lat, reshape=True, variable='var')
                px_xr_aux_dmi = xr.DataArray(px_aux_dmi[0, :, :])
                py_xr_aux_dmi = xr.DataArray(py_aux_dmi[0, :, :])



            print('Plot...')
            PlotReg(data=aux_n34_wodmi*MakerMaskSig(aux_corr_n34), data_cor=aux_corr_n34,
                    levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
                    title=title_var[v_count] + '_' + s +
                          '_' + str(p[0] + y1) + '_' + str(p[1]) + '_Niño3.4 -{DMI}',
                    name_fig=v + '_' + s + str(p[0] + y1) + '_' + str(p[1]) + '_N34_wodmi',
                    save=save, sig=False,
                    two_variables=True, data2=aux_n34_wodmi,
                    sig2=False, levels2=scales[v_count],
                    SA=SA[v_count], step=1,
                    color_map='#grey',
                    color_sig='k', sig_point=True, r_crit=r_crit,
                    waf=True, px=px_aux_n34 * weights_arr, py=py_aux_n34 * weights_arr, data_waf=data_sf,
                    aux_waf_scale=1 / 1000)

            PlotReg(data=aux_dmi_won34*MakerMaskSig(aux_corr_dmi), data_cor=aux_corr_dmi,
                    levels=scales[v_count], cmap=cmap[v_count], dpi=200,
                    title=title_var[v_count] + '_' + s +
                          '_' + str(p[0] + y1) + '_' + str(p[1]) + '_DMI -{N34}',
                    name_fig=v + '_' + s + str(p[0] + y1) + '_' + str(p[1]) + '_DMI_woN34',
                    save=save, sig=False,
                    two_variables=True, data2=aux_dmi_won34,
                    sig2=False, levels2=scales[v_count],
                    SA=SA[v_count], step=1,
                    color_map='#grey',
                    color_sig='k', sig_point=True, r_crit=r_crit,
                    waf=True, px=px_aux_dmi * weights_arr, py=py_aux_dmi * weights_arr, data_waf=data_sf,
                    aux_waf_scale=1 / 1000)

            del aux_n34_wodmi, aux_dmi_won34, aux_corr_dmi, aux_corr_n34, \
                aux_n34_wodmi_2, aux_dmi_won34_2, aux_corr_dmi_2, aux_corr_n34_2

            s_count += 1
    v_count += 1
