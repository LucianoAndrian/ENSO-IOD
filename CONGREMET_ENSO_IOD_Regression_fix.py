"""
ENSO vs IOD Regression
T, PP y HGT200
"""
########################################################################################################################
from itertools import groupby
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
import xarray as xr
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import statsmodels.formula.api as sm
import os
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import colors
from ENSO_IOD_Funciones import Nino34CPC
from ENSO_IOD_Funciones import DMI
from ENSO_IOD_Funciones import MakeMask
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

########################################################################################################################
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/CONGREMET/'
hgt_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'

########################################################################################################################
save = True
dpi = 300
full_season = False
text = False
# Functions ############################################################################################################

def OpenDataSet(name, interp=False, lat_interp=None, lon_interp=None):

    if name == 'pp_gpcc':
        # GPCC2018
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/' +
                              'pp_gpcc_d_w_c_1950-2020_1.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(-60, 15))
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
    elif name == 'cmap':
        print('Nop')
        # aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_no_detrend/' +
        #                       'precip.mon.mean.nc')
        # pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        # if interp:
        #     pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        #return pp_gpcc
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


def LinearReg(xrda, dim, deg=1):
    # liner reg along a single dimension
    aux = xrda.polyfit(dim=dim, deg=deg, skipna=True)
    return aux

def LinearReg1_D(dmi, n34):
    import statsmodels.formula.api as smf

    df = pd.DataFrame({'dmi': dmi.values, 'n34': n34.values})

    result = smf.ols(formula='n34~dmi', data=df).fit()
    n34_pred_dmi = result.params[1] * dmi.values + result.params[0]

    result = smf.ols(formula='dmi~n34', data=df).fit()
    dmi_pred_n34 = result.params[1] * n34.values + result.params[0]

    return n34 - n34_pred_dmi, dmi - dmi_pred_n34

def is_months(month, mmin, mmax):
    return (month >= mmin) & (month <= mmax)

def RegWEffect(n34, dmi,data=None, data2=None, m=9,two_variables=False):
    var_reg_n34_2=0
    var_reg_dmi_2=1

    data['time'] = n34
     #print('Full Season')
    aux = LinearReg(data.groupby('month')[m], 'time')
    # aux = xr.polyval(data.groupby('month')[m].time, aux.var_polyfit_coefficients[0]) + \
    #       aux.var_polyfit_coefficients[1]
    var_reg_n34 = aux.var_polyfit_coefficients[0]

    data['time'] = dmi
    aux = LinearReg(data.groupby('month')[m], 'time')
    var_reg_dmi = aux.var_polyfit_coefficients[0]

    if two_variables:
        print('Two Variables')

        data2['time'] = n34
        #print('Full Season data2, m ignored')
        aux = LinearReg(data2.groupby('month')[m], 'time')
        var_reg_n34_2 = aux.var_polyfit_coefficients[0]

        data['time'] = dmi
        aux = LinearReg(data2.groupby('month')[m], 'time')
        var_reg_dmi_2 = aux.var_polyfit_coefficients[0]

    return var_reg_n34, var_reg_dmi, var_reg_n34_2, var_reg_dmi_2

def RegWOEffect(n34, n34_wo_dmi, dmi, dmi_wo_n34, m=9, datos=None):

    datos['time'] = n34

    aux = LinearReg(datos.groupby('month')[m], 'time')
    aux = xr.polyval(datos.groupby('month')[m].time, aux.var_polyfit_coefficients[0]) +\
          aux.var_polyfit_coefficients[1]

    #wo n34
    var_regdmi_won34 = datos.groupby('month')[m]-aux

    var_regdmi_won34['time'] = dmi_wo_n34.groupby('time.month')[m] #index wo influence
    var_dmi_won34 = LinearReg(var_regdmi_won34,'time')

    #-----------------------------------------#

    datos['time'] = dmi
    aux = LinearReg(datos.groupby('month')[m], 'time')
    aux = xr.polyval(datos.groupby('month')[m].time, aux.var_polyfit_coefficients[0]) + \
          aux.var_polyfit_coefficients[1]

    #wo dmi
    var_regn34_wodmi = datos.groupby('month')[m]-aux

    var_regn34_wodmi['time'] = n34_wo_dmi.groupby('time.month')[m] #index wo influence
    var_n34_wodmi = LinearReg(var_regn34_wodmi,'time')

    return var_n34_wodmi.var_polyfit_coefficients[0],\
           var_dmi_won34.var_polyfit_coefficients[0],\
           var_regn34_wodmi,var_regdmi_won34

def Corr(datos, index, time_original, m=9):
    aux_corr1 = xr.DataArray(datos.groupby('month')[m]['var'],
                             coords={'time': time_original.groupby('time.month')[m].values,
                                     'lon': datos.lon.values, 'lat': datos.lat.values},
                             dims=['time', 'lat', 'lon'])
    aux_corr2 = xr.DataArray(index.groupby('time.month')[m],
                             coords={'time': time_original.groupby('time.month')[m]},
                             dims={'time'})

    return xr.corr(aux_corr1, aux_corr2, 'time')

def PlotReg(data, data_cor, levels=np.linspace(-100,100,2), cmap='RdBu_r'
            , dpi=100, save=False, title='\m/', name_fig='fig_PlotReg', sig=True
            ,two_variables = False, data2=None, data_cor2=None, levels2 = np.linspace(-100,100,2)
            , sig2=True, step=1,SA=False, color_map = '#d9d9d9', color_sig='magenta', sig_point=False, SESA=False):

    levels_contour = levels.copy()
    if isinstance(levels_contour, np.ndarray):
        levels_contour = levels_contour[levels_contour != 0]
    else:
        levels_contour.remove(0)
    if SA:
        SESA=False
        fig = plt.figure(figsize=(5, 6), dpi=dpi)
    elif SESA:
        fig = plt.figure(figsize=(5, 6), dpi=dpi)
    else:
        fig = plt.figure(figsize=(7, 3.5), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    if SA:
        ax.set_extent([270,330, -60,20], crs=crs_latlon)
    elif SESA:
        ax.set_extent([295, 315, -40, -12], crs=crs_latlon)
    else:
        ax.set_extent([30, 340, -80, 20], crs=crs_latlon)

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

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.7)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='white', edgecolor=color_map)
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    # ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

    if SA:
        ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
        ax.set_yticks(np.arange(-60, 20, 20), crs=crs_latlon)
        ax.add_feature(cartopy.feature.COASTLINE)
    elif SESA:
        ax.set_xticks(np.arange(295, 315, 5), crs=crs_latlon)
        ax.set_yticks(np.arange(-40, -10, 5), crs=crs_latlon)
        ax.add_feature(cartopy.feature.COASTLINE)
        ax.add_feature(cartopy.feature.BORDERS,  alpha=0.7)
    else:
        ax.set_xticks(np.arange(30, 340, 30), crs=crs_latlon)
        ax.set_yticks(np.arange(-80, 20, 10), crs=crs_latlon)
        ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
        ax.coastlines(color=color_map, linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.minorticks_on()
    ax.tick_params(labelsize=8)
    plt.title(title, fontsize=10)
    plt.tight_layout()

    if save:
        print('save: ' + out_dir + name_fig + '.jpg')
        plt.savefig(out_dir + name_fig + '.jpg', dpi=dpi)
        plt.close()

    else:
        plt.show()


def ComputeWithEffect(data=None, data2=None, n34=None, dmi=None,
                     two_variables=False, full_season=False,
                     time_original=None,m=9):
    print('Reg...')
    print('#-- With influence --#')
    aux_n34, aux_dmi, aux_n34_2, aux_dmi_2 = RegWEffect(data=data, data2=data2,
                                                       n34=n34.__mul__(1 / n34.std('time')),
                                                       dmi=dmi.__mul__(1 / dmi.std('time')),
                                                       m=m, two_variables=two_variables)
    if full_season:
        print('Full Season')
        n34 = n34.rolling(time=5, center=True).mean()
        dmi = dmi.rolling(time=5, center=True).mean()

    print('Corr...')
    aux_corr_n34 = Corr(datos=data, index=n34, time_original=time_original, m=m)
    aux_corr_dmi = Corr(datos=data, index=dmi, time_original=time_original, m=m)

    aux_corr_dmi_2 = 0
    aux_corr_n34_2 = 0
    if two_variables:
        print('Corr2..')
        aux_corr_n34_2 = Corr(datos=data2, index=n34, time_original=time_original, m=m)
        aux_corr_dmi_2 = Corr(datos=data2, index=dmi, time_original=time_original, m=m)

    return aux_n34, aux_corr_n34, aux_dmi, aux_corr_dmi, aux_n34_2, aux_corr_n34_2, aux_dmi_2, aux_corr_dmi_2

def ComputeWithoutEffect(data, n34, dmi, m):
    # -- Without influence --#
    print('# -- Without influence --#')
    print('Reg...')
    # dmi wo n34 influence and n34 wo dmi influence
    dmi_wo_n34, n34_wo_dmi = LinearReg1_D(n34.__mul__(1 / n34.std('time')),
                                          dmi.__mul__(1 / dmi.std('time')))

    # Reg WO
    aux_n34_wodmi, aux_dmi_won34, data_n34_wodmi, data_dmi_won34 = \
        RegWOEffect(n34=n34.__mul__(1 / n34.std('time')),
                   n34_wo_dmi=n34_wo_dmi,
                   dmi=dmi.__mul__(1 / dmi.std('time')),
                   dmi_wo_n34=dmi_wo_n34,
                   m=m, datos=data)

    print('Corr...')
    aux_corr_n34 = Corr(datos=data_n34_wodmi, index=n34_wo_dmi, time_original=time_original,m=m)
    aux_corr_dmi = Corr(datos=data_dmi_won34, index=dmi_wo_n34, time_original=time_original,m=m)

    return aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi

########################################################################################################################
variables = ['pp_gpcc_or', 't_cru', 'hgt200_mer_d_w']#, 'pp_prec']
name_var = ['precip']
title_var = ['PP GPCC', 'T Cru', 'HGT200 ERA5', 'pp_prec']
seasons = [7, 10] # main month
seasons_name = [['SON']]
interp = False
two_variables=False
SA = [False, False, False, False]
SESA = [True, True, False, True]
sig = True

scales = [np.linspace(-15, 15, 13),   #pp
          [-.6,-.4,-.2,-.1,-.05,0,0.05,0.1,0.2,0.4,0.6], #t
          [-150,-100,-75,-50,-25,-15,0,15,25,50,75,100,150], #hgt200
          np.linspace(-15, 15, 13)] #pp


cbar = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar.set_over('#641B00')
cbar.set_under('#012A52')
cbar.set_bad(color='white')


cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC', '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07', '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')

cmap = [cbar_pp, cbar, cbar, cbar_pp]

periodos = [[1920,2020], [1950,2020]]
t_critic = 1.66 # es MUY similar (2 digitos) para ambos períodos

dmi_or = DMI(filter_bwa=False, start_per='1920', end_per='2020')[2]
n34_or = Nino34CPC(xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc"),
                   start=1920, end=2020)[0]

v_count = 0
for v in variables:
    if v != variables[2]:
        plt.rcParams['hatch.linewidth'] = 2
    else:
        plt.rcParams['hatch.linewidth'] = 0.5

    p = periodos[1]
    r_crit = np.sqrt(1 / (((np.sqrt((p[1] - p[0]) - 2) / t_critic) ** 2) + 1))

    # indices: --------------------------------------------------------------------------------------------------------#
    dmi = dmi_or.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))
    n34 = n34_or.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))

    if v == variables[2]:
        data = xr.open_dataset(hgt_dir + v + '.nc')
    else:
        data = OpenDataSet(name=v, interp=False)


    data = data.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))
    time_original = data.time

    if v == 'pp_gpcc' or v == 'pp_gpcc_or' :
        try:
            data = data.rename({name_var[v_count]: 'var'})
        except:
            pass
        dmi = dmi.sel(time=slice(str(p[0]) + '-01-01', str(time_original[-1].dt.year.values) + '-12-01'))
        n34 = n34.sel(time=slice(str(p[0]) + '-01-01', str(time_original[-1].dt.year.values) + '-12-01'))

    # Anomaly ---------------------------------------------------------------------------------------------------------#
    data = data.groupby('time.month') - data.groupby('time.month').mean('time', skipna=True)

    # 3-month running mean --------------------------------------------------------------------------------------------#
    data = data.rolling(time=3, center=True).mean(skipna=True)

    # Seasons ---------------------------------------------------------------------------------------------------------#
    s_count = 1
    #for s in seasons_name:
    s = 'SON'
    aux_n34, aux_corr_n34, aux_dmi, \
    aux_corr_dmi, aux_n34_2, aux_corr_n34_2, \
    aux_dmi_2, aux_corr_dmi_2 = ComputeWithEffect(data=data, data2=None, n34=n34, dmi=dmi,
                                                  two_variables=two_variables, m=seasons[s_count],
                                                  full_season=False, time_original=time_original)
    if v == variables[2]:
        mask = 1
    else:
        mask = MakeMask(aux_n34, 'var')['var']
    print('Plot')
    PlotReg(data=aux_n34*mask, data_cor=aux_corr_n34*mask,
            levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
            title=title_var[v_count] + '_' + s +
                  '_' + str(p[0]) + '_' + str(p[1]) + '_Niño3.4',
            name_fig=v + '_' + s + '_' + str(p[0]) +
                     '_' + str(p[1]) + '_N34',
            save=save, sig=True,
            two_variables=False,
            SA=SA[v_count], step=1,
            color_map='#4B4B4B',
            color_sig='k', sig_point=True,
            SESA=SESA[v_count])

    PlotReg(data=aux_dmi*mask, data_cor=aux_corr_dmi*mask,
            levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
            title=title_var[v_count] + '_' + s +
                  '_' + str(p[0]) + '_' + str(p[1]) + '_DMI',
            name_fig=v + '_' + s + '_' + str(p[0]) +
                     '_' + str(p[1]) + '_DMI',
            save=save, sig=True,
            two_variables=False,
            SA=SA[v_count], step=1,
            color_map='#4B4B4B',
            color_sig='k', sig_point=True,
            SESA=SESA[v_count])

    del aux_n34, aux_dmi, aux_n34_2, aux_dmi_2, aux_corr_dmi, aux_corr_n34, \
        aux_corr_dmi_2, aux_corr_n34_2

    aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi = \
        ComputeWithoutEffect(data, n34, dmi, seasons[s_count])

    aux_n34_wodmi_2 = 0
    aux_corr_n34_2 = 0
    aux_dmi_won34_2 = 0
    aux_corr_dmi_2 = 0

    print('Plot...')
    PlotReg(data=aux_n34_wodmi*mask, data_cor=aux_corr_n34*mask,
            levels=scales[v_count], cmap=cmap[v_count], dpi=dpi,
            title=title_var[v_count] + '_' + s +
                  '_' + str(p[0]) + '_' + str(p[1]) + '_Niño3.4 -{DMI}',
            name_fig=v + '_' + s + str(p[0]) + '_' + str(p[1]) + '_N34_wodmi',
            save=save, sig=True,
            two_variables=False,
            SA=SA[v_count], step=1,
            color_map='#4B4B4B',
            color_sig='k', sig_point=True,
            SESA=SESA[v_count])

    PlotReg(data=aux_dmi_won34*mask, data_cor=aux_corr_dmi*mask,
            levels=scales[v_count], cmap=cmap[v_count], dpi=200,
            title=title_var[v_count] + '_' + s +
                  '_' + str(p[0]) + '_' + str(p[1]) + '_DMI -{N34}',
            name_fig=v + '_' + s + str(p[0]) + '_' + str(p[1]) + '_DMI_woN34',
            save=save, sig=True,
            two_variables=False,
            SA=SA[v_count], step=1,
            color_map='#4B4B4B',
            color_sig='k', sig_point=True,
            SESA=SESA[v_count])

    del aux_n34_wodmi, aux_dmi_won34, aux_corr_dmi, aux_corr_n34, \
        aux_n34_wodmi_2, aux_dmi_won34_2, aux_corr_dmi_2, aux_corr_n34_2

    v_count += 1


