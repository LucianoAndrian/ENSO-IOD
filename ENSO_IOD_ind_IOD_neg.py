"""
Anonmalías de HGT y Ks (estado medio a partir del viento) de eventos
individuales negativos "DMI_un_neg" con dmi standard y dmi true-dipole

Se pueden incorporar otras estaciones y cases
"""
################################################################################
data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/1940_2020/'
################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature
import cartopy.feature
import cartopy.crs as ccrs
import pandas as pd
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import sys
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
from ENSO_IOD_Funciones import CompositeSimple, CaseComp
#from ENSO_IOD_Funciones import CompositeSimple

import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
################################################################################
#------------------------------------------------------------------------------#
save = False
seasons = ['SON']
min_max_months = [[9, 11]]
true_dipole = [False, True]

cases = ['DMI_un_neg', 'DMI_un_pos']
title_case = ['IOD pure negative phase ', 'IOD pure positive phase ']
#------------------------------------------------------------------------------#

if save:
    dpi = 300
else:
    dpi = 50

if len(seasons)>1:
    for i in seasons:
        if i != 'SON':
            print('DATOS KS SOLO EN SON!!!!')
            sys.exit(1)
step = 1
v = 'UV200'
v2 = ['HGT200', 'HGT750']
# Functions ####################################################################
def WhatDMI(true_dipole):
    if true_dipole:
        out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/paper1/' \
                  '1940_2020/composite/dmi_true_dipole/'
        nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                      'nc_composites_dates/'
        dmi_index='DMI_true_dipole'
    else:
        out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/paper1/' \
                  '1940_2020/composite/dmi_standard/'
        nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                      'nc_composites_dates_no_ind_sst_anom/'
        dmi_index = 'DMI_standard'
    return out_dir, nc_date_dir, dmi_index


def Subplots(data_array, level, cmap, title, save, dpi, name_fig, out_dir,
             contour=True):

    extent= [0, 359, -80, 10]
    crs_latlon = ccrs.PlateCarree()
    if contour:
        levels_contour = level.copy()
        if isinstance(levels_contour, np.ndarray):
            levels_contour = levels_contour[levels_contour != 0]
        else:
            levels_contour.remove(0)

    time_values = data_array.time.values
    time_steps = len(time_values)
    num_cols = 3
    # cantidad de filas necesarias
    num_rows = np.ceil(time_steps/num_cols).astype(int)

    fig, axes = plt.subplots(
        num_rows, num_cols, figsize=(20, 3*num_rows),
        subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180) })

    for i, (ax, time_val) in enumerate(zip(axes.flatten(), time_values)):

        aux = data_array.sel(time=time_val)
        try:
            aux_var = aux['var']
        except:
            aux_var = aux
        im = ax.contourf(aux.lon, aux.lat, aux_var, levels=level,
                         transform=crs_latlon, cmap=cmap, extend='both')
        if contour:
            ax.contour(aux.lon, aux.lat, aux_var, linewidths=.5, alpha=0.5,
                       levels=levels_contour, transform=crs_latlon,
                       colors='black')

        ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey',
                       edgecolor='k')
        ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
        ax.coastlines(color='k', linestyle='-', alpha=1)
        ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
        ax.set_extent(extent, crs=crs_latlon)
        ax.set_title(f'Time: {pd.Timestamp(time_val).year}')

    # Eliminar los lugares en blanco que existan
    for i in range(time_steps, num_rows * num_cols):
        fig.delaxes(axes.flatten()[i])

    pos = fig.add_axes([0.2, 0.05, 0.6, 0.03])
    cb = fig.colorbar(im, cax=pos, pad=0.2, orientation='horizontal')
    cb.ax.tick_params(labelsize=10)

    fig.suptitle(title, fontsize=16, y=0.98)

    if save:
        plt.savefig(out_dir + name_fig + '.jpg', dpi=dpi, bbox_inches='tight')
        plt.close()
    else:
        plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)
        plt.show()

def extractlevel(x):
    aux = [int(s) for s in x if s.isdigit()]
    return''.join(str(num) for num in aux)

################################################################################
cbar = colors.ListedColormap(['#9B1C00', '#B9391B', '#CD4838', '#E25E55',
                                '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3',
                                '#2064AF', '#014A9B'][::-1])
cbar.set_over('#691800')
cbar.set_under('#013774')
cbar.set_bad(color='white')

cbar_ks = colors.ListedColormap(['#FFE541', '#46B9F4'])
cbar_ks.set_under('#E07700')
cbar_ks.set_over('#016ABA')
cbar_ks.set_bad(color='white')

################################################################################
#for test
td = False
c = 'DMI_un_neg'
c_count=0
s = 'SON'
s_count=0

scale_hgt = [-400, -350, -300, -250, -200, -150, -100, -50,
             0,
             50, 100, 150, 200, 250, 300, 350, 400]
scale_ks = [2, 3, 4]

for td in true_dipole:
    print( 'True-dipole: ' + str(td) + '--------------------------------------')
    out_dir, nc_date_dir, dmi_index = WhatDMI(td)

    # open data ---------------------------------------------------------------#
    try:
        u = xr.open_dataset(data_dir + 'u_' + v + '_w_detrend.nc')
        eta_grad_y = xr.open_dataset(data_dir + 'etay_' + v + '.nc')
        eta_grad_y = eta_grad_y.rename(
            {'meridional_gradient_of_absolute_vorticity': 'var'})

        u = u.interp(lon=eta_grad_y.lon.values, lat=eta_grad_y.lat.values)
        ntime, nlats, nlons = u['var'].shape
    except:
        print('Error en la apertura de datos')

    for c_count, c in enumerate(cases):
        print(c)
        for s_count, s in enumerate(seasons):
            print(s)
            for v2_v in v2:
                level = extractlevel(v2_v)
                aux = xr.open_dataset(nc_date_dir + '1920_2020_' + s + '.nc')
                data_hgt = xr.open_dataset(
                    data_dir + v2_v + '_' + s + '_mer_d_w.nc')

                # Esto termina haciendo lo mismo que el neutro de CaseComp
                # los mmin y mmax ya no son necesarios, de todos modos
                # NO HACEN NADA. Eran de cuando la data entraba sin el rolling
                # y todos los trimestres juntos
                neutro = data_hgt.sel(
                    time=data_hgt.time.dt.year.isin(aux.Neutral)).mean('time')
                data_case = data_hgt.sel(
                    time=data_hgt.time.dt.year.isin(aux[c]))
                data_hgt = data_case - neutro

                print('Plot HGT ----------------------------------------------')
                title = 'z' + level + 'hPa - ' + title_case[c_count]
                name_fig = 'single_events_' + v2_v + '_' + c
                # data_hgt['lon'] = data_hgt.lon.values + 180
                Subplots(data_hgt, scale_hgt, cbar, title, save, dpi, name_fig,
                         out_dir)

                if level == '200':
                    print('Ks ------------------------------------------------')
                    data_y_eta = eta_grad_y. \
                        sel(time=eta_grad_y.time.dt.year.isin(aux[c]))
                    data_y_u = u.sel(time=u.time.dt.year.isin(aux[c]))

                    Rt2 = np.transpose(
                        np.tile(6.378e6 * np.cos(u.lat.values * np.pi / 180),
                                [359, 1]), [1, 0])

                    ks = np.sqrt(
                        data_y_eta['var'][0, :, :] / data_y_u['var']) * Rt2

                    ks = xr.where(ks < 0, np.nan, ks)
                    print('Plot Ks--------------------------------------------')
                    title = 'Ks - ' + title_case[c_count]
                    name_fig = 'single_events_KS_' + c
                    Subplots(ks, scale_ks, cbar_ks, title, save, dpi, name_fig,
                             out_dir, contour=False)

################################################################################
