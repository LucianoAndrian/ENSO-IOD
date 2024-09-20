"""
Figuras 12-13 paper ENSO-IOD
BinsByCases funciona de manera retorcida, por eso la parte mas compleja
de este codigo es para reordenar la salida de BinsByCases para poder graficar
con la misma metodologia que el resto de las figuras.
"""
################################################################################
save = True
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/paper1/figuras_final/'
# import #######################################################################
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import xarray as xr
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from matplotlib import colors
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore")

from ENSO_IOD_Funciones import BinsByCases, PlotFinal_Figs12_13

################################################################################
data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/1940_2020/'
index_dir = '/pikachu/datos/luciano.andrian/DMI_N34_Leads_r/'

nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
              'nc_composites_dates_no_ind_sst_anom/' #fechas
data_dir_t_pp = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                'data_obs_d_w_c/' #T y PP ya procesados
sig_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_quantiles/' \
          'DMIbase/' # resultados de MC

cases_dir = "/pikachu/datos/luciano.andrian/cases_fields/"
dates_dir = '/pikachu/datos/luciano.andrian/DMI_N34_Leads_r/'

if save:
    dpi = 300
else:
    dpi = 100

levels = [-375, -275, -175, -75, -25, 0, 25, 75, 175, 275, 375]
scale_hgt_comp = [-375, -275, -175, -75, -25, 0, 25, 75, 175, 275, 375]
scale_hgt_snr = [-1, -.8, -.6, -.5, -.1, 0, 0.1, 0.5, 0.6, 0.8, 1]

cbar = colors.ListedColormap(['#9B1C00', '#B9391B', '#CD4838', '#E25E55',
                              '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3',
                              '#2064AF', '#014A9B'][::-1])
cbar.set_over('#641B00')
cbar.set_under('#012A52')
cbar.set_bad(color='white')


cbar_snr = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#6FFE9B',
                                  '#FFFFFF',
                                  '#FFFFFF', '#FFFFFF',
                                  '#FEB77E','#CA3E72','#782281','#251255'])
cbar_snr.set_over('#251255')
cbar_snr.set_under('#070B4F')
cbar_snr.set_bad(color='white')

################################################################################
cases = ['dmi_puros_pos', 'dmi_puros_neg',
        'n34_puros_pos', 'n34_puros_neg',
        'sim_pos', 'sim_neg',
        'dmi_neg_n34_pos', 'dmi_pos_n34_neg',
        'neutros']

row_titles = [None, None, 'Strong EN', None, None, None,
              None, 'Moderate EN', None, None, 'Neutro ENSO',
              None, None, None, None, 'Moderate LN', None, None,
              None, None, 'Strong LN']

col_titles = [None, None, 'Neutro IOD', 'Moderate IOD+',
              'Strong IOD+', None, None, None, None, None,
              'Strong IOD-', 'Moderate IOD-']

lat = np.linspace(-80, 20, 101)
lon = np.linspace(0, 359, 360)
################################################################################
# por motivos esteticos...
clim = xr.open_dataset('/pikachu/datos/luciano.andrian/'
                       'val_clim_cfsv2/hindcast_cfsv2_meanclim_son.nc')
clim = clim.rename({'hgt':'var'})
clim = clim.drop('P')
################################################################################

# Orden de ploteo ------------------------------------------------------------ #
cases_magnitude = [None, None, 's_en', 's_en-m_iodp', 's_en-s_iodp',
                   None, None, 'm_en', 'm_en-m_iodp', 'm_en-s_iodp',
                   's_iodn', 'm_iodn', 'clim', 'm_iodp', 's_iodp',
                   'm_ln-s_iodn', 'm_ln-m_iodn', 'm_ln', None, None,
                   's_ln-s_iodn', 's_ln-m_iodn', 's_ln', None, None]

bin_limits = [[-4.5,-1], #0 s
              [-1, -0.5], #1 m
              [-0.5, 0.5], #2 -
              [0.5, 1], #3  m
              [1, 4.5]] #4 s

bins_by_cases_dmi = [[3, 4], [0, 1],
                     [2], [2],
                     [3, 4], [0, 1],
                     [0, 1],[3, 4],
                     [2]]

bins_by_cases_n34 = [[2], [2],
                     [3, 4], [0,1],
                     [3, 4], [0, 1],
                     [3, 4], [0, 1],
                     [2]]

# Descubriendo el orden de computo original ---------------------------------- #
# y asignando un nombre para poder usar cases_magnitud definido arriba
bin_names = ['s', 'm', '', 'm', 's']
cases_names = []
for c_count, c  in enumerate(cases):
    aux_h = '-'
    for d in bins_by_cases_dmi[c_count]:
        d_aux = sum(bin_limits[d])
        d_aux_mag_name = bin_names[d]
        d_aux_h = '_'
        if d_aux>0:
            d_aux_name = 'iodp'
        elif d_aux<0:
            d_aux_name = 'iodn'
        elif d_aux==0:
            d_aux_name = ''
            d_aux_mag_name = ''
            d_aux_h = ''
            aux_h = ''

        iod_name = f"{d_aux_mag_name}{d_aux_h}{d_aux_name}"

        for n in bins_by_cases_n34[c_count]:
            n_aux = sum(bin_limits[n])
            n_aux_mag_name = bin_names[n]
            n_aux_h = '_'

            if n_aux > 0:
                n_aux_name = 'en'
            elif n_aux < 0:
                n_aux_name = 'ln'
            elif n_aux == 0:
                n_aux_name = ''
                n_aux_mag_name=''
                n_aux_h = ''
                aux_h = ''

            enso_name = f"{n_aux_mag_name}{n_aux_h}{n_aux_name}"
            case_name = f"{enso_name}{aux_h}{iod_name}"
            cases_names.append(case_name)

print('#######################################################################')
print('Figure12')
print('#######################################################################')
# Computo en el orden orignal ------------------------------------------------ #
aux_comps = {}
aux_num_comps = {}
n_count = 0
for c_count, c  in enumerate(cases):
    cases_bin, num_bin, auxx = BinsByCases(v='hgt', v_name='hgt',
                                           fix_factor=9.8, s='SON', mm=10, c=c,
                                           c_count=c_count,
                                           bin_limits=bin_limits,
                                           bins_by_cases_dmi=bins_by_cases_dmi,
                                           bins_by_cases_n34=bins_by_cases_n34,
                                           snr=False, cases_dir=cases_dir,
                                           dates_dir=dates_dir,
                                           neutro_clim=True)

    bins_aux_dmi = bins_by_cases_dmi[c_count]
    bins_aux_n34 = bins_by_cases_n34[c_count]

    for b_dmi in range(0, len(bins_aux_dmi)):
        for b_n34 in range(0, len(bins_aux_n34)):
            aux_comps[cases_names[n_count]] = cases_bin[b_dmi][b_n34]
            aux_num_comps[cases_names[n_count]] = num_bin[b_dmi][b_n34]
            n_count += 1
# ---------------------------------------------------------------------------- #
# Reordenando para plotear --------------------------------------------------- #
cases_ordenados = []

aux_num=[]
for c in cases_magnitude:
    try:
        aux = aux_comps[c]
        aux_num.append(aux_num_comps[c])
    except:
        aux = aux_comps[cases_magnitude[2]]*0
        aux_num.append('')


    da = xr.DataArray(aux, dims=["lat", "lon"], coords={"lat": lat, "lon": lon},
                      name="var")
    cases_ordenados.append(da)

cases_ordenados = xr.concat(cases_ordenados, dim='plots')
# ---------------------------------------------------------------------------- #

PlotFinal_Figs12_13(cases_ordenados, scale_hgt_comp, cbar, aux_num, 'f12',
                    'hs', save, dpi, out_dir,
                    data_ctn=cases_ordenados, color_ctn='#505050',
                    row_titles=row_titles, col_titles=col_titles,
                    clim_plot=clim, clim_cbar='Reds',
                    clim_levels=np.linspace(11000,12500,7), high=0.55)

print('#######################################################################')
print('Figure13')
print('#######################################################################')
# Computo en el orden orignal ------------------------------------------------ #
aux_comps = {}
aux_num_comps = {}
n_count = 0
for c_count, c  in enumerate(cases):
    cases_bin, num_bin, auxx = BinsByCases(v='hgt', v_name='hgt',
                                           fix_factor=9.8, s='SON', mm=10, c=c,
                                           c_count=c_count,
                                           bin_limits=bin_limits,
                                           bins_by_cases_dmi=bins_by_cases_dmi,
                                           bins_by_cases_n34=bins_by_cases_n34,
                                           snr=True, cases_dir=cases_dir,
                                           dates_dir=dates_dir,
                                           neutro_clim=True)

    bins_aux_dmi = bins_by_cases_dmi[c_count]
    bins_aux_n34 = bins_by_cases_n34[c_count]

    for b_dmi in range(0, len(bins_aux_dmi)):
        for b_n34 in range(0, len(bins_aux_n34)):
            aux_comps[cases_names[n_count]] = cases_bin[b_dmi][b_n34]
            aux_num_comps[cases_names[n_count]] = num_bin[b_dmi][b_n34]
            n_count += 1
# ---------------------------------------------------------------------------- #
# Reordenando para plotear --------------------------------------------------- #
cases_ordenados = []

aux_num=[]
for c in cases_magnitude:
    try:
        aux = aux_comps[c]
        aux_num.append(aux_num_comps[c])
    except:
        aux = aux_comps[cases_magnitude[2]]*0
        aux_num.append('')


    da = xr.DataArray(aux['var'], dims=["lat", "lon"],
                      coords={"lat": lat, "lon": lon},
                      name="var")
    cases_ordenados.append(da)

cases_ordenados = xr.concat(cases_ordenados, dim='plots')
# ---------------------------------------------------------------------------- #

PlotFinal_Figs12_13(cases_ordenados, scale_hgt_snr, cbar_snr, aux_num,
                    'f13',
                    'hs', save, dpi, out_dir,
                    data_ctn=cases_ordenados, color_ctn='#505050',
                    row_titles=row_titles, col_titles=col_titles,
                    clim_plot=clim, clim_cbar='Reds',
                    clim_levels=np.linspace(11000,12500,7), high=0.55)
print('#######################################################################')
print('Done')
print('#######################################################################')