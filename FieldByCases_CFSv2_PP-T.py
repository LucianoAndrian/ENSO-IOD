"""
BinsByCases funciona de manera retorcida, por eso la parte mas compleja
de este codigo es para reordenar la salida de BinsByCases para poder graficar
con la misma metodologia que el resto de las figuras.
"""
################################################################################
save = False
out_dir = ('/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Composites/'
           'fieldbycases/')
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

from ENSO_IOD_Funciones import BinsByCases, PlotFinal_Figs12_13, MakeMask

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

scale_pp = np.linspace(-30, 30, 13)
scale_t = [-1, -.75, -.5, -.25, -.1, 0, .1, .25, .5, .75, 1]
scale_snr = [-1, -.8, -.6, -.5, -.1, 0, 0.1, 0.5, 0.6, 0.8, 1]

cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC',
                                 '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07',
                                 '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')

cbar_pp_snr = colors.ListedColormap(['#002A3D','#074D4F', '#1E6D5A' ,'#52C39D',
                                     '#6FFE9B',
                                  '#FFFFFF',
                                  '#DCBC75', '#995D13','#6A3D07','#543005',
                                     '#3F2404'][::-1])
cbar_pp_snr.set_under('#3F2404')
cbar_pp_snr.set_over('#002A3D')
cbar_pp_snr.set_bad(color='white')


cbar_t = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89',
                                '#FFCECC', 'white', '#B3DBFF', '#83B9EB',
                                '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_t.set_over('#9B1C00')
cbar_t.set_under('#014A9B')
cbar_t.set_bad(color='white')

cbar_snr_t = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#52C39D',
                                  '#6FFE9B', '#FFFFFF', '#FEB77E', '#FB8761',
                                  '#CA3E72','#782281','#251255'])
cbar_snr_t.set_over('#251255')
cbar_snr_t.set_under('#070B4F')
cbar_snr_t.set_bad(color='white')


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

lat = np.arange(-60, 15+1)
lon = np.arange(275, 330+1)
lat = np.arange(-60, 15+1)
lon = np.arange(275, 330+1)
################################################################################
# por motivos esteticos...
clim = (xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/'
                       'data_no_detrend/pp_gpcc_v2020_0.25.nc').
        sel(lat=slice(lat[-1], lat[0]+1), lon=slice(lon[0], lon[-1]+1)))
clim = clim.sel(time=clim.time[0])*np.nan
clim = clim.rename({'precip':'var'})
# clim = clim.rolling(time=3, center=True).mean()
# clim = clim.sel(time=clim.time.dt.month.isin(10))
# clim = clim.rename({'precip':'var'})
# clim = clim.mean(['time'])

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
print('PREC ##################################################################')
print('#######################################################################')

# Computo en el orden orignal ------------------------------------------------ #
aux_comps = {}
aux_num_comps = {}
n_count = 0
for c_count, c  in enumerate(cases):
    cases_bin, num_bin, auxx = BinsByCases(v='prec', v_name='prec',
                                           fix_factor=30, s='SON', mm=10, c=c,
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

            try:
                aux_comps[cases_names[n_count]] = cases_bin['var']
                aux_num_comps[cases_names[n_count]] = np.nan
            except:
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

    try:
        da = xr.DataArray(aux, dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")
    except:
        da = xr.DataArray(aux['var'], dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")

    cases_ordenados.append(da)

cases_ordenados = xr.concat(cases_ordenados, dim='plots')
# ---------------------------------------------------------------------------- #

cases_ordenados = cases_ordenados*MakeMask(cases_ordenados, 'var')['var']
PlotFinal_Figs12_13(cases_ordenados, scale_pp, cbar_pp, aux_num,
                    'fieldbycases_prec_SON-CFSv2',
                    'SA', save, dpi, out_dir,
                    data_ctn=cases_ordenados, color_ctn='k',
                    row_titles=row_titles, col_titles=col_titles,
                    clim_plot=clim, clim_cbar='YlGnBu',
                    clim_levels=np.linspace(0,500,10), high=2, width=7,
                    plot_step=1, contourf_clim=False)

print('PREC - SNR#############################################################')
# Computo en el orden orignal ------------------------------------------------ #
aux_comps = {}
aux_num_comps = {}
n_count = 0
for c_count, c  in enumerate(cases):
    cases_bin, num_bin, auxx = BinsByCases(v='prec', v_name='prec',
                                           fix_factor=30, s='SON', mm=10, c=c,
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

            try:
                aux_comps[cases_names[n_count]] = cases_bin['var']
                aux_num_comps[cases_names[n_count]] = np.nan
            except:
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

    try:
        da = xr.DataArray(aux, dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")
    except:
        da = xr.DataArray(aux['var'], dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")

    cases_ordenados.append(da)

cases_ordenados = xr.concat(cases_ordenados, dim='plots')
# ---------------------------------------------------------------------------- #

cases_ordenados = cases_ordenados*MakeMask(cases_ordenados, 'var')['var']
PlotFinal_Figs12_13(cases_ordenados, scale_snr, cbar_pp_snr, aux_num,
                    'fieldbycases_prec_SON-CFSv2',
                    'SA', save, dpi, out_dir,
                    data_ctn=cases_ordenados, color_ctn='k',
                    row_titles=row_titles, col_titles=col_titles,
                    clim_plot=clim, clim_cbar='YlGnBu',
                    clim_levels=np.linspace(0,500,10), high=2, width=7,
                    plot_step=1, contourf_clim=False)

print('#######################################################################')
print('Tref ##################################################################')
print('#######################################################################')
lat = np.arange(-60, 15+1)
# Computo en el orden orignal ------------------------------------------------ #
aux_comps = {}
aux_num_comps = {}
n_count = 0
for c_count, c  in enumerate(cases):
    cases_bin, num_bin, auxx = BinsByCases(v='tref', v_name='tref',
                                           fix_factor=1, s='SON', mm=10, c=c,
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
            try:
                aux_comps[cases_names[n_count]] = cases_bin['var']
                aux_num_comps[cases_names[n_count]] = np.nan
            except:
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

    try:
        da = xr.DataArray(aux, dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")
    except:
        da = xr.DataArray(aux['var'], dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")

    cases_ordenados.append(da)

cases_ordenados = xr.concat(cases_ordenados, dim='plots')
# ---------------------------------------------------------------------------- #

cases_ordenados = cases_ordenados*MakeMask(cases_ordenados, 'var')['var']
PlotFinal_Figs12_13(cases_ordenados, scale_t, cbar_t, aux_num,
                    'fieldbycases_tref_SON-CFSv2',
                    'SA', save, dpi, out_dir,
                    data_ctn=cases_ordenados, color_ctn='k',
                    row_titles=row_titles, col_titles=col_titles,
                    clim_plot=clim, clim_cbar='YlGnBu',
                    clim_levels=np.linspace(0,500,10), high=2, width=7,
                    plot_step=1, contourf_clim=False)

print('Tref - SNR ############################################################')
# Computo en el orden orignal ------------------------------------------------ #
aux_comps = {}
aux_num_comps = {}
n_count = 0
for c_count, c  in enumerate(cases):
    cases_bin, num_bin, auxx = BinsByCases(v='tref', v_name='tref',
                                           fix_factor=1, s='SON', mm=10, c=c,
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
            try:
                aux_comps[cases_names[n_count]] = cases_bin['var']
                aux_num_comps[cases_names[n_count]] = np.nan
            except:
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

    try:
        da = xr.DataArray(aux, dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")
    except:
        da = xr.DataArray(aux['var'], dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")

    cases_ordenados.append(da)

cases_ordenados = xr.concat(cases_ordenados, dim='plots')
# ---------------------------------------------------------------------------- #

cases_ordenados = cases_ordenados*MakeMask(cases_ordenados, 'var')['var']
PlotFinal_Figs12_13(cases_ordenados, scale_snr, cbar_snr_t, aux_num,
                    'fieldbycases_tref_SON-CFSv2',
                    'SA', save, dpi, out_dir,
                    data_ctn=cases_ordenados, color_ctn='k',
                    row_titles=row_titles, col_titles=col_titles,
                    clim_plot=clim, clim_cbar='YlGnBu',
                    clim_levels=np.linspace(0,500,10), high=2, width=7,
                    plot_step=1, contourf_clim=False)


print('#######################################################################')
print('Tsigma ################################################################')
print('#######################################################################')
lat = np.arange(-60, 15+1)[::-1]
# Computo en el orden orignal ------------------------------------------------ #
aux_comps = {}
aux_num_comps = {}
n_count = 0
for c_count, c  in enumerate(cases):
    cases_bin, num_bin, auxx = BinsByCases(v='tsigma', v_name='TMP',
                                           fix_factor=1, s='SON', mm=10, c=c,
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
            try:
                aux_comps[cases_names[n_count]] = cases_bin['var']
                aux_num_comps[cases_names[n_count]] = np.nan
            except:
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


    try:
        da = xr.DataArray(aux, dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")
    except:
        da = xr.DataArray(aux['var'], dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")

    cases_ordenados.append(da)

cases_ordenados = xr.concat(cases_ordenados, dim='plots')
# ---------------------------------------------------------------------------- #

cases_ordenados = cases_ordenados*MakeMask(cases_ordenados, 'var')['var']
PlotFinal_Figs12_13(cases_ordenados, scale_t, cbar_t, aux_num,
                    'fieldbycases_tsigma_SON-CFSv2',
                    'SA', save, dpi, out_dir,
                    data_ctn=cases_ordenados, color_ctn='k',
                    row_titles=row_titles, col_titles=col_titles,
                    clim_plot=clim, clim_cbar='YlGnBu',
                    clim_levels=np.linspace(0,500,10), high=2, width=7,
                    plot_step=1, contourf_clim=False)

print('Tsigma - SNR ##########################################################')
# Computo en el orden orignal ------------------------------------------------ #
aux_comps = {}
aux_num_comps = {}
n_count = 0
for c_count, c  in enumerate(cases):
    cases_bin, num_bin, auxx = BinsByCases(v='tsigma', v_name='TMP',
                                           fix_factor=1, s='SON', mm=10, c=c,
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
            try:
                aux_comps[cases_names[n_count]] = cases_bin['var']
                aux_num_comps[cases_names[n_count]] = np.nan
            except:
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

    try:
        da = xr.DataArray(aux, dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")
    except:
        da = xr.DataArray(aux['var'], dims=["lat", "lon"],
                          coords={"lat": lat, "lon": lon},
                          name="var")

    cases_ordenados.append(da)

cases_ordenados = xr.concat(cases_ordenados, dim='plots')
# ---------------------------------------------------------------------------- #

cases_ordenados = cases_ordenados*MakeMask(cases_ordenados, 'var')['var']
PlotFinal_Figs12_13(cases_ordenados, scale_snr, cbar_snr_t, aux_num,
                    'fieldbycases_tsigma_SON-CFSv2',
                    'SA', save, dpi, out_dir,
                    data_ctn=cases_ordenados, color_ctn='k',
                    row_titles=row_titles, col_titles=col_titles,
                    clim_plot=clim, clim_cbar='YlGnBu',
                    clim_levels=np.linspace(0,500,10), high=2, width=7,
                    plot_step=1, contourf_clim=False)


print('#######################################################################')
print('Done ##################################################################')
print('#######################################################################')
