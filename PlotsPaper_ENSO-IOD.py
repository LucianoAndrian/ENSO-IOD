"""
Figuras paper ENSO-IOD
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

from ENSO_IOD_Funciones import (Nino34CPC, DMI, DMI2, ComputeWithEffect, WAF,
                                ComputeWithoutEffect, SetDataToPlotFinal,
                                PlotFinal, CaseComp)
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

if save:
    dpi = 300
else:
    dpi = 100
# Funciones ####################################################################

def auxScatter(n34, n34_3, dmi, dmi_3, s):
    dmi_todos = dmi_3.sel(time=dmi_3.time.dt.month.isin([s]))
    dmi_criteria_y = dmi.where((dmi.Mes == s)).Años.dropna().values

    n34_todos = n34.sel(time=n34.time.dt.month.isin([s]))
    n34_criteria_y = n34_3.where((n34_3.Mes == s)).Años.dropna().values

    sim_y = np.intersect1d(n34_criteria_y, dmi_criteria_y)

    dmi_sim = dmi_todos.sel(time=dmi_todos.time.dt.year.isin(sim_y))
    n34_sim = n34_todos.sel(time=n34_todos.time.dt.year.isin(sim_y))

    dmi_sim_pos = dmi_sim.where(dmi_sim > 0)
    n34_sim_pos = n34_sim.where(n34_sim > 0)

    dmi_sim_neg = dmi_sim.where(dmi_sim < 0)
    n34_sim_neg = n34_sim.where(n34_sim < 0)

    dmi_pos_n34_neg = dmi_sim_pos.where(~np.isnan(n34_sim_neg.values))
    dmi_neg_n34_pos = dmi_sim_neg.where(~np.isnan(n34_sim_pos.values))

    dmi_dates_ref = dmi_todos.time.dt.year
    mask = np.in1d(dmi_dates_ref, dmi_criteria_y)
    aux_dmi = dmi_todos.sel(
        time=dmi_todos.time.dt.year.isin(dmi_dates_ref[mask]))

    n34_dates_ref = n34_todos.time.dt.year
    mask = np.in1d(n34_dates_ref, n34_criteria_y)
    aux_n34 = n34_todos.sel(
        time=n34_todos.time.dt.year.isin(n34_dates_ref[mask]))

    aux_dates_ref = aux_dmi.time.dt.year
    mask = np.in1d(aux_dates_ref, sim_y, invert=True)
    dmi_un = aux_dmi.sel(
        time=aux_dmi.time.dt.year.isin(aux_dates_ref[mask]))

    dmi_un_pos = dmi_un.where(dmi_un > 0)
    dmi_un_pos_n34_values = n34_todos.sel(
        time=n34_todos.time.isin(dmi_un_pos.time))

    dmi_un_neg = dmi_un.where(dmi_un < 0)
    dmi_un_neg_n34_values = n34_todos.sel(
        time=n34_todos.time.isin(dmi_un_neg.time))

    aux_dates_ref = aux.time.dt.year
    mask = np.in1d(aux_dates_ref, sim_y, invert=True)
    n34_un = aux_n34.sel(
        time=aux_n34.time.dt.year.isin(aux_dates_ref[mask]))

    n34_un_pos = n34_un.where(n34_un > 0)
    n34_un_pos_dmi_values = dmi_todos.sel(
        time=dmi_todos.time.isin(n34_un_pos.time))
    n34_un_neg = n34_un.where(n34_un < 0)
    n34_un_neg_dmi_values = dmi_todos.sel(
        time=dmi_todos.time.isin(n34_un_neg.time))

    return dmi_un_pos, dmi_un_pos_n34_values, dmi_un_neg, \
           dmi_un_neg_n34_values, n34_un_pos, n34_un_pos_dmi_values,\
           n34_un_neg, n34_un_neg_dmi_values, dmi_sim_pos, n34_sim_pos, \
           dmi_sim_neg, n34_sim_neg, dmi_todos, n34_todos, dmi_pos_n34_neg, \
           dmi_neg_n34_pos

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

def SigDivMask(data, thr):
    return xr.where(np.abs(data) < thr, np.nan, 1)

def MakerMaskSig(data):
    mask_sig = data.where((data < -1 * r_crit) | (data > r_crit))
    mask_sig = mask_sig.where(np.isnan(mask_sig), 1)

    return mask_sig

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

    weights = np.transpose(
        np.tile(-2 *
                np.cos(data_clim.lat.values * 1 * np.pi / 180) + 2.1, (359, 1)))
    weights_arr = np.zeros_like(px)
    weights_arr[0, :, :] = weights
    px *= weights_arr
    py *= weights_arr

    return px, py


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

def SelectParIndex(dmi_case, n34_case, sd_dmi_s, sd_n34_s, s, by_r=True,
                   open_n34=True, open_dmi=False):
    aux = xr.open_dataset(
        cases_dir + 'dmi_values_' + dmi_case + '_' + s + '.nc')\
        .__mul__(1 / sd_dmi_s)
    aux2 = xr.open_dataset(
        cases_dir + 'N34_values_' + n34_case + '_' + s + '.nc')\
        .__mul__(1 / sd_n34_s)

    if by_r:
        if open_n34:
            n34_s = xr.open_dataset(
                index_dir + 'N34_' + s + '_Leads_r_CFSv2.nc')\
                .__mul__(1/sd_n34_s)
            aux3 = n34_s.sel(r=aux.r, time=aux.time)

            if len(np.where(aux3.L.values==aux.L.values)[0]):
                #return aux.sst.values, aux3.sst.values
                return aux.sst.values.round(2), aux3.sst.values.round(2)
            else:
                print('Error: CASES')
                return [], []

        if open_dmi:
            dmi_s = xr.open_dataset(
                index_dir + 'DMI_' + s + '_Leads_r_CFSv2.nc')\
                .__mul__(1/sd_dmi_s)
            aux3 = dmi_s.sel(r=aux2.r, time=aux2.time)

            if len(np.where(aux3.L.values == aux2.L.values)[0]):
                return aux3.sst.values.round(2), aux2.sst.values.round(2)
            else:
                print('Error: CASES')
                return [], []
    else:
        aux2 = aux2.sel(time=aux2.time.isin([aux.time.values]))

        if len(aux2.time) == len(aux.time):
            return aux.sst.values.round(2), aux2.sst.values.round(2)
        else:
            print('Error: CASES')
            return [], []

def Weights(data):
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180),
                                   (len(data.lon), 1)))
    data_w = data * weights
    return data_w

# Escalas y barras de colores ##################################################
# Regresion ------------------------------------------------------------------ #
scale_vp = [-3e6, -2.5e6, -2e6, -1.5e6, -1e6, -0.5e6, 0, 0.5e6, 1e6, 1.5e6,
            2e6, 2.5e6, 3e6]
scale_sst = [-1, -.5, -.1, 0, .1, .5, 1]
scale_div = [-4.33e-07, 4.33e-07]
scale_hgt = [-150, -100, -75, -50, -25, -15, 0, 15, 25, 50, 75, 100, 150]

# Composite ------------------------------------------------------------------ #
scale_vp_comp = np.linspace(-4.5e6, 4.5e6, 13)
scale_div_comp = [-1.6e-06, 1.6e-06]
scale_sst_comp = [-1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5]
scale_hgt_comp = [-375, -275, -175, -75, -25, 0, 25, 75, 175, 275, 375]

# CFSv2
scale_hgt_snr = [-1, -.8, -.6, -.5, -.1, 0, 0.1, 0.5, 0.6, 0.8, 1]
# ---------------------------------------------------------------------------- #
cbar_sst = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89',
                                  '#FFCECC', 'white', '#B3DBFF', '#83B9EB',
                                  '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_sst.set_over('#9B1C00')
cbar_sst.set_under('#014A9B')
cbar_sst.set_bad(color='white')


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

# Variables comunes ---------------------------------------------------------- #
#sig = True
periodos = [[1940, 2020]]
y1 = 0
p = periodos[0]
s = 'SON'
s_count = 0

# scatter
in_label_size = 18
label_legend_size = 18
tick_label_size = 15
scatter_size_fix = 3

################################################################################
subtitulos_regre = [r"$Ni\tilde{n}o\ 3.4$",  r"$DMI$",
                    r"$Ni\tilde{n}o\ 3.4|_{DMI}$",
                    r"$DMI|_{Ni\tilde{n}o\ 3.4}$"]

# En orden de ploteo
cases = ['N34_un_pos', 'N34_un_neg', 'DMI_un_pos', 'DMI_un_neg',
         'DMI_sim_pos', 'DMI_sim_neg']

cases_cfsv2 =['n34_puros_pos', 'n34_puros_neg',
              'dmi_puros_pos', 'dmi_puros_neg',
              'sim_pos', 'sim_neg']

title_case = ['Pure El Niño', 'Pure La Niña',
              'Pure positive IOD ', 'Pure negative IOD',
              'El Niño - positive IOD', 'La Niña - negative IOD']
################################################################################
print('#######################################################################')
print('Figure1')
print('#######################################################################')
dmi, dmi_2, dmi_3 = DMI2(filter_bwa=False, start_per='1920', end_per='2020',
                         sst_anom_sd=False, opposite_signs_criteria=False)
dmi_3 = dmi_3.sel(time=slice('1940-01-01', '2020-12-01'))
dmi = dmi.where(dmi.Años>=1940).dropna()
dmi_3 = dmi_3 / dmi_3.std('time')

aux = xr.open_dataset("/pikachu/datos/luciano.andrian/verif_2019_2023/"
                      "sst.mnmean.nc")
n34, n34_2, n34_3 = Nino34CPC(aux, start=1920)
n34 = n34.sel(time=slice('1940-01-01', '2020-12-01'))
n34_3 = n34_3.where(n34_3.Años>=1940).dropna()
n34 = n34 / n34.std('time')

dmi_un_pos, dmi_un_pos_n34_values, dmi_un_neg, dmi_un_neg_n34_values, \
n34_un_pos, n34_un_pos_dmi_values, n34_un_neg, n34_un_neg_dmi_values, \
dmi_sim_pos, n34_sim_pos, dmi_sim_neg, n34_sim_neg, dmi_todos, n34_todos, \
dmi_pos_n34_neg, dmi_neg_n34_pos = auxScatter(n34, n34_3, dmi, dmi_3, 10)

in_label_size = 18
label_legend_size = 18
tick_label_size = 15
scatter_size_fix = 3

# ---------------------------------------------------------------------------- #
fig, ax = plt.subplots(dpi=dpi, figsize=(10, 10))
# todos
ax.scatter(x=dmi_todos, y=n34_todos, marker='o', label='Niño3.4 vs DMI',
           s=50*scatter_size_fix, edgecolor='k', color='dimgray', alpha=1)
# dmi puros
ax.scatter(x=dmi_un_pos.values, y=dmi_un_pos_n34_values.values, marker='o',
           s=70*scatter_size_fix, edgecolor='k', facecolor='firebrick',
           alpha=1, label='IOD+')
ax.scatter(x=dmi_un_neg.values, y=dmi_un_neg_n34_values.values, marker='o',
           s=70*scatter_size_fix, facecolor='limegreen', edgecolor='k',
           alpha=1, label='IOD-')
# n34 puros
ax.scatter(y=n34_un_pos.values, x=n34_un_pos_dmi_values.values, marker='o',
           s=70*scatter_size_fix, edgecolors='k', facecolor='navy', alpha=1,
           label='El Niño')
ax.scatter(y=n34_un_neg.values, x=n34_un_neg_dmi_values.values, marker='o',
           s=70*scatter_size_fix, edgecolors='k', facecolor='deeppink',
           alpha=1, label='La Niña')
# sim
ax.scatter(x=dmi_sim_pos.values, y=n34_sim_pos.values, marker='o',
           s=50*scatter_size_fix, edgecolor='k', color='#FF5B12', alpha=1,
           label='El Niño & IOD+')
ax.scatter(x=dmi_sim_neg.values, y=n34_sim_neg.values, marker='o',
           s=50*scatter_size_fix, edgecolor='k', color='deepskyblue',
           alpha=1, label='La Niña & IOD-')
# sim opp. sing
ax.scatter(x=dmi_pos_n34_neg.values, y=n34_sim_neg.values, marker='o',
           s=50*scatter_size_fix, edgecolor='k', color='orange', alpha=1,
           label='La Niña & IOD+')
ax.scatter(x=dmi_neg_n34_pos.values, y=n34_sim_pos.values, marker='o',
           s=50*scatter_size_fix, edgecolor='k', color='gold', alpha=1,
           label='El Niño & IOD-')

ax.legend(loc=(.01, .57), fontsize=label_legend_size)
ax.tick_params(axis='both', which='major', labelsize=tick_label_size)
ax.set_ylim((-5, 5))
ax.set_xlim((-5, 5))
ax.axhspan(-.5, .5, alpha=0.2, color='black', zorder=0)
ax.axvspan(-.5, .5, alpha=0.2, color='black', zorder=0)
ax.set_xlabel('DMI', size=in_label_size)
ax.set_ylabel('Niño 3.4', size=in_label_size)
ax.text(-4.9, 4.6, 'EN/IOD-', dict(size=in_label_size))
ax.text(-.2, 4.6, 'EN', dict(size=in_label_size))
ax.text(+3.7, 4.6, 'EN/IOD+', dict(size=in_label_size))
ax.text(+4.2, -.1, 'IOD+', dict(size=in_label_size))
ax.text(+3.78, -4.9, 'LN/IOD+', dict(size=in_label_size))
ax.text(-.2, -4.9, 'LN', dict(size=in_label_size))
ax.text(-4.9, -4.9, ' LN/IOD-', dict(size=in_label_size))
ax.text(-4.9, -.1, 'IOD-', dict(size=in_label_size))
plt.tight_layout()

if save:
    plt.savefig(out_dir + 'figure1.jpg')
else:
    plt.show()

print('#######################################################################')
print('Figure2')
print('#######################################################################')
variables = ['div_UV200', 'vp_from_UV200_w']
t_critic = 1.66 # es MUY similar (2 digitos) para ambos períodos
dmi_or = DMI(filter_bwa=False, start_per='1920', end_per='2020')[2]
n34_or = Nino34CPC(xr.open_dataset(
    "/pikachu/datos/luciano.andrian/verif_2019_2023/sst.mnmean.nc"),
                   start=1920, end=2020)[0]

r_crit = np.sqrt(1 / (((np.sqrt((p[1] - p[0]) - 2) / t_critic) ** 2) + 1))

# indices: ---------------------------------------------------------------------
dmi = dmi_or.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))
n34 = n34_or.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))

# data: ------------------------------------------------------------------------
data_div = xr.open_dataset(data_dir + variables[0] + '.nc')
data_div = data_div.sel(variable='var')
data_div = data_div.rename({'divergence':'var'})
data_div = data_div.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))

data_vp = xr.open_dataset(data_dir + variables[1] + '.nc')
data_vp = data_vp.sel(variable='var')
data_vp = data_vp.rename({'velocity_potential':'var'})
data_vp = data_vp.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))

data_sst = xr.open_dataset(
    "/pikachu/datos/luciano.andrian/verif_2019_2023/sst.mnmean.nc")
data_sst = data_sst.rename({'sst':'var'})
data_sst = data_sst.drop('time_bnds')
data_sst = data_sst.rolling(time=3, center=True).mean()
data_sst = data_sst.sel(time=data_sst.time.dt.month.isin(10))
data_sst = Detrend(data_sst, 'time')
data_sst = data_sst.sel(time=slice('1940-01-01', '2020-12-01'))

time_original = data_sst.time

data_div = data_div - data_div.mean('time')
data_vp = data_vp - data_vp.mean('time')
data_sst = data_sst - data_sst.mean('time')

#------------------------------------------------------------------------------#
div_n34, div_corr_n34, div_dmi, div_corr_dmi, vp_n34, vp_corr_n34, \
    vp_dmi, vp_corr_dmi \
    = ComputeWithEffect(data=data_div, data2=data_vp,
                        n34=n34.sel(time=n34.time.dt.month.isin(10)),
                        dmi=dmi.sel(time=dmi.time.dt.month.isin(10)),
                        two_variables=True, m=10, full_season=False,
                        time_original=time_original)

div_n34_wodmi, div_corr_n34_wodmi, div_dmi_won34, div_corr_dmi_won34 = \
    ComputeWithoutEffect(data_div,
                         n34.sel(time=n34.time.dt.month.isin(10)),
                         dmi.sel(time=dmi.time.dt.month.isin(10)),
                         10, time_original)

vp_n34_wodmi, vp_corr_n34_wodmi, vp_dmi_won34, vp_corr_dmi_won34 = \
    ComputeWithoutEffect(data_vp,
                         n34.sel(time=n34.time.dt.month.isin(10)),
                         dmi.sel(time=dmi.time.dt.month.isin(10)),
                         10, time_original)
#------------------------------------------------------------------------------#
sst_n34, sst_corr_n34, sst_dmi, sst_corr_dmi, aux, aux, aux, aux \
    = ComputeWithEffect(data=data_sst, data2=None,
                        n34=n34.sel(time=n34.time.dt.month.isin(10)),
                        dmi=dmi.sel(time=dmi.time.dt.month.isin(10)),
                        two_variables=False, m=10, full_season=False,
                        time_original=time_original)

sst_n34_wodmi, sst_corr_n34_wodmi, sst_dmi_won34, sst_corr_dmi_won34 \
    = ComputeWithoutEffect(data_sst,
                         n34.sel(time=n34.time.dt.month.isin(10)),
                         dmi.sel(time=dmi.time.dt.month.isin(10)),
                         10, time_original)
#------------------------------------------------------------------------------#
aux_sst = SetDataToPlotFinal(sst_n34 * MakerMaskSig(sst_corr_n34),
                             sst_dmi * MakerMaskSig(sst_corr_dmi),
                             sst_n34_wodmi * MakerMaskSig(sst_corr_n34_wodmi),
                             sst_dmi_won34 * MakerMaskSig(sst_corr_dmi_won34))

aux_vp = SetDataToPlotFinal(vp_n34 * MakerMaskSig(vp_corr_n34),
                            vp_dmi * MakerMaskSig(vp_corr_dmi),
                            vp_n34_wodmi * MakerMaskSig(vp_corr_n34_wodmi),
                            vp_dmi_won34 * MakerMaskSig(vp_corr_dmi_won34))

aux_div = SetDataToPlotFinal(div_n34, div_dmi, div_n34_wodmi, div_dmi_won34)

PlotFinal(data=aux_sst, levels=scale_sst, cmap=cbar_sst,
          titles=subtitulos_regre, namefig='figure2', map='hs',
          save=save, dpi=dpi, out_dir=out_dir,
          data_ctn=aux_vp, levels_ctn=scale_vp, color_ctn='k',
          data_ctn2=aux_div, levels_ctn2=scale_div,
          color_ctn2=['#D300FF', '#00FF5D'])
################################################################################
# HGT y WAF
################################################################################
print('#######################################################################')
print('Figure3')
print('#######################################################################')
variables = ['HGT200_SON_mer_d_w', 'HGT750_SON_mer_d_w']
#seasons = [10]  # main month
plt.rcParams['hatch.linewidth'] = 0.5

dmi_or = DMI(filter_bwa=False, start_per='1920', end_per='2020')[2]
n34_or = Nino34CPC(
    xr.open_dataset(
        "/pikachu/datos/luciano.andrian/verif_2019_2023/sst.mnmean.nc"),
    start=1920, end=2020)[0]

for v, v_count, hpalevel in zip(variables, [0, 1], [200, 750]):
    # indices: --------------------------------------------------------------- #
    dmi = dmi_or.sel(
        time=slice(str(p[0] + y1) + '-01-01', str(p[1]) + '-12-01'))
    n34 = n34_or.sel(
        time=slice(str(p[0] + y1) + '-01-01', str(p[1]) + '-12-01'))

    data = xr.open_dataset(data_dir + v + '.nc')
    data = data.sel(time=slice(str(p[0] + y1) + '-01-01', str(p[1]) + '-12-01'))
    time_original = data.time

    data = data - data.mean('time')

    data_sf = xr.open_dataset(data_dir + 'sf_from_UV' + str(hpalevel) + '_w.nc')
    data_sf = data_sf.sel(variable='var')  # WTF
    data_sf = data_sf.rename({'streamfunction': 'var'})
    data_sf = data_sf.sel(
        time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))

    # Climatology SON -------------------------------------------------------- #
    data_clim = data_sf.sel(time=data_sf.time.dt.month.isin(10)).mean('time')
    data_sf = data_sf.sel(time=data_sf.time.dt.month.isin(10))

    # Anomaly ---------------------------------------------------------------- #
    data_sf = data_sf - data_sf.mean('time')

    hgt_n34, hgt_corr_n34, hgt_dmi, hgt_corr_dmi, sf_n34, sf_corr_n34, \
        sf_dmi, sf_corr_dmi = \
        ComputeWithEffect(data=data, data2=data_sf,
                          n34=n34.sel(time=n34.time.dt.month.isin(10)),
                          dmi=dmi.sel(time=dmi.time.dt.month.isin(10)),
                          two_variables=True, m=10, full_season=False,
                          time_original=time_original)

    px_n34, py_n34 = ComputeWaf(sf_n34, data_clim, hpalevel)
    px_dmi, py_dmi = ComputeWaf(sf_dmi, data_clim, hpalevel)

   # ------------------------------------------------------------------------- #
    hgt_n34_wodmi, hgt_corr_n34_wodmi, hgt_dmi_won34, hgt_corr_dmi_won34 = \
        ComputeWithoutEffect(data,
                             n34.sel(time=n34.time.dt.month.isin(10)),
                             dmi.sel(time=dmi.time.dt.month.isin(10)),
                             10, time_original)

    sf_n34_wodmi, sf_corr_n34_wodmi, sf_dmi_won34, sf_corr_dmi_won34 = \
        ComputeWithoutEffect(data_sf,
                             n34.sel(time=n34.time.dt.month.isin(10)),
                             dmi.sel(time=dmi.time.dt.month.isin(10)),
                             10, time_original)

    px_n34_wodmi, py_n34_wodmi = ComputeWaf(sf_n34_wodmi, data_clim, hpalevel)
    px_dmi_won34, py_dmi_won34 = ComputeWaf(sf_dmi_won34, data_clim, hpalevel)

    # ------------------------------------------------------------------------ #
    aux_hgt = SetDataToPlotFinal(hgt_n34 * MakerMaskSig(hgt_corr_n34),
                                 hgt_dmi * MakerMaskSig(hgt_corr_dmi),
                                 hgt_n34_wodmi *
                                 MakerMaskSig(hgt_corr_n34_wodmi),
                                 hgt_dmi_won34 *
                                 MakerMaskSig(hgt_corr_dmi_won34))

    aux_hgt_wo_corr = SetDataToPlotFinal(hgt_n34, hgt_dmi, hgt_n34_wodmi,
                                         hgt_dmi_won34)

    wafx = SetDataToPlotFinal(px_n34[0, :, :], px_dmi[0, :, :],
                              px_n34_wodmi[0, :, :], px_dmi_won34[0, :, :])

    wafy = SetDataToPlotFinal(py_n34[0, :, :], py_dmi[0, :, :],
                              py_n34_wodmi[0, :, :], py_dmi_won34[0, :, :])

    # ------------------------------------------------------------------------ #
    PlotFinal(data=aux_hgt, levels=scale_hgt, cmap=cbar,
              titles=subtitulos_regre, namefig=f"figure{v_count+3}", map='hs',
              save=save, dpi=dpi, out_dir=out_dir,
              data_ctn=aux_hgt_wo_corr, levels_ctn=None, color_ctn='k',
              data_waf=data_sf, wafx=wafx, wafy=wafy,
              waf_scale=None, waf_label=10e-6, waf_step=4)

################################################################################
# Composites
################################################################################
print('#######################################################################')
print('Figure 5')
print('#######################################################################')
plt.rcParams['hatch.linewidth'] = 1.5


# Vp y Div ------------------------------------------------------------------- #
#----------------------------------------------------------------------------- #
v_from_w = ['div_UV200', 'vp_from_UV200_w'] # creadas a partir de windphsere

data1 = xr.open_dataset(data_dir + v_from_w[0] + '.nc')
data1 = Detrend(
    OrdenarNC_wTime_fromW(data1.rename({'divergence':'var'})), 'time')

data2 = xr.open_dataset(data_dir+ v_from_w[1] + '.nc')
data2 = Detrend(
    OrdenarNC_wTime_fromW(data2.rename({'velocity_potential':'var'})), 'time')

data3 = xr.open_dataset("/pikachu/datos/luciano.andrian/verif_2019_2023/"
                        "sst.mnmean.nc")
data3 = data3.sel(time=data3.time.dt.year.isin(range(1940,2021)))
data3 = data3.rename({'sst':'var'})
data3 = Detrend(data3, 'time')

aux_div = []
aux_vp = []
aux_sst = []
for c_count, c  in enumerate(cases):

    comp1, num_case = CaseComp(data1, s, mmonth=[9, 11],
                               c=c, two_variables=False, data2=None,
                               nc_date_dir=nc_date_dir)

    comp2, num_case = CaseComp(data2, s, mmonth=[9, 11],
                               c=c, two_variables=False, data2=None,
                               nc_date_dir=nc_date_dir)

    comp3, num_case = CaseComp(data3, s, mmonth=[9, 11],
                               c=c, two_variables=False, data2=None,
                               nc_date_dir=nc_date_dir)

    aux_div.append(comp1)
    aux_vp.append(comp2)
    aux_sst.append(comp3)

aux_div = xr.concat(aux_div, dim='plots')

aux_vp = xr.concat(aux_vp, dim='plots')

aux_sst = xr.concat(aux_sst, dim='plots')

PlotFinal(data=aux_sst, levels=scale_sst_comp, cmap=cbar_sst,
          titles=title_case, namefig='figure5', map='hs',
          save=save, dpi=dpi, out_dir=out_dir,
          data_ctn=aux_vp, levels_ctn=scale_vp_comp, color_ctn='k',
          data_ctn2=aux_div, levels_ctn2=scale_div_comp,
          color_ctn2=['#D300FF', '#00FF5D'])

print('#######################################################################')
print('Figure 6-7')
print('#######################################################################')

for v, hpalevel, f_count in zip(['HGT200_SON_mer_d_w', 'HGT750_SON_mer_d_w'],
                                [200,750], [6,7]):

    data = xr.open_dataset(data_dir + v + '.nc')

    data_sf = xr.open_dataset(data_dir + 'sf_from_UV' + str(hpalevel) + '_w.nc')
    data_sf = data_sf.rename({'streamfunction': 'var'})

    aux_hgt = []
    aux_hgt_no_sig = []
    aux_wafx = []
    aux_wafy = []
    for c_count, c in enumerate(cases):
        comp1, num_case, neutro_comp = CaseComp(data, s, mmonth=[9,11],
                                                c=c, two_variables=False,
                                                data2=None,
                                                return_neutro_comp=True,
                                                nc_date_dir=nc_date_dir)

        comp_sf, aux, neutro_comp = CaseComp(data_sf, s, mmonth=[9,11],
                                             c=c, two_variables=False,
                                             data2=None,
                                             return_neutro_comp=True,
                                             nc_date_dir=nc_date_dir)

        px, py = WAF(neutro_comp, comp_sf, comp_sf.lon, comp_sf.lat,
                     reshape=True, variable='var', hpalevel=hpalevel)
        weights = np.transpose(np.tile(-2 * np.cos(
            comp_sf.lat.values * 1 * np.pi / 180) + 2.1, (359, 1)))
        weights_arr = np.zeros_like(px)
        weights_arr[0, :, :] = weights
        px *= weights_arr
        py *= weights_arr

        data_sig = xr.open_dataset(sig_dir + v.split('_')[0] + '_' +
                                   c + '1940_2020_' + s + '.nc')
        sig = comp1.where((comp1 < data_sig['var'][0]) |
                          (comp1 > data_sig['var'][1]))
        sig = sig.where(np.isnan(sig['var']), 1)


        aux_hgt.append(comp1 * sig)
        aux_hgt_no_sig.append(comp1)
        aux_wafx.append(px[0,:,:])
        aux_wafy.append(py[0,:,:])

    aux_hgt = xr.concat(aux_hgt, dim='plots')
    aux_hgt_no_sig = xr.concat(aux_hgt_no_sig, dim='plots')

    aux_wafx = SetDataToPlotFinal(aux_wafx[0], aux_wafx[1], aux_wafx[2],
                                  aux_wafx[3], aux_wafx[4], aux_wafx[5])
    aux_wafy = SetDataToPlotFinal(aux_wafy[0], aux_wafy[1], aux_wafy[2],
                                  aux_wafy[3], aux_wafy[4], aux_wafy[5])

    PlotFinal(data=aux_hgt, levels=scale_hgt_comp, cmap=cbar,
              titles=title_case, namefig=f"figure{f_count}", map='hs',
              save=save, dpi=dpi, out_dir=out_dir,
              data_ctn=aux_hgt_no_sig, levels_ctn=None, color_ctn='k',
              data_waf=data_sf, wafx=aux_wafx, wafy=aux_wafy,
              waf_scale=None, waf_label=10e-5, waf_step=4)
################################################################################
print('#######################################################################')
print('Figure8')
print('#######################################################################')

sd_dmi_s = xr.open_dataset(index_dir + 'DMI_' + s + '_Leads_r_CFSv2.nc').std()
sd_n34_s = xr.open_dataset(index_dir + 'N34_' + s + '_Leads_r_CFSv2.nc').std()

# Sim_pos
c = 'sim_pos'
dmi_sim_pos, n34_sim_pos = SelectParIndex(c, c, sd_dmi_s, sd_n34_s, s)
# sim_neg
c = 'sim_neg'
dmi_sim_neg, n34_sim_neg = SelectParIndex(c, c, sd_dmi_s, sd_n34_s, s)
# dmi_puros_pos
dmi_puros_pos, n34_in_dmi_puros_pos = SelectParIndex('dmi_puros_pos', 'neutros',
                                                     sd_dmi_s, sd_n34_s, s,
                                                     by_r=True, open_dmi=False,
                                                     open_n34=True)
# dmi_puros_neg
dmi_puros_neg, n34_in_dmi_puros_neg = SelectParIndex('dmi_puros_neg', 'neutros',
                                                     sd_dmi_s, sd_n34_s, s,
                                                     by_r=True, open_dmi=False,
                                                     open_n34=True)
# n34_puros_pos
dmi_in_n34_puros_pos, n34_puros_pos = SelectParIndex('neutros', 'n34_puros_pos',
                                                     sd_dmi_s, sd_n34_s, s,
                                                     by_r=True, open_dmi=True,
                                                     open_n34=False)
# n34_puros_neg
dmi_in_n34_puros_neg, n34_puros_neg = SelectParIndex('neutros', 'n34_puros_neg',
                                                     sd_dmi_s, sd_n34_s, s,
                                                     by_r=True, open_dmi=True,
                                                     open_n34=False)
# dmi_pos_n34_neg
dmi_pos_n34_neg, n34_in_dmi_pos_n34_neg = SelectParIndex('dmi_pos_n34_neg',
                                                         'dmi_pos_n34_neg',
                                                         sd_dmi_s, sd_n34_s, s,
                                                         by_r=False,
                                                         open_dmi=False,
                                                         open_n34=False)
# dmi_neg_n34_pos
dmi_neg_n34_pos, n34_in_dmi_neg_n34_pos = SelectParIndex('dmi_neg_n34_pos',
                                                         'dmi_neg_n34_pos',
                                                         sd_dmi_s, sd_n34_s, s,
                                                         by_r=False,
                                                         open_dmi=False,
                                                         open_n34=False)
# neutros
dmi_neutro, n34_neutro = SelectParIndex('neutros', 'neutros',
                                        sd_dmi_s, sd_n34_s, s,
                                        by_r=False)


fig, ax = plt.subplots(dpi=dpi, figsize=(10, 10))

# todos
ax.scatter(y=n34_neutro, x=dmi_neutro, marker='o', label='Niño3.4 vs DMI',
           s=50*scatter_size_fix, edgecolor='k', color='dimgray', alpha=1)

# dmi puros
ax.scatter(x=dmi_puros_pos, y=n34_in_dmi_puros_pos, marker='o',
           s=70*scatter_size_fix, edgecolor='k', color='firebrick',
           alpha=1, label='IOD+')
ax.scatter(x=dmi_puros_neg, y=n34_in_dmi_puros_neg, marker='o',
           s=70*scatter_size_fix, edgecolor='k', color='limegreen',
           alpha=1, label='IOD-')

# n34 puros
ax.scatter(x=dmi_in_n34_puros_pos, y=n34_puros_pos, marker='o',
            s=70*scatter_size_fix, edgecolor='k', color='navy', alpha=1,
           label='El Niño')
ax.scatter(x=dmi_in_n34_puros_neg, y=n34_puros_neg, marker='o',
           s=70*scatter_size_fix, edgecolor='k', color='deeppink',
           alpha=1, label='La Niña')

# sim
ax.scatter(x=dmi_sim_pos, y=n34_sim_pos, marker='o', s=50*scatter_size_fix,
           edgecolor='k', color='#FF5B12', alpha=1, label='El Niño & IOD+')
ax.scatter(x=dmi_sim_neg, y=n34_sim_neg, marker='o', s=50*scatter_size_fix,
           edgecolor='k', color='deepskyblue', alpha=1,
           label='La Niña & IOD-')

# sim opp. sing
ax.scatter(x=dmi_pos_n34_neg, y=n34_in_dmi_pos_n34_neg, marker='o',
           s=50*scatter_size_fix, edgecolor='k', color='orange', alpha=1,
           label='La Niña & IOD+')
ax.scatter(x=dmi_neg_n34_pos, y=n34_in_dmi_neg_n34_pos, marker='o',
            s=50*scatter_size_fix, edgecolor='k', color='gold', alpha=1,
           label='El Niño & IOD-')

ax.legend(loc=(.01, .57), fontsize=label_legend_size)
ax.tick_params(axis='both', which='major', labelsize=tick_label_size)
ax.set_ylim((-5, 5))
ax.set_xlim((-5, 5))
ax.axhspan(-.5, .5, alpha=0.2, color='black', zorder=0)
ax.axvspan(-.5, .5, alpha=0.2, color='black', zorder=0)
ax.set_xlabel('DMI', size=in_label_size)
ax.set_ylabel('Niño 3.4', size=in_label_size)
ax.text(-4.9, 4.6, 'EN/IOD-', dict(size=in_label_size))
ax.text(-.2, 4.6, 'EN', dict(size=in_label_size))
ax.text(+3.7, 4.6, 'EN/IOD+', dict(size=in_label_size))
ax.text(+4.2, -.1, 'IOD+', dict(size=in_label_size))
ax.text(+3.78, -4.9, 'LN/IOD+', dict(size=in_label_size))
ax.text(-.2, -4.9, 'LN', dict(size=in_label_size))
ax.text(-4.9, -4.9, ' LN/IOD-', dict(size=in_label_size))
ax.text(-4.9, -.1, 'IOD-', dict(size=in_label_size))
plt.tight_layout()

if save:
    plt.savefig(out_dir + 'figure8.jpg')
else:
    plt.show()
################################################################################
print('#######################################################################')
print('Figure 9')
print('#######################################################################')
neutro = xr.open_dataset(cases_dir + 'neutros_SON.nc').rename({'sst': 'var'})

aux_sst = []
for c in cases_cfsv2:
    case = xr.open_dataset(cases_dir + c + '_SON.nc').rename({'sst': 'var'})
    #num_case = len(case.time)
    comp = case.mean('time') - neutro.mean('time')
    aux_sst.append(comp)

aux_sst = xr.concat(aux_sst, dim='plots')

PlotFinal(data=aux_sst, levels=scale_sst_comp, cmap=cbar_sst,
          titles=title_case, namefig=f"figure9", map='tr',
          save=save, dpi=dpi, out_dir=out_dir)
################################################################################
print('#######################################################################')
print('Figure 10-11')
print('#######################################################################')
neutro = xr.open_dataset(cases_dir + 'hgt_neutros_SON.nc')\
    .rename({'hgt': 'var'})
neutro = Weights(neutro.__mul__(9.80665))

aux_hgt = []
aux_hgt_snr = []
for c in cases_cfsv2:
    case = xr.open_dataset(cases_dir + 'hgt_' + c + '_SON.nc').\
        rename({'hgt': 'var'})
    case = Weights(case.__mul__(9.80665))
    #num_case = len(case.time)
    comp = case.mean('time') - neutro.mean('time')

    spread = case - comp
    spread = spread.std('time')
    snr = comp / spread

    aux_hgt.append(comp)
    aux_hgt_snr.append(snr)

aux_hgt = xr.concat(aux_hgt, dim='plots')
aux_hgt_snr = xr.concat(aux_hgt_snr, dim='plots')

PlotFinal(data=aux_hgt, levels=scale_hgt_comp, cmap=cbar,
          titles=title_case, namefig=f"figure10", map='hs',
          save=save, dpi=dpi, out_dir=out_dir,
          data_ctn=aux_hgt, color_ctn='k')

PlotFinal(data=aux_hgt_snr, levels=scale_hgt_snr, cmap=cbar_snr,
          titles=title_case, namefig=f"figure11", map='hs',
          save=save, dpi=dpi, out_dir=out_dir,
          data_ctn=aux_hgt_snr, color_ctn='k')
################################################################################
print('#######################################################################')
print('Figure 14')
print('#######################################################################')