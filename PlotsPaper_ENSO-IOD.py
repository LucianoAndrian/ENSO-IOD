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

from ENSO_IOD_Funciones import (Nino34CPC, DMI, ComputeWithEffect, WAF,
                                ComputeWithoutEffect, SetDataToPlotFinal,
                                PlotFinal)
################################################################################
data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/1940_2020/'
if save:
    dpi = 300
else:
    dpi = 100
# Funciones ####################################################################
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
################################################################################
print('#######################################################################')
print('Figure2')
print('#######################################################################')

variables = ['div_UV200', 'vp_from_UV200_w']
sig = True

scale_vp = [-3e6, -2.5e6, -2e6, -1.5e6, -1e6, -0.5e6, 0, 0.5e6, 1e6, 1.5e6,
            2e6, 2.5e6, 3e6]
scale_sst = [-1, -.5, -.1, 0, .1, .5, 1]
scale_div = [-4.33e-07, 4.33e-07]

cbar = colors.ListedColormap(['#CD4838', '#E25E55', '#F28C89',
                              'white',
                              '#83B9EB', '#5E9AD7', '#3C7DC3'][::-1])
cbar.set_over('#CD4838')
cbar.set_under('#3C7DC3')
cbar.set_bad(color='white')

p = [1940, 2020]

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

s = 'SON'
s_count = 0
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

aux_div = SetDataToPlotFinal(div_n34,# * MakerMaskSig(div_corr_n34),
                             div_dmi,# * MakerMaskSig(div_corr_dmi),
                             div_n34_wodmi,# * MakerMaskSig(div_corr_n34_wodmi),
                             div_dmi_won34)# * MakerMaskSig(div_corr_dmi_won34))

subtitulos = [r"$Ni\tilde{n}o\ 3.4$",  r"$DMI$", r"$Ni\tilde{n}o\ 3.4|_{DMI}$",
              r"$DMI|_{Ni\tilde{n}o\ 3.4}$"]

PlotFinal(data=aux_sst, levels=scale_sst, cmap=cbar,
          titles=subtitulos, namefig='figure2', map='hs',
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
sig = True

# variables = ['hgt200_HS_mer_d_w', 'hgt750_mer_d_w']
variables = ['HGT200_SON_mer_d_w', 'HGT750_SON_mer_d_w']
name_var = ['precip']
title_var = ['HGT200 ERA5', 'HGT750 ERA5']
seasons = [10]  # main month
seasons_name = ['SON']
interp = False
two_variables = False
SA = [False, False]

scale_hgt =   [-150, -100, -75, -50, -25, -15, 0, 15, 25, 50, 75, 100, 150]

cbar = colors.ListedColormap(['#9B1C00', '#B9391B', '#CD4838', '#E25E55',
                              '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3',
                              '#2064AF', '#014A9B'][::-1])
cbar.set_over('#641B00')
cbar.set_under('#012A52')
cbar.set_bad(color='white')

plt.rcParams['hatch.linewidth'] = 0.5

dmi_or = DMI(filter_bwa=False, start_per='1920', end_per='2020')[2]
n34_or = Nino34CPC(
    xr.open_dataset(
        "/pikachu/datos/luciano.andrian/verif_2019_2023/sst.mnmean.nc"),
    start=1920, end=2020)[0]

periodos = [[1940, 2020]]
y1 = 0
p = periodos[0]

for v, v_count, hpalevel in zip(variables, [0, 1], [200, 750]):
    # indices: -----------------------------------------------------------------
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
    # Climatology SON ----------------------------------------------------------
    data_clim = data_sf.sel(time=data_sf.time.dt.month.isin(10)).mean('time')
    data_sf = data_sf.sel(time=data_sf.time.dt.month.isin(10))

    # Anomaly --------------------------------------------------------------
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
                             seasons[s_count], time_original)

    sf_n34_wodmi, sf_corr_n34_wodmi, sf_dmi_won34, sf_corr_dmi_won34 = \
        ComputeWithoutEffect(data_sf,
                             n34.sel(time=n34.time.dt.month.isin(10)),
                             dmi.sel(time=dmi.time.dt.month.isin(10)),
                             seasons[s_count], time_original)

    px_n34_wodmi, py_n34_wodmi = ComputeWaf(sf_n34_wodmi, data_clim, hpalevel)
    px_dmi_won34, py_dmi_won34 = ComputeWaf(sf_dmi_won34, data_clim, hpalevel)

    # ------------------------------------------------------------------------- #
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

    subtitulos = [r"$Ni\tilde{n}o\ 3.4$", r"$DMI$",
                  r"$Ni\tilde{n}o\ 3.4|_{DMI}$",
                  r"$DMI|_{Ni\tilde{n}o\ 3.4}$"]

    PlotFinal(data=aux_hgt, levels=scale_hgt, cmap=cbar,
              titles=subtitulos, namefig=f"figure{v_count+3}", map='hs',
              save=save, dpi=dpi, out_dir=out_dir,
              data_ctn=aux_hgt_wo_corr, levels_ctn=None, color_ctn='k',
              data_waf=data_sf, wafx=wafx, wafy=wafy,
              waf_scale=None, waf_label=10e-6, waf_step=4)

################################################################################
################################################################################









