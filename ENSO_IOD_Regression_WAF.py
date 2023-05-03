"""
ENSO vs IOD Regression. intento WAF
"""
########################################################################################################################
from itertools import groupby
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
import xarray as xr
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import os
from matplotlib import colors
from ENSO_IOD_Funciones import Nino34CPC
from ENSO_IOD_Funciones import DMI

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

from ENSO_IOD_Funciones import ComputeWithEffect, ComputeWithoutEffect, PlotReg
########################################################################################################################
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/regression/fix/OV/'
data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'

save = True
dpi = 100
full_season = False
text = False
# Functions ############################################################################################################
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


def CompositeOfWAF_EDIT(original_data, index, mmin, mmax, neutro_comp):
    def is_months(month, mmin, mmax):
        return (month >= mmin) & (month <= mmax)

    if len(index) != 0:
        comp_field = original_data.sel(time=original_data.time.dt.year.isin([index]))
        comp_field = comp_field.sel(
            time=is_months(month=comp_field['time.month'], mmin=mmin, mmax=mmax))

        #WAF
        comp_field = comp_field.groupby('time.year').mean()
        comp_field = comp_field - neutro_comp

        y_count=0
        for y in comp_field.year.values:
            comp = comp_field.sel(year=y)

            px, py = WAF(neutro_comp, comp, neutro_comp.lon, neutro_comp.lat, reshape=True, variable='var')

            if y_count==0:
                px_xr = xr.DataArray(px)
                py_xr = xr.DataArray(py)
            else:
                px_xr_next = xr.DataArray(px)
                py_xr_next = xr.DataArray(py)
                px_xr = xr.concat([px_xr, px_xr_next], dim='dim_0')
                py_xr = xr.concat([py_xr, py_xr_next], dim='dim_0')

            del px, py
            y_count += 1


        #composite:
        px_final = px_xr
        #px_final = np.expand_dims(px_final, axis=0) #para graficar
        py_final = py_xr
        #py_final = np.expand_dims(py_final, axis=0)  # para graficar

        return px_final, py_final




def ReshapeData(data):
    data = xr.Dataset(
        data_vars=dict(
            var=(['y', 'x'], data['var'][0, :, :].values)

        ),
        coords=dict(
            lon=(['x'], data.lon),
            lat=(['y'], data.lat),
        )
    )
    return data
########################################################################################################################

# variables = ['div_from_w', 'vp_from_w']
# name_var = ['div_from_w', 'vp_from_w']
#title_var = []
seasons = [7, 10] # main month
seasons_name = ['JJA', 'SON']
interp = False
two_variables=False
SA = False
sig = False

scales = [[-5e-06, -4.33e-07, 0, 4.33e-07, 5e-06], [-3e6, -2.5e6, -2e6, -1.5e6, -1e6, -0.5e6, 0, 0.5e6, 1e6, 1.5e6, 2e6, 2.5e6, 3e6]]
scale_sst = [-1, -.5, -.1, 0, .1, .5, 1]

cbar = colors.ListedColormap(['#CD4838', '#E25E55', '#F28C89',
                              'white',
                              '#83B9EB', '#5E9AD7', '#3C7DC3'][::-1])
cbar.set_over('#CD4838')
cbar.set_under('#3C7DC3')
cbar.set_bad(color='white')

p = [1950, 2020]





t_critic = 1.66 # es MUY similar (2 digitos) para ambos períodos
dmi_or = DMI(filter_bwa=False, start_per='1920', end_per='2020')[2]
n34_or = Nino34CPC(xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc"),
                   start=1920, end=2020)[0]

r_crit = np.sqrt(1 / (((np.sqrt((p[1] - p[0]) - 2) / t_critic) ** 2) + 1))

# indices: ------------------------------------------------------------------------------------------------------------#
dmi = dmi_or.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))
n34 = n34_or.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))

# data: ---------------------------------------------------------------------------------------------------------------#
data = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/sf.nc')
data = data.rename({'streamfunction':'var'})

# Climatology ---------------------------------------------------------------------------------------------------------#
data_clim = data.groupby('time.month').mean('time', skipna=True)
data_clim = data_clim.rolling(month=3, center=True).mean() #esto es malo en los bordes, enero y diciembre, pero no los uso
# Anomaly -------------------------------------------------------------------------------------------------------------#
data = data.groupby('time.month') - data.groupby('time.month').mean('time', skipna=True)

# 3-month running mean ------------------------------------------------------------------------------------------------#
data = data.rolling(time=3, center=True).mean()

# WAF -----------------------------------------------------------------------------------------------------------------#
from ENSO_IOD_Funciones import WAF

dt=data.time.values[1]
dt_seasons = data.time.values
months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] * 71

for dt, dt_count in zip(dt_seasons, range(0, len(dt_seasons))):
    data_sel = data.sel(time=dt)

    data_sel = ReshapeData(data_sel)
    data_clim_season = ReshapeData(data_clim.sel(month=months[dt_count]))

    px, py = WAF(data_clim_season, data_sel, data_clim.lon, data_clim.lat, reshape=True, variable='var')

    if dt_count == 0:
        px_xr = xr.DataArray(px)
        py_xr = xr.DataArray(py)
    else:
        px_xr_next = xr.DataArray(px)
        py_xr_next = xr.DataArray(py)
        px_xr = xr.concat([px_xr, px_xr_next], dim='dim_0')
        py_xr = xr.concat([py_xr, py_xr_next], dim='dim_0')

    del px, py

# <xarray.Dataset>
# Dimensions:  (lat: 89, lon: 180, time: 852)
# Coordinates:
#   * lat      (lat) float32 88.0 86.0 84.0 82.0 80.0 ... -82.0 -84.0 -86.0 -88.0
#   * lon      (lon) float32 0.0 2.0 4.0 6.0 8.0 ... 350.0 352.0 354.0 356.0 358.0
#   * time     (time) datetime64[ns] 1950-01-01 1950-02-01 ... 2020-12-01
# Data variables:
#     var      (time, lat, lon) float64 nan nan nan nan nan ... nan nan nan nan

px_xr_f = xr.Dataset(
        data_vars=dict(
            var=(['time','lat', 'lon'], px_xr)

        ),
        coords=dict(
            lon=(['lon'], data_clim_season.lon),
            lat=(['lat'], data_clim_season.lat),
            time=(['time'], pd.DatetimeIndex(dt_seasons))
        )
    )
px_xr_f.coords['month'] = ('time', data.month)

py_xr_f = xr.Dataset(
        data_vars=dict(
            var=(['time','lat', 'lon'], py_xr)

        ),
        coords=dict(
            lon=(['lon'], data_clim_season.lon),
            lat=(['lat'], data_clim_season.lat),
            time=(['time'], dt_seasons)
        )
    )
py_xr_f.coords['month'] = ('time', data.month)

x_aux_n34, x_aux_corr_n34, x_aux_dmi, \
x_aux_corr_dmi, x_aux_n34_2, x_aux_corr_n34_2, \
x_aux_dmi_2, x_aux_corr_dmi_2 = ComputeWithEffect(data=px_xr_f, data2=None, n34=n34, dmi=dmi,
                                              two_variables=False, m=10,
                                              full_season=False, time_original=dmi.time)

y_aux_n34, y_aux_corr_n34, y_aux_dmi, \
y_aux_corr_dmi, y_aux_n34_2, y_aux_corr_n34_2, \
y_aux_dmi_2, y_aux_corr_dmi_2 = ComputeWithEffect(data=py_xr_f, data2=None, n34=n34, dmi=dmi,
                                              two_variables=False, m=10,
                                              full_season=False, time_original=dmi.time)

from ENSO_IOD_Funciones import PlotWAFCountours

weights = np.transpose(np.tile(-2 * np.cos(np.arange(-90, 89) * 1 * np.pi / 180) + 2.1, (359, 1)))
weights_arr = np.zeros_like(px)
weights_arr[0, :, :] = weights

PlotWAFCountours(x_aux_dmi, x_aux_dmi, title="title", name_fig="name_fig",
                 save=False, dpi=100, levels=np.linspace(-1.2e-10, 1.2e-10, 4),
                 contour=False, cmap='RdBu', number_events=0,
                 waf=True, px=x_aux_n34 * weights_arr[0,:,:], py=y_aux_n34 * weights_arr[0,:,:], text=False, waf_scale=1,
                 two_variables=False, comp2=None, step=1, step_waf=5,
                 levels2=np.linspace(-450, 450, 13), contour0=False)

fig = plt.figure(figsize=(9, 3.5), dpi=dpi)
from numpy import ma
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([0, 359, -80, 20], crs=crs_latlon)
Q60 = np.percentile(np.sqrt(np.add(np.power(px, 2), np.power(py, 2))), 0)
M = np.sqrt(np.add(np.power(px, 2), np.power(py, 2))) < Q60
# mask array
px_mask = ma.array(px, mask=M)
py_mask = ma.array(py, mask=M)
# plot vectors
lons, lats = np.meshgrid(px_xr_f.lon.values, px_xr_f.lat.values)
ax.quiver(lons[::step_waf, ::step_waf],
          lats[::step_waf, ::step_waf],
          px_mask[::step_waf, ::step_waf],
          py_mask[::step_waf, ::step_waf], transform=crs_latlon, pivot='tail',
          width=0.0020, headwidth=4., alpha=1, color=color_arrow, scale=0.001, scale_units=None)
# , scale=1/10)#, width=1.5e-3, headwidth=3.1,  # headwidht (default3)
# headlength=2.2)  # (default5))
plt.show()