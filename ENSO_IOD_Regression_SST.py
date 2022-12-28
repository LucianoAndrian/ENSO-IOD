"""
ENSO vs IOD Regression Others Variables
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

save = False
dpi = 150
full_season = False
text = False
# Functions ############################################################################################################
# def Detrend(xrda, dim):
#     aux = xrda.polyfit(dim=dim, deg=1)
#     try:
#         trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
#     except:
#         trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
#     dt = xrda - trend
#     return dt
#
# def OrdenarNC_wTime_fromW(data):
#     newdata = xr.Dataset(
#         data_vars=dict(
#             var=(['time', 'lat', 'lon'], data['var'][0, :, :, :].values)
#
#         ),
#         coords=dict(
#             lon=(['lon'], data.lon),
#             lat=(['lat'], data.lat),
#             time=(['time'], data.time)
#         )
#     )
#     return newdata
#
def SigMask(data, thr):
    return xr.where(np.abs(data) < thr, np.nan, 1)
########################################################################################################################

variables = ['div_from_w', 'vp_from_w']
name_var = ['div_from_w', 'vp_from_w']
#title_var = []
seasons = [7, 10] # main month
seasons_name = ['JJA', 'SON']
interp = [True, False]
two_variables=False
SA = False
sig = True

scale = [-1, -.75, -.5, -.25, -.1, 0, .1, .25, .5, .75, 1]
cbar = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55', '#F28C89',
                              'white',
                              '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar.set_over('#641B00')
cbar.set_under('#012A52')
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
def Detrend(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
    dt = xrda - trend
    return dt

data = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
#data = data.sel(lon=slice(50, 270),lat=slice(20, -20))
data = data.rename({'sst':'var'})
data = data.drop('time_bnds')
data_50_20 = data.sel(time=slice('1950-01-01', '2020-12-01'))

data_sst = Detrend(data_50_20, 'time')
time_original = data_sst.time

# Anomaly -------------------------------------------------------------------------------------------------------------#
# data_div = data_div.groupby('time.month') - data_div.groupby('time.month').mean('time', skipna=True)
# data_vp = data_vp.groupby('time.month') - data_vp.groupby('time.month').mean('time', skipna=True)
data_sst = data_sst.groupby('time.month') - data_sst.groupby('time.month').mean('time', skipna=True)

# 3-month running mean ------------------------------------------------------------------------------------------------#
# data_div = data_div.rolling(time=3, center=True).mean()
# data_vp = data_vp.rolling(time=3, center=True).mean()

s_count=0
for s in seasons_name:
    aux_n34, aux_corr_n34, aux_dmi, \
    aux_corr_dmi, aux_n34_2, aux_corr_n34_2, \
    aux_dmi_2, aux_corr_dmi_2 = ComputeWithEffect(data=data_sst, data2=None, n34=n34, dmi=dmi,
                                                  two_variables=two_variables, m=seasons[s_count],
                                                  full_season=False, time_original=time_original)
    v_count=0
    print('Plot')
    PlotReg(data=aux_n34, data_cor=aux_corr_n34*SigMask(aux_n34, 0.1),
            levels=scale, cmap=cbar, dpi=dpi,
            title='SST' + '\n' + s +
                          ' - ' + str(p[0]) + '-' + str(p[1]) + ' Niño3.4',
            name_fig='vp_div_' + s + str(p[0]) + '_' + str(p[1]) + '_N34',
            save=save, sig=True, sig_point=True, sig2=False, sig_point2=False,
            two_variables=False,
            SA=False, step=1,
            color_map='grey', color_sig='k', r_crit=r_crit, out_dir=out_dir)

    PlotReg(data=aux_dmi, data_cor=aux_corr_dmi*SigMask(aux_dmi, .1),
            levels=scale, cmap=cbar, dpi=dpi,
            title='SST' + '\n' + s +
                          ' - ' + str(p[0]) + '-' + str(p[1]) + ' DMI',
            name_fig='vp_div_' + s + str(p[0]) + '_' + str(p[1]) + '_DMI',
            save=save, sig=True, sig_point=True, sig2=False, sig_point2=False,
            two_variables=False,
            SA=False, step=1,
            color_map='grey', color_sig='k', r_crit=r_crit, out_dir=out_dir)

    aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi = \
        ComputeWithoutEffect(data_sst, n34, dmi, seasons[s_count], time_original)


    PlotReg(data=aux_n34_wodmi, data_cor=aux_corr_n34*SigMask(aux_n34_wodmi, .1),
            levels=scale, cmap=cbar, dpi=dpi,
            title='SST' + '\n' + s +
                          ' - ' + str(p[0]) + '-' + str(p[1]) + ' Niño3.4 -{DMI}',
            name_fig='vp_div_' + s + str(p[0]) + '_' + str(p[1]) + '_N34_wodmi',
            save=save, sig=True, sig_point=True, sig2=False, sig_point2=False,
            two_variables=False,
            SA=False, step=1,
            color_map='grey', color_sig='k', r_crit=r_crit, out_dir=out_dir)

    PlotReg(data=aux_dmi_won34, data_cor=aux_corr_dmi*SigMask(aux_dmi_won34, .1),
            levels=scale, cmap=cbar, dpi=dpi,
            title='SST' + '\n' + s +
                          ' - ' + str(p[0]) + '-' + str(p[1]) + ' DMI -{Niño3.4}',
            name_fig='vp_div_' + s + str(p[0]) + '_' + str(p[1]) + '_DMI_won34',
            save=save, sig=True, sig_point=True, sig2=False, sig_point2=False,
            two_variables=False,
            SA=False, step=1,
            color_map='grey', color_sig='k',  r_crit=r_crit, out_dir=out_dir)

    s_count+=1

########################################################################################################################
