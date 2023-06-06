"""
ENSO vs IOD Regression Others Variables
"""
################################################################################
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
################################################################################
#out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/regression/fix/OV/'
#data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'
data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/1940_2020/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/paper1/1940_2020/' \
          'regression/'

save = True
dpi = 300
full_season = False
text = False
# Functions ####################################################################
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
################################################################################

variables = ['div_UV200', 'vp_from_UV200_w']
#name_var = ['div_from_w', 'vp_from_w']
#title_var = []
seasons = [10] # main month
seasons_name = ['SON']
interp = [True, False]
two_variables=True
SA = False
sig = True

scales = [[-5e-06, -4.33e-07, 0, 4.33e-07, 5e-06],
          [-3e6, -2.5e6, -2e6, -1.5e6, -1e6, -0.5e6, 0, 0.5e6, 1e6, 1.5e6, 2e6,
           2.5e6, 3e6]]

scale_sst = [-1, -.5, -.1, 0, .1, .5, 1]

cbar = colors.ListedColormap(['#CD4838', '#E25E55', '#F28C89',
                              'white',
                              '#83B9EB', '#5E9AD7', '#3C7DC3'][::-1])
cbar.set_over('#CD4838')
cbar.set_under('#3C7DC3')
cbar.set_bad(color='white')

p = [1940, 2020]

t_critic = 1.66 # es MUY similar (2 digitos) para ambos períodos
dmi_or = DMI(filter_bwa=False, start_per='1920', end_per='2020')[2]
n34_or = Nino34CPC(xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc"),
                   start=1920, end=2020)[0]

r_crit = np.sqrt(1 / (((np.sqrt((p[1] - p[0]) - 2) / t_critic) ** 2) + 1))

# indices: ---------------------------------------------------------------------
dmi = dmi_or.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))
n34 = n34_or.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-01'))

# data: ------------------------------------------------------------------------
data_div = xr.open_dataset(data_dir + variables[0] + '.nc')
data_div = data_div.sel(variable='var')
data_div = data_div.rename({'divergence':'var'})
# data_div = Detrend(OrdenarNC_wTime_fromW(data_div.rename({'divergence':'var'})),
#                    'time')
data_div = data_div.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))

data_vp = xr.open_dataset(data_dir + variables[1] + '.nc')
# data_vp = Detrend(
#     OrdenarNC_wTime_fromW(data_vp.rename({'velocity_potential':'var'})), 'time')
data_vp = data_vp.sel(variable='var')
data_vp = data_vp.rename({'velocity_potential':'var'})
data_vp = data_vp.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))

data_sst = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
data_sst = data_sst.rename({'sst':'var'})
data_sst = data_sst.drop('time_bnds')
data_sst = data_sst.rolling(time=3, center=True).mean()
data_sst = data_sst.sel(time=data_sst.time.dt.month.isin(10))
data_sst = Detrend(data_sst, 'time')
data_sst = data_sst.sel(time=slice('1940-01-01', '2020-12-01'))

time_original = data_sst.time

# Anomaly ---------------------------------------------------------------------
# data_div = data_div.groupby('time.month') - \
#            data_div.groupby('time.month').mean('time', skipna=True)
#
# data_vp = data_vp.groupby('time.month') - \
#           data_vp.groupby('time.month').mean('time', skipna=True)
#
# data_sst = data_sst.groupby('time.month') - \
#            data_sst.groupby('time.month').mean('time', skipna=True)
#esto no hace nada
data_div = data_div - data_div.mean('time')
data_vp = data_vp - data_vp.mean('time')
data_sst = data_sst - data_sst.mean('time')

# 3-month running mean ---------------------------------------------------------
# data_div = data_div.rolling(time=3, center=True).mean()
# data_vp = data_vp.rolling(time=3, center=True).mean()
# data_sst = data_sst.rolling(time=3, center=True).mean()

for s, s_count in zip(seasons_name, [0,1]):
    aux_n34, aux_corr_n34, aux_dmi, \
    aux_corr_dmi, aux_n34_2, aux_corr_n34_2, \
    aux_dmi_2, aux_corr_dmi_2 = ComputeWithEffect(
        data=data_div, data2=data_vp,
        n34=n34.sel(time=n34.time.dt.month.isin(10)),
        dmi=dmi.sel(time=dmi.time.dt.month.isin(10)),
        two_variables=two_variables, m=seasons[s_count],
        full_season=False, time_original=time_original)

    aux_n34_sst, aux_corr_n34_sst, aux_dmi_sst, \
    aux_corr_dmi_sst, aux_n34_2_sst, aux_corr_n34_2_sst, \
    aux_dmi_2_sst, aux_corr_dmi_2_sst = ComputeWithEffect(
        data=data_sst, data2=None,
        n34=n34.sel(time=n34.time.dt.month.isin(10)),
        dmi=dmi.sel(time=dmi.time.dt.month.isin(10)),
        two_variables=False, m=seasons[s_count],
        full_season=False, time_original=time_original)

    v_count=0
    #terminar esto
    # nueva mascara para la corr. una q sea posta. graficar solo los sig.
    #solucionar el problema de la div., es muy ruidosa para graficar contorno
    # probar otro step o intercambiar con sst
    print('Plot')
    PlotReg(data=aux_n34_sst*SigDivMask(aux_corr_n34_sst, r_crit),
            data_cor=aux_corr_n34_sst*SigDivMask(aux_corr_n34_sst, r_crit),
            levels=scale_sst, cmap=cbar, dpi=dpi,
            title='SST[shading], Velocity Potential [black cont.], '
                  'divegence [green cont.] - 200hPa' + '\n' + s +
                          ' - ' + str(p[0]) + '-' + str(p[1]) + ' Niño3.4',
            name_fig='vp_div_' + s + str(p[0]) + '_' + str(p[1]) + '_N34',
            save=save, sig=False, sig_point=True, sig2=False, sig_point2=False,
            two_variables=True,
            SA=False, step=1,
            color_map='grey', color_sig='k', color_sig2='grey',
            data2=aux_n34_2, data_cor2=aux_corr_n34_2, levels2=scales[1],
            r_crit=r_crit, out_dir=out_dir,
            third_variable=True, data3=aux_n34, levels3=[-4.33e-07, 4.33e-07])

    PlotReg(data=aux_dmi_sst*SigDivMask(aux_corr_dmi_sst, r_crit),
            data_cor=aux_corr_dmi*SigDivMask(aux_dmi, scales[0][3]),
            levels=scale_sst, cmap=cbar, dpi=dpi,
            title='SST[shading], Velocity Potential [black cont.], '
                  'divegence [green cont.] - 200hPa' + '\n' + s +
                          ' - ' + str(p[0]) + '-' + str(p[1]) + ' DMI',
            name_fig='vp_div_' + s + str(p[0]) + '_' + str(p[1]) + '_DMI',
            save=save, sig=False, sig_point=True, sig2=False, sig_point2=False,
            two_variables=True,
            SA=False, step=1,
            color_map='grey', color_sig='k', color_sig2='grey',
            data2=aux_dmi_2, data_cor2=aux_corr_dmi_2, levels2=scales[1],
            r_crit=r_crit, out_dir=out_dir,
            third_variable=True, data3=aux_dmi, levels3=[-4.33e-07, 4.33e-07])

    aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi = \
        ComputeWithoutEffect(data_div,
                             n34.sel(time=n34.time.dt.month.isin(10)),
                             dmi.sel(time=dmi.time.dt.month.isin(10)),
                             seasons[s_count], time_original)
    aux_n34_wodmi_2, aux_corr_n34_2, aux_dmi_won34_2, aux_corr_dmi_2 = \
        ComputeWithoutEffect(data_vp,
                             n34.sel(time=n34.time.dt.month.isin(10)),
                             dmi.sel(time=dmi.time.dt.month.isin(10)),
                             seasons[s_count], time_original)
    aux_n34_wodmi_sst, aux_corr_n34_sst, aux_dmi_won34_sst, aux_corr_dmi_sst = \
        ComputeWithoutEffect(data_sst,
                             n34.sel(time=n34.time.dt.month.isin(10)),
                             dmi.sel(time=dmi.time.dt.month.isin(10)),
                             seasons[s_count], time_original)


    PlotReg(data=aux_n34_wodmi_sst*SigDivMask(aux_corr_n34_sst, r_crit),
            data_cor=aux_corr_n34_sst*SigDivMask(aux_n34_wodmi_sst, r_crit),
            levels=scale_sst, cmap=cbar, dpi=dpi,
            title='SST[shading], Velocity Potential [black cont.], '
                  'divegence [green cont.] - 200hPa' + '\n' + s +
                          ' - ' + str(p[0]) + '-' + str(p[1]) + ' Niño3.4 -{DMI}',
            name_fig='vp_div_' + s + str(p[0]) + '_' + str(p[1]) + '_N34_wodmi',
            save=save, sig=False, sig_point=True, sig2=False, sig_point2=False,
            two_variables=True,
            SA=False, step=1,
            color_map='grey', color_sig='k', color_sig2='grey',
            data2=aux_n34_wodmi_2, data_cor2=aux_corr_n34_2, levels2=scales[1],
            r_crit=r_crit, out_dir=out_dir,
            third_variable=True, data3=aux_n34_wodmi, levels3=[-4.33e-07, 4.33e-07])

    PlotReg(data=aux_dmi_won34_sst*SigDivMask(aux_corr_dmi_sst, r_crit),
            data_cor=aux_corr_dmi*SigDivMask(aux_dmi_won34, scales[0][3]),
            levels=scale_sst, cmap=cbar, dpi=dpi,
            title='SST[shading], Velocity Potential [black cont.], '
                  'divegence [green cont.] - 200hPa' + '\n' + s +
                          ' - ' + str(p[0]) + '-' + str(p[1]) + ' DMI -{Niño3.4}',
            name_fig='vp_div_' + s + str(p[0]) + '_' + str(p[1]) + '_DMI_won34',
            save=save, sig=False, sig_point=True, sig2=False, sig_point2=False,
            two_variables=True,
            SA=False, step=1,
            color_map='grey', color_sig='k', color_sig2='grey',
            data2=aux_dmi_won34_2, data_cor2=aux_corr_dmi_2, levels2=scales[1],
            r_crit=r_crit, out_dir=out_dir,
            third_variable=True, data3=aux_dmi_won34, levels3=[-4.33e-07, 4.33e-07])

#
# ################################################################################
# #--------------------------- RWS ----------------------------------------------#
# ################################################################################
# # data: ------------------------------------------------------------------------
# #data_rws = xr.open_dataset(data_dir + 'RWS.nc')
# data_rws = xr.open_dataset(data_dir + 'RWSUV200.nc')
# data_rws = Detrend(OrdenarNC_wTime_fromW(data_rws),'time')
# data_rws = data_rws.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))
# time_original = data_rws.time
#
# # Anomaly ----------------------------------------------------------------------
# data_rws = data_rws.groupby('time.month') - \
#            data_rws.groupby('time.month').mean('time', skipna=True)
#
# # 3-month running mean ---------------------------------------------------------
# data_rws = data_rws.rolling(time=3, center=True).mean()
#
# cbar_rws = colors.ListedColormap(['#4FCDF4','white', '#F4C847'])
# cbar_rws.set_over('#F49927')
# cbar_rws.set_under('#3489F4')
# cbar_rws.set_bad(color='white')
#
# s_count=0
# for s in seasons_name:
#     aux_n34, aux_corr_n34, aux_dmi, \
#     aux_corr_dmi, aux_n34_2, aux_corr_n34_2, \
#     aux_dmi_2, aux_corr_dmi_2 = ComputeWithEffect(
#         data=data_rws, data2=None, n34=n34, dmi=dmi, two_variables=False,
#         m=seasons[s_count], full_season=False, time_original=time_original)
#
#     v_count=0
#     print('Plot')
#     PlotReg(data=aux_n34, data_cor=aux_corr_n34,
#             levels=[-1.2e-10, -0.1e-10, 0, 0.1e-10,1.2e-10], cmap=cbar_rws,
#             dpi=dpi,
#             title='RWS ' + s + ' - ' + str(p[0]) + '-' + str(p[1]) + ' Niño3.4',
#             name_fig='rws_' + s + str(p[0]) + '_' + str(p[1]) + '_N34',
#             save=save, sig=True, sig_point=True, sig2=True, sig_point2=False,
#             two_variables=False,
#             SA=False, step=2,
#             color_map='grey', color_sig='k', r_crit=r_crit, out_dir=out_dir)
#
#     PlotReg(data=aux_dmi, data_cor=aux_corr_dmi,
#             levels=[-1.2e-10, -0.1e-10, 0, 0.1e-10,1.2e-10], cmap=cbar_rws,
#             dpi=dpi,
#             title='RWS ' + s + ' - ' + str(p[0]) + '-' + str(p[1]) + ' DMI',
#             name_fig='rws_' + s + str(p[0]) + '_' + str(p[1]) + '_DMI',
#             save=save, sig=True, sig_point=True, sig2=True, sig_point2=False,
#             two_variables=False,
#             SA=False, step=2,
#             color_map='grey', color_sig='k', r_crit=r_crit, out_dir=out_dir)
#
#     aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi = \
#         ComputeWithoutEffect(data_rws, n34, dmi, seasons[s_count], time_original)
#     # aux_n34_wodmi_2, aux_corr_n34_2, aux_dmi_won34_2, aux_corr_dmi_2 = \
#     #     ComputeWithoutEffect(data_vp, n34, dmi, seasons[s_count], time_original)
#
#
#     PlotReg(data=aux_n34_wodmi, data_cor=aux_corr_n34,
#             levels=[-1.2e-10, -0.1e-10, 0, 0.1e-10,1.2e-10], cmap=cbar_rws,
#             dpi=dpi,
#             title='RWS ' + s + ' - ' + str(p[0]) + '-' + str(p[1]) + ' Niño3.4 -{DMI}',
#             name_fig='rws_' + s + str(p[0]) + '_' + str(p[1]) + '_N34_wodmi',
#             save=save, sig=True, sig_point=True, sig2=True, sig_point2=False,
#             two_variables=False,
#             SA=False, step=2,
#             color_map='grey', color_sig='k', r_crit=r_crit, out_dir=out_dir)
#
#     PlotReg(data=aux_dmi_won34, data_cor=aux_corr_dmi,
#             levels=[-1.2e-10, -0.1e-10, 0, 0.1e-10,1.2e-10], cmap=cbar_rws,
#             dpi=dpi,
#             title='RWS ' + s + ' - ' + str(p[0]) + '-' + str(p[1]) + ' DMI -{Niño3.4}',
#             name_fig='rws_' + s + str(p[0]) + '_' + str(p[1]) + '_DMI_won34',
#             save=save, sig=True, sig_point=True, sig2=True, sig_point2=False,
#             two_variables=False,
#             SA=False, step=2,
#             color_map='grey', color_sig='k', r_crit=r_crit, out_dir=out_dir)
#
#     s_count+=1


########################################################################################################################
# stream function e intento WAF
########################################################################################################################
# # data: ---------------------------------------------------------------------------------------------------------------#
# cbar = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
#                               'white',
#                               '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
# cbar.set_over('#641B00')
# cbar.set_under('#012A52')
# cbar.set_bad(color='white')
#
#
# data_sf = xr.open_dataset(data_dir + 'sf.nc')
# data_sf = Detrend(OrdenarNC_wTime_fromW(data_sf.rename({'streamfunction':'var'})), 'time')
# data_sf = data_sf.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))
# time_original = data_sf.time
#
# # Anomaly -------------------------------------------------------------------------------------------------------------#
# data_sf = data_sf.groupby('time.month') - data_sf.groupby('time.month').mean('time', skipna=True)
#
# # 3-month running mean ------------------------------------------------------------------------------------------------#
# data_sf = data_sf.rolling(time=3, center=True).mean()
# #data_sf = data_sf.sel(variable='var').drop('variable')
# for s, s_count in zip(seasons_name, [0,1]):
#     aux_n34, aux_corr_n34, aux_dmi, \
#     aux_corr_dmi, aux_n34_2, aux_corr_n34_2, \
#     aux_dmi_2, aux_corr_dmi_2 = ComputeWithEffect(data=data_sf, data2=None, n34=n34, dmi=dmi,
#                                                   two_variables=False, m=seasons[s_count],
#                                                   full_season=False, time_original=time_original)
#
#     v_count=0
#     print('Plot')
#     PlotReg(data=aux_n34*MakerMaskSig(aux_corr_n34), data_cor=None,
#             levels= [-3e6, -2.5e6, -2e6, -1.5e6, -1e6, -0.5e6, -0.25e6, 0, 0.25e6, 0.5e6, 1e6, 1.5e6, 2e6, 2.5e6, 3e6],
#             cmap=cbar, dpi=dpi,
#             title='sf200' + '\n' + s + ' - ' + str(p[0]) + '-' + str(p[1]) + ' Niño3.4',
#             name_fig='sf_' + s + str(p[0]) + '_' + str(p[1]) + '_N34',
#             save=save, sig=False, sig_point=False, sig2=False, sig_point2=False,
#             two_variables=True,
#             SA=False, step=1,
#             color_map='grey', color_sig='k', color_sig2='grey',
#             data2=aux_n34, data_cor2=None,
#             levels2=[-3e6, -2e6, -1e6, -0.5e6, -0.25e6, 0, 0.25e6, 0.5e6, 1e6, 2e6, 3e6],
#             r_crit=r_crit, out_dir=out_dir)
#
#     PlotReg(data=aux_dmi*MakerMaskSig(aux_corr_dmi), data_cor=None,
#             levels= [-3e6, -2.5e6, -2e6, -1.5e6, -1e6, -0.5e6, -0.25e6, 0, 0.25e6, 0.5e6, 1e6, 1.5e6, 2e6, 2.5e6, 3e6],
#             cmap=cbar, dpi=dpi,
#             title='sf200' + '\n' + s + ' - ' + str(p[0]) + '-' + str(p[1]) + ' DMI',
#             name_fig='sf_' + s + str(p[0]) + '_' + str(p[1]) + '_DMI',
#             save=save, sig=False, sig_point=False, sig2=False, sig_point2=False,
#             two_variables=True,
#             SA=False, step=1,
#             color_map='grey', color_sig='k', color_sig2='grey',
#             data2=aux_dmi, data_cor2=None,
#             levels2=[-3e6, -2e6, -1e6, -0.5e6, -0.25e6, 0, 0.25e6, 0.5e6, 1e6, 2e6, 3e6],
#             r_crit=r_crit, out_dir=out_dir)
#
#     aux_n34_wodmi, aux_corr_n34, aux_dmi_won34, aux_corr_dmi = \
#         ComputeWithoutEffect(data_sf, n34, dmi, seasons[s_count], time_original)
#
#     PlotReg(data=aux_n34_wodmi*MakerMaskSig(aux_corr_n34), data_cor=None,
#             levels= [-3e6, -2.5e6, -2e6, -1.5e6, -1e6, -0.5e6, -0.25e6, 0, 0.25e6, 0.5e6, 1e6, 1.5e6, 2e6, 2.5e6, 3e6],
#             cmap=cbar, dpi=dpi,
#             title='sf200' + '\n' + s + ' - ' + str(p[0]) + '-' + str(p[1]) +  ' Niño3.4 -{DMI}',
#             name_fig='sf_' + s + str(p[0]) + '_' + str(p[1]) + '_N34_wodmi',
#             save=save, sig=False, sig_point=False, sig2=False, sig_point2=False,
#             two_variables=True,
#             SA=False, step=1,
#             color_map='grey', color_sig='k', color_sig2='grey',
#             data2=aux_n34_wodmi, data_cor2=None,
#             levels2=[-3e6, -2e6, -1e6, -0.5e6, -0.25e6, 0, 0.25e6, 0.5e6, 1e6, 2e6, 3e6],
#             r_crit=r_crit, out_dir=out_dir)
#
#     PlotReg(data=aux_dmi_won34*MakerMaskSig(aux_corr_dmi), data_cor=None,
#             levels= [-3e6, -2.5e6, -2e6, -1.5e6, -1e6, -0.5e6, -0.25e6, 0, 0.25e6, 0.5e6, 1e6, 1.5e6, 2e6, 2.5e6, 3e6],
#             cmap=cbar, dpi=dpi,
#             title='sf200' + '\n' + s + ' - ' + str(p[0]) + '-' + str(p[1]) +  ' '
#                                                                               'DMI -{Noño3.4}',
#             name_fig='sf_' + s + str(p[0]) + '_' + str(p[1]) + '_DMI_won34',
#             save=save, sig=False, sig_point=False, sig2=False, sig_point2=False,
#             two_variables=True,
#             SA=False, step=1,
#             color_map='grey', color_sig='k', color_sig2='grey',
#             data2=aux_dmi_won34, data_cor2=None,
#             levels2=[-3e6, -2e6, -1e6, -0.5e6, -0.25e6, 0, 0.25e6, 0.5e6, 1e6, 2e6, 3e6],
#             r_crit=r_crit, out_dir=out_dir)

#
# from ENSO_IOD_Funciones import WAF
# data_sf = xr.open_dataset(data_dir + 'sf.nc')
# data_sf = Detrend(OrdenarNC_wTime_fromW(data_sf.rename({'streamfunction':'var'})), 'time')
# data_sf = data_sf.sel(time=slice(str(p[0]) + '-01-01', str(p[1]) + '-12-31'))
# time_original = data_sf.time
# data_clim = data_sf.sel(time=data_sf.time.dt.month.isin(10)).mean('time')
#
# aux_xr_dmi_won34 = xr.Dataset(
#         data_vars=dict(
#             var=(['lat', 'lon'], aux_dmi_won34)
#
#         ),
#         coords=dict(
#             lon=(['lon'], aux_dmi_won34.lon),
#             lat=(['lat'], aux_dmi_won34.lat)
#         )
#     )
# px, py = WAF(data_clim, aux_xr_dmi_won34, data_clim.lon, data_clim.lat, reshape=True, variable='var')
# px_xr = xr.DataArray(px[0,:,:])
# py_xr = xr.DataArray(py[0,:,:])
# from ENSO_IOD_Funciones import PlotWAFCountours
#
# weights = np.transpose(np.tile(-2 * np.cos(np.arange(-90, 89) * 1 * np.pi / 180) + 2.1, (359, 1)))
# weights_arr = np.zeros_like(px)
# weights_arr[0, :, :] = weights
#
# PlotWAFCountours(aux_xr_dmi_won34, aux_xr_dmi_won34['var'], title='title', name_fig='name_fig',
#                  save=save, dpi=dpi,
#                  levels=[-3e6, -2.5e6, -2e6, -1.5e6, -1e6, -0.5e6, -0.25e6, 0, 0.25e6, 0.5e6, 1e6, 1.5e6, 2e6, 2.5e6, 3e6],
#                  contour=True, cmap=cbar, number_events=None,
#                  waf=True, px=px * weights_arr, py=py * weights_arr, text=False, waf_scale=1/4000,
#                  two_variables=False, step_waf=5)
