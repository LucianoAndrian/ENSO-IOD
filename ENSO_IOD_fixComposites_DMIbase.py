"""
Composites OBSERVADOS
DMI standard
"""
################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os
#from ENSO_IOD_Funciones import OpenDatasets
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
warnings.filterwarnings("ignore")

from ENSO_IOD_Funciones import WAF, CaseComp, PlotComposite_wWAF
################################################################################
nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
              'nc_composites_dates_no_ind_sst_anom/' #fechas
data_dir_t_pp = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                'data_obs_d_w_c/' #T y PP ya procesados
data_dir_era5 = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                'ERA5/1940_2020/' # ERA5 ya procesados
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/paper1/1940_2020/' \
          'composite/dmi_standard/'

#Plot
save = False
dpi = 70
sig = False
waf = True
sig_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_quantiles/' \
          'DMIbase/' # resultados de MC
################################################################################
# colorbars
cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC',
                                 '#B4E2DB',
                                 'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07',
                                 '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')


cbar_t = colors.ListedColormap(['#9B1C00', '#B9391B', '#CD4838', '#E25E55',
                                '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3',
                                '#2064AF', '#014A9B'][::-1])
cbar_t.set_over('#691800')
cbar_t.set_under('#013774')
cbar_t.set_bad(color='white')

cbar_t_r = colors.ListedColormap(['#9B1C00', '#B9391B', '#CD4838', '#E25E55',
                                  '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3',
                                  '#2064AF', '#014A9B'])
cbar_t_r.set_under('#691800')
cbar_t_r.set_over('#013774')
cbar_t_r.set_bad(color='white')

cbar_sst = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89',
                                  '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3',
                                  '#2064AF'][::-1])
cbar_sst.set_over('#9B1C00')
cbar_sst.set_under('#014A9B')
cbar_sst.set_bad(color='white')

cmap_t_pp = [cbar_pp, cbar_t]
cmap_era5 = [cbar_t, cbar_t_r]
################################################################################
# HGT + WAF --------------------------------------------------------------------
plt.rcParams['hatch.linewidth'] = 1.5
tw=[False, False] # por ahora sin vp
sig2 = [True, False]
steps = [1, 1]
contours1 = [True, False]
sig_v = [False, False]

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_un_pos',
         'DMI_un_neg', 'N34_un_pos', 'N34_un_neg']
title_case = ['DMI-ENSO simultaneous positive phase ',
              'DMI-ENSO simultaneous negative phase ',
              'DMI pure positive phase ',
              'DMI pure negative phase ',
              'ENSO pure positive phase ',
              'ENSO pure negative phase ']
seasons = ['SON']
min_max_months = [[9,11]]

v_count = 0
for v, hpalevel in zip(['HGT200_SON_mer_d_w', 'HGT750_SON_mer_d_w'], [200,750]):

    data = xr.open_dataset(data_dir_era5 + v + '.nc')

    if waf:
        data_sf = xr.open_dataset(
            data_dir_era5 + 'sf_from_UV' + str(hpalevel) + '_w.nc')
        data_sf = data_sf.rename({'streamfunction': 'var'})

    c_count = 0
    for c in cases:
        s_count = 0
        for s in seasons:
            comp1, num_case, neutro_comp = \
                CaseComp(data, s, mmonth=min_max_months[s_count], c=c,
                         two_variables=False, data2=None,
                         return_neutro_comp=True,
                         nc_date_dir=nc_date_dir)

            if waf:
                print(v)
                print(hpalevel)
                comp_sf, aux, neutro_comp = \
                    CaseComp(data_sf, s, mmonth=min_max_months[s_count], c=c,
                             two_variables=False, data2=None,
                             return_neutro_comp=True, nc_date_dir=nc_date_dir)


                px, py = WAF(neutro_comp, comp_sf, comp_sf.lon, comp_sf.lat,
                             reshape=True, variable='var', hpalevel=hpalevel)
                weights = np.transpose(
                    np.tile(-2 * np.cos(comp_sf.lat.values * 1 * np.pi / 180) +
                            2.1, (359, 1)))
                weights_arr = np.zeros_like(px)
                weights_arr[0, :, :] = weights
                px *= weights_arr
                py *= weights_arr
            else:
                px, py = None, None
                data_sf = None
                comp_sf = None

            if sig_v[v_count]:
                data_sig = xr.open_dataset(sig_dir + v.split('_')[0] +
                                           '_' + c + '1940_2020_' + s + '.nc')
                #comp1_i = comp1.interp(lon=data_sig.lon.values, lat=data_sig.lat.values)
                sig = comp1.where((comp1 < data_sig['var'][0]) |
                                  (comp1 > data_sig['var'][1]))
                sig = sig.where(np.isnan(sig['var']), 1)

            else:
                data_sig = None
                sig = 1

            PlotComposite_wWAF(comp=comp1 * sig, cmap=cbar_t, step1=1,
                               levels=[-300, -250, -200, -150, -100, -50, -25,
                                       0, 25, 50, 100, 150, 200, 250, 300],
                               contour1=True, two_variables=True, comp2=comp1,
                               linewidht2=1,
                               levels2=[-300, -200, -100, -50, 0, 50, 100, 200,
                                        300], #Levels del contour
                               mapa='hs', significance=False, # No usa los puntos para marcar significancia
                               title=v.split('_')[0] + ' - ' +
                                     title_case[c_count] + '\n' + s +
                                     ' - Events: ' + str(num_case),
                               name_fig=v.split('_')[0] + s + '_' +
                                        cases[c_count] + '_mer_d_w_DMIST',
                               out_dir=out_dir,
                               dpi=dpi, save=save, comp_sig=sig, color_sig='k',
                               color_map='grey',
                               waf=waf, px=px, py=py, data_waf=comp_sf,
                               waf_scale=None, step_waf=4) #WAF

            s_count += 1
        c_count += 1
    v_count += 1

########################################################################################################################
# Vp y Div ------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
v_from_w = ['div_UV200', 'vp_from_UV200_w'] # creadas a partir de windphsere
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

data1 = xr.open_dataset(data_dir_era5 + v_from_w[0] + '.nc')
data1 = Detrend(
    OrdenarNC_wTime_fromW(data1.rename({'divergence':'var'})), 'time')

data2 = xr.open_dataset(data_dir_era5 + v_from_w[1] + '.nc')
data2 = Detrend(
    OrdenarNC_wTime_fromW(data2.rename({'velocity_potential':'var'})), 'time')

data3 = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
data3 = data3.rename({'sst':'var'})
data3 = Detrend(data3, 'time')


c_count = 0
for c in cases:
    s_count = 0
    for s in seasons:
        comp1, num_case = CaseComp(data1, s, mmonth=min_max_months[s_count],
                                   c=c, two_variables=False, data2=None,
                                   nc_date_dir=nc_date_dir)

        comp2, num_case = CaseComp(data2, s, mmonth=min_max_months[s_count],
                                   c=c, two_variables=False, data2=None,
                                   nc_date_dir=nc_date_dir)

        comp3, num_case = CaseComp(data3, s, mmonth=min_max_months[s_count],
                                   c=c, two_variables=False, data2=None,
                                   nc_date_dir=nc_date_dir)

        PlotComposite_wWAF(comp=comp3,
                           levels=[-1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5],
                           cmap=cbar_sst, step1=1, contour1=True,
                           two_variables=True, comp2=comp2,
                           levels2=np.linspace(-4.5e6, 4.5e6, 13),
                           significance=False,
                           mapa='HS',
                           title=
                           'SST, Divergence and potential velocity - 200hPa' +
                           '\n' + title_case[c_count] + '\n' + s +
                           ' - Events: ' + str(num_case),
                           name_fig=
                           'SSTDivVP_' + s + '_' + cases[c_count] + '_d_DMIST',
                           color_map='grey', dpi=dpi, save=save, linewidht2=.8,
                           out_dir=out_dir, third_variable=True, comp3=comp1,
                           levels_contour3=[-1.6e-06, 1.6e-06])
        s_count += 1
    c_count += 1

########################################################################################################################
########################################################################################################################
# T y PP HS y SA
def OpenObsDataSet(name, sa=True,
                   dir='/pikachu/datos/luciano.andrian/observado/ncfiles'
                        '/data_obs_d_w_c/'):

    aux = xr.open_dataset(dir + name + '.nc')
    if sa:
        aux2 = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if len(aux2.lat) > 0:
            return aux2
        else:
            aux2 = aux.sel(lon=slice(270, 330), lat=slice(-60, 15))
            return aux2
    else:
        return aux
########################################################################################################################
# #T y PP con contornos de HGT200
variables_tpp = ['ppgpcc_w_c_d_1', 'tcru_w_c_d_0.25']
title_var = ['PP GPCC', 'T Cru']

variables_tpp = ['ppgpcc_w_c_d_1']
title_var = ['PP GPCC']

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_un_pos',
         'DMI_un_neg', 'N34_un_pos', 'N34_un_neg']
title_case = ['DMI-ENSO simultaneous positive phase ',
              'DMI-ENSO simultaneous negative phase ',
              'DMI pure positive phase ',
              'DMI pure negative phase ',
              'ENSO pure positive phase ',
              'ENSO pure negative phase ']

scales = [np.linspace(-30, 30, 13), # pp
          [-1, -.75, -.5, -.25, -.1, 0, .1, .25, .5, .75, 1]] #t
#seasons = [7, 10] # main month
# min_max_months = [[6,8],[9,11]]
# seasons_name = ['JJA', 'SON']
# SA = borders = [True, False]

min_max_months = [[9,11]]
seasons_name = ['SON']
SA = borders = [True]

aux_name = ['HS', 'SA'] # esto funciona con sa, [True] = 'SA'
aux_name=['hs_ex', 'sa']


for v, v_count in zip(variables_tpp, [0,1]):
    for sa in SA:
        if sa:
            plt.rcParams['hatch.linewidth'] = 2
            hatches = '...'
            color_sig = 'k'
        else:
            plt.rcParams['hatch.linewidth'] = 1
            hatches = '/////'
            color_sig = 'magenta'

        for s, s_count in zip(seasons_name, [0,1]):
            data = OpenObsDataSet(name=v + '_' + s, sa=sa)
            data = data.sel(time=slice('1940-01-01','2020-12-31'))

            # sel ?
            for c, c_count in zip(cases, range(0, len(cases))):
                comp1, num_case = CaseComp(data, s,
                                           mmonth=min_max_months[s_count], c=c,
                                           nc_date_dir=nc_date_dir)


                # Significancia MC
                data_sig = xr.open_dataset(
                    sig_dir + v.split('_')[0] + '_' + c + '1940_2020_' +
                    s + '.nc')
                comp1_i=comp1.interp(lon=data_sig.lon.values,
                                     lat=data_sig.lat.values)
                sig = comp1_i.where((comp1_i < data_sig['var'][0]) |
                                    (comp1_i > data_sig['var'][1]))
                sig = sig.where(np.isnan(sig['var']), 0)

                PlotComposite_wWAF(comp=comp1, levels=scales[v_count],
                                   cmap=cmap_t_pp[v_count], step1=1,
                                   contour1=False,
                                   two_variables=False,
                                   mapa=aux_name[sa], significance=True,
                                   title=title_var[v_count] + '\n' +
                                         title_case[c_count] + '\n' + s +
                                         ' - Events: ' + str(num_case) ,
                                   name_fig=v + '_' + s + '_' +
                                            cases[c_count] + '_mer_d_w_' +
                                            aux_name[sa] + 'AAGG_DMIST',
                                   dpi=dpi, save=save, comp_sig=sig,
                                   color_sig=color_sig, color_map='k',
                                   borders=sa, hatches=hatches,
                                   out_dir=out_dir)
################################################################################