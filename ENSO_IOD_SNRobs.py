"""
SNR de ENSO, IOD y ENSO-IOD observados
Calculados de manera similar al del CFSv2
"""
################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
warnings.filterwarnings("ignore")

from ENSO_IOD_Funciones import CaseSNR, PlotComposite_wWAF, MakeMask
################################################################################
DMIbase = True # si es False usa DMI true-dipole

# Plot
save = True
save_nc = True
dpi = 300
plot = False
################################################################################
# rutas segun criterio del DMI
if DMIbase:
    # fechas
    nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                  'nc_composites_dates_no_ind_sst_anom/'

    # salidas
    out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/snr/dmi_standard/'
    out_dir2 = '/pikachu/datos/luciano.andrian/esquemas/'

    aux_name = 'DMI_standard'
else:
    # fechas
    nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                  'nc_composites_dates/'
    # salidas
    out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/snr/dmi_true_dipole/'
    print('set save_nc = False')
    save_nc = False
    aux_name = 'DMI_true_dipole'

# T y PP
data_dir_t_pp = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                'data_obs_d_w_c/'
# ERA5
data_dir_era5 = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                'ERA5/1940_2020/'

################################################################################
def OpenObsDataSet(name, sa=True, dir=data_dir_t_pp):

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
################################################################################
# colorbar
cbar_snr_hgt = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#6FFE9B',
                                  '#FFFFFF',
                                  '#FFFFFF', '#FFFFFF',
                                  '#FEB77E','#CA3E72','#782281','#251255'])
cbar_snr_hgt.set_over('#251255')
cbar_snr_hgt.set_under('#070B4F')
cbar_snr_hgt.set_bad(color='white')

cbar_snr = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#6FFE9B',
                                  '#FFFFFF',
                                  '#FFFFFF', '#FFFFFF',
                                  '#FEB77E','#CA3E72','#782281','#251255'][::-1])
cbar_snr.set_under('#251255')
cbar_snr.set_over('#070B4F')
cbar_snr.set_bad(color='white')

cbar_snr_t = colors.ListedColormap(['#003069','#4142B0', '#67B3FE' ,'#A9FED6',
                                  '#FFFFFF',
                                  '#FFFFFF', '#FFFFFF',
                                  '#FEF785','#FE8051','#B22E3C','#710733'])
cbar_snr_t.set_under('#003069')
cbar_snr_t.set_over('#710733')
cbar_snr_t.set_bad(color='white')

cbars = [cbar_snr, cbar_snr_t]
################################################################################
# general
cases = ['N34_pos', 'N34_neg', 'DMI_pos', 'DMI_neg',
         'DMI_sim_pos', 'DMI_sim_neg', 'DMI_un_pos',
         'DMI_un_neg', 'N34_un_pos', 'N34_un_neg']

title_case = ['ENSO positive phase',
              'ENSO negative phase',
              'IOD positive phase',
              'IOD negative phase',
              'DMI-ENSO simultaneous positive phase ',
              'DMI-ENSO simultaneous negative phase ',
              'DMI pure positive phase ',
              'DMI pure negative phase ',
              'ENSO pure positive phase ',
              'ENSO pure negative phase ']

seasons = ['JJA','SON']
min_max_months = [[6,8], [9,11]]

if DMIbase:
    print('DMI Standard')
else:
    print('DMI true-dipole')

print('-------------------------------------------------')
if save:
    print('Plots will be saved with dpi = ' + str(dpi))
else:
    print('Plots with dpi =  ' + str(dpi))

print('-------------------------------------------------')
print('Computing and plot SNR for: ')

################################################################################
print('HGT 200 y 750hPa')
# HGT
# seasons
for s, s_count in zip(seasons, range(0, len(seasons))):
    # esto no es muy necesario, se podria abrir un solo archivo
    # la selección de fechas se hace en CaseSNR junto con case

    # variable (level)
    for v, hpalevel, v_count in zip(['HGT200_' + s + '_mer_d_w',
                                     'HGT750_' + s + '_mer_d_w'],
                                    [200,750], [0,1]):

        data = xr.open_dataset(data_dir_era5 + v + '.nc')

        # Cases
        for c, c_count in zip(cases, range(0, len(cases))):
            snr, case_num = CaseSNR(data, s, min_max_months[s_count],
                                    c, nc_date_dir)

            if save_nc & (s=='SON') & (hpalevel==200):
                snr.to_netcdf(out_dir2 + v.split('_')[0] + s + '_' +
                              cases[c_count] + '_mer_d_w_' + aux_name + '.nc')
            if plot:
                PlotComposite_wWAF(comp=snr, cmap=cbar_snr_hgt, step1=1,
                                   levels=
                                   [-1, -.8, -.6, -.5, -.1, 0, 0.1, 0.5, 0.6,
                                    0.8, 1],
                                   contour1=True, two_variables=False,
                                   mapa='hs', significance=False,

                                   title='Signal-to-Noise ratio ' +
                                         v.split('_')[0] + ' - ' +
                                         title_case[c_count] +
                                         '\n' + s + ' - Events: ' +
                                         str(case_num),

                                   name_fig=v.split('_')[0] + s + '_' +
                                            cases[c_count] + '_mer_d_w_' +
                                            aux_name,

                                   out_dir=out_dir, dpi=dpi, save=save,
                                   color_map='grey')

print('Done HGT')
################################################################################
print('PP y T')
# PP y T
variables_tpp = ['ppgpcc_w_c_d_1', 'tcru_w_c_d_0.25']
title_var = ['PP GPCC', 'T Cru']

# seasons
for s, s_count in zip(seasons, range(0, len(seasons))):
    # variable
    for v, v_count in zip(variables_tpp, [0,1]):

        data = OpenObsDataSet(name=v + '_' + s, sa=True)
        data = data.sel(time=slice('1940-01-01', '2020-12-31'))

        # Cases
        for c, c_count in zip(cases, range(0, len(cases))):
            snr, case_num = CaseSNR(data, s, min_max_months[s_count],
                                    c, nc_date_dir)

            snr = snr*MakeMask(snr, 'var')

            PlotComposite_wWAF(comp=snr, cmap=cbars[v_count], step1=1,
                               levels =
                               [-1,-.8,-.6,-.5,-.1,0,0.1,0.5,0.6,0.8,1],
                               contour1=True, two_variables=False,
                               mapa='sa', significance=False,

                               title='Signal-to-Noise ratio ' +
                                     title_var[v_count] + ' - ' +
                                     title_case[c_count] +
                                     '\n' + s + ' - Events: ' +
                                     str(case_num),

                               name_fig=v.split('_')[0] + s + '_' +
                                        cases[c_count] + '_mer_d_w_' +
                                        aux_name,

                               out_dir=out_dir, dpi=dpi, save=save,
                               color_map='grey')
print('Done PP y T')

################################################################################
print('-------------------------------------------------')
print('done')
################################################################################
















