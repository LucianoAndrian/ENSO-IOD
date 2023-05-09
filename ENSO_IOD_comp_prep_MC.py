########################################################################################################################
"""
Igual que ENSO-IOD.py pero sólo prepara las fechas (años) de los eventos, sim, un y all. NO COMPOSITE
"""
########################################################################################################################
import xarray as xr
import pandas as pd
pd.options.mode.chained_assignment = None
import os
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
from ENSO_IOD_Funciones import Nino34CPC
from ENSO_IOD_Funciones import DMI
from ENSO_IOD_Funciones import MultipleComposite

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

w_dir = '/home/luciano.andrian/doc/salidas/'
out_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates/'
pwd = '/datos/luciano.andrian/ncfiles/'

########################################################################################################################
#
# seasons = [6, 7, 8, 9, 10]
# seasons_name = ['MJJ', 'JJA', 'JAS','ASO', 'SON']
#
# full_season2 = False
# bwa = False

seasons = [7, 10]
seasons_name = ['JJA', 'SON']

full_season2 = False
bwa = False

# Usando T: Beic y PP:GPCC se cubre 1920-2020
# Se puede usar junto los ERA-Frankenstein... o solo ERA5 o ERA5 + ERA5_50-78

#start = [1920,1950,1980] #*ERA5 va desde 1979 pero es una molestia en Nino34CPC y su climatologia movil.
start = 1940
end = 2020
########################################################################################################################

# for i in start:
#     for fs in full_season2:

i = 1920
fs = False
dmi, aux, dmi_aux = DMI(filter_bwa=False, start_per=i, end_per=end)
del aux
aux = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
n34 = Nino34CPC(aux, start=i, end=end)[2]

if fs:
    full_season = True
    s = 666
    print('Full Season')
    Neutral, DMI_sim_pos, DMI_sim_neg, DMI_un_pos, DMI_un_neg, N34_un_pos, N34_un_neg, DMI_pos, \
    DMI_neg, N34_pos, N34_neg, DMI_pos_N34_neg, DMI_neg_N34_pos = \
        MultipleComposite(var=dmi_aux, n34=n34, dmi=dmi, season=s, start=i, full_season=full_season,
                          compute_composite=False)

    ds = xr.Dataset(
        data_vars={
            'Neutral': (Neutral),
            "DMI_sim_pos": (DMI_sim_pos),
            "DMI_sim_neg": (DMI_sim_neg),
            "DMI_un_pos": (DMI_un_pos),
            "DMI_un_neg": (DMI_un_neg),
            "N34_un_pos": (N34_un_pos),
            "N34_un_neg": (N34_un_neg),
            "DMI_pos": (DMI_pos),
            "DMI_neg": (DMI_neg),
            "N34_pos": (N34_pos),
            "N34_neg": (N34_neg),
        }
    )
    ds.to_netcdf(pwd + 'nc_composites_dates/' + 'Composite_' +
                 str(i) + '_' + str(end) + '_Full_Season.nc')

else:
    full_season = False
    s_count=0
    for s in seasons:
        print(seasons_name[s_count])
        Neutral, DMI_sim_pos, DMI_sim_neg, DMI_un_pos, DMI_un_neg, N34_un_pos, N34_un_neg, DMI_pos, \
        DMI_neg, N34_pos, N34_neg, DMI_pos_N34_neg, DMI_neg_N34_pos =\
            MultipleComposite(var=dmi_aux, n34=n34, dmi=dmi, season=s - 1, start=i, full_season=full_season,
                              compute_composite=False)

        ds = xr.Dataset(
            data_vars={
                'Neutral': (Neutral),
                "DMI_sim_pos": (DMI_sim_pos),
                "DMI_sim_neg": (DMI_sim_neg),
                "DMI_un_pos": (DMI_un_pos),
                "DMI_un_neg": (DMI_un_neg),
                "N34_un_pos": (N34_un_pos),
                "N34_un_neg": (N34_un_neg),
                "DMI_pos": (DMI_pos),
                "DMI_neg": (DMI_neg),
                "N34_pos": (N34_pos),
                "N34_neg": (N34_neg),
            }
        )
        ds.to_netcdf(out_dir +
                     str(i) + '_' + str(end) + '_' + seasons_name[s_count] + '.nc')

        s_count += 1
########################################################################################################################
########################################################################################################################
# DMI sin tener en cuenta magnitud ni signo de las anomalias de sst en iode y iodw, solo dmi.
########################################################################################################################
out_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates_no_ind_sst_anom/'

i = 1920
fs = False
dmi, aux, dmi_aux = DMI(filter_bwa=False, start_per=i, end_per=end,
                        sst_anom_sd=False, opposite_signs_criteria=False)
del aux
aux = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
n34 = Nino34CPC(aux, start=i, end=end)[2]

if fs:
    full_season = True
    s = 666
    print('Full Season')
    Neutral, DMI_sim_pos, DMI_sim_neg, DMI_un_pos, DMI_un_neg, N34_un_pos, N34_un_neg, DMI_pos, \
    DMI_neg, N34_pos, N34_neg, DMI_pos_N34_neg, DMI_neg_N34_pos =\
        MultipleComposite(var=dmi_aux, n34=n34, dmi=dmi, season=s, start=i, full_season=full_season,
                          compute_composite=False)

    ds = xr.Dataset(
        data_vars={
            'Neutral': (Neutral),
            "DMI_sim_pos": (DMI_sim_pos),
            "DMI_sim_neg": (DMI_sim_neg),
            "DMI_un_pos": (DMI_un_pos),
            "DMI_un_neg": (DMI_un_neg),
            "N34_un_pos": (N34_un_pos),
            "N34_un_neg": (N34_un_neg),
            "DMI_pos": (DMI_pos),
            "DMI_neg": (DMI_neg),
            "N34_pos": (N34_pos),
            "N34_neg": (N34_neg),
        }
    )
    ds.to_netcdf(pwd + 'nc_composites_dates/' + 'Composite_' +
                 str(i) + '_' + str(end) + '_Full_Season.nc')

else:
    full_season = False
    s_count=0
    for s in seasons:
        print(seasons_name[s_count])
        Neutral, DMI_sim_pos, DMI_sim_neg, DMI_un_pos, DMI_un_neg, N34_un_pos, N34_un_neg, DMI_pos, \
        DMI_neg, N34_pos, N34_neg, DMI_pos_N34_neg, DMI_neg_N34_pos = \
            MultipleComposite(var=dmi_aux, n34=n34, dmi=dmi, season=s - 1, start=i, full_season=full_season,
                              compute_composite=False)

        ds = xr.Dataset(
            data_vars={
                'Neutral': (Neutral),
                "DMI_sim_pos": (DMI_sim_pos),
                "DMI_sim_neg": (DMI_sim_neg),
                "DMI_un_pos": (DMI_un_pos),
                "DMI_un_neg": (DMI_un_neg),
                "N34_un_pos": (N34_un_pos),
                "N34_un_neg": (N34_un_neg),
                "DMI_pos": (DMI_pos),
                "DMI_neg": (DMI_neg),
                "N34_pos": (N34_pos),
                "N34_neg": (N34_neg),
            }
        )
        ds.to_netcdf(out_dir +
                     str(i) + '_' + str(end) + '_' + seasons_name[s_count] + '.nc')

        s_count += 1
########################################################################################################################