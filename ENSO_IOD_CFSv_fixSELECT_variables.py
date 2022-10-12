"""
Seleccion de los campos de las variables para cada caso de eventos IOD y ENSO
En lugar de usar los datos de LoadLeads*.py, las fechas y los sst_*.nc de ENSO_IOD_fixCFSv2_DMI_N34.py para asegurar
correspondencia entre los eventos de los índices y los campos seleccionados. (esto fue uno de los problemas iniciales)

SST ~20min
HGT ~>2hs
PP ~20min
"""
########################################################################################################################
import xarray as xr
import numpy as np
from multiprocessing.pool import ThreadPool
from ENSO_IOD_Funciones import SelectVariables
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
########################################################################################################################
cases_date_dir = '/pikachu/datos/luciano.andrian/cases_dates/'
data_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
########################################################################################################################
# variables = ['sst']
seasons = ['JJA', 'JAS', 'SON', 'ASO']
cases = ['dmi_puros_pos', 'dmi_puros_neg', 'n34_puros_pos', 'n34_puros_neg', 'sim_pos', 'sim_neg',
         'neutros', 'dmi_neg_n34_pos', 'dmi_pos_n34_neg']

#cases = ['dmi_neg_n34_pos', 'dmi_pos_n34_neg']


# SST ##################################################################################################################
def SelectEvents(s):
    for c in cases:
        aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc')\
            .rename({'__xarray_dataarray_variable__': 'index'})

        data_sst_s = xr.open_dataset(data_dir + 'sst_' + s.lower() + '.nc' )

        case_events = SelectVariables(aux_cases, data_sst_s)

        case_events.to_netcdf(out_dir + c + '_' + s + '.nc')

# Multiprocess #########################################################################################################
pool = ThreadPool(4) # uno por season
pool.map(SelectEvents, [s for s in seasons])
########################################################################################################################

# HGT ##################################################################################################################
def SelectEventsHGT(s):
    for c in cases:
        try:
            aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
                .rename({'__xarray_dataarray_variable__': 'index'})
        except:
            aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
                .rename({'sst': 'index'})

        data_hgt_s = xr.open_dataset(data_dir + 'hgt_' + s.lower() + '.nc' )

        case_events = SelectVariables(aux_cases, data_hgt_s)

        case_events.to_netcdf(out_dir + 'hgt_' + c + '_' + s + '.nc')

# Multiprocess #########################################################################################################
pool = ThreadPool(4) # uno por season
pool.map(SelectEventsHGT, [s for s in seasons])
########################################################################################################################

# PP ###################################################################################################################
def SelectEventsPP(s):
    for c in cases:
        try:
            aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
                .rename({'__xarray_dataarray_variable__': 'index'})
        except:
            aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
                .rename({'sst': 'index'})

        data_pp_s = xr.open_dataset(data_dir + 'prec_' + s.lower() + '.nc' )

        case_events = SelectVariables(aux_cases, data_pp_s)

        case_events.to_netcdf(out_dir + 'prec_' + c + '_' + s + '.nc')

# Multiprocess #########################################################################################################
pool = ThreadPool(4) # uno por season
pool.map(SelectEventsPP, [s for s in seasons])
########################################################################################################################

# T ####################################################################################################################
def SelectEventsT(s):
    for c in cases:
        try:
            aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
                .rename({'__xarray_dataarray_variable__': 'index'})
        except:
            aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
                .rename({'sst': 'index'})

        data_tref_s = xr.open_dataset(data_dir + 'tref_' + s.lower() + '.nc' )

        case_events = SelectVariables(aux_cases, data_tref_s)

        case_events.to_netcdf(out_dir + 'tref_' + c + '_' + s + '.nc')

# Multiprocess #########################################################################################################
# por alguna razon esto no anda...
# el uso momentaneo de casi todos los nucleos lo detiene? --> probar en otro momento
pool = ThreadPool(4) # uno por season
pool.map(SelectEventsT, [s for s in seasons])

for s in seasons:
    SelectEventsT(s)
########################################################################################################################

# PP - NO DETREND ######################################################################################################
def SelectEventsPP(s):
    for c in cases:
        aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc')\
            .rename({'__xarray_dataarray_variable__': 'index'})

        data_prec_s = xr.open_dataset(data_dir + 'prec_' + s.lower() + '_nodetrend.nc' )

        case_events = SelectVariables(aux_cases, data_prec_s)

        case_events.to_netcdf(out_dir + 'prec_' + c + '_' + s + '_nodetrend.nc')

# Multiprocess #########################################################################################################
pool = ThreadPool(4) # uno por season
pool.map(SelectEventsPP, [s for s in seasons])
########################################################################################################################

# DMI ##################################################################################################################
cases_data_dir_dmi = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/'
def SelectEventsDMI(s):
    for c in cases:
        try:
            aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
                .rename({'__xarray_dataarray_variable__': 'index'})
        except:
            aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
                .rename({'sst': 'index'})

        data_dmi_s = xr.open_dataset(cases_data_dir_dmi + 'DMI_' + s + '_Leads_r_CFSv2.nc')

        case_events = SelectVariables(aux_cases, data_dmi_s)

        case_events.to_netcdf(out_dir + 'dmi_values_' + c + '_' + s + '.nc')

# Multiprocess #########################################################################################################
pool = ThreadPool(4) # uno por season
pool.map(SelectEventsDMI, [s for s in seasons])
########################################################################################################################

# N34 ##################################################################################################################
cases_data_dir_dmi = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/'
def SelectEventsN34(s):
    for c in cases:
        try:
            aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
                .rename({'__xarray_dataarray_variable__': 'index'})
        except:
            aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
                .rename({'sst': 'index'})

        data_n34_s = xr.open_dataset(cases_data_dir_dmi + 'N34_' + s + '_Leads_r_CFSv2.nc')

        case_events = SelectVariables(aux_cases, data_n34_s)

        case_events.to_netcdf(out_dir + 'N34_values_' + c + '_' + s + '.nc')

# Multiprocess #########################################################################################################
pool = ThreadPool(4) # uno por season
pool.map(SelectEventsN34, [s for s in seasons])
########################################################################################################################