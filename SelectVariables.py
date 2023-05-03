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

# variables = ['sst']
cases = ['dmi_puros_pos', 'dmi_puros_neg', #DMI puros
         'n34_puros_pos', 'n34_puros_neg', #N34 puros
         'sim_pos', 'sim_neg', #sim misma fase
         'neutros', #neutros
         'dmi_pos', 'dmi_neg', 'n34_pos', 'n34_neg'] #todos de cada caso para validación

########################################################################################################################
def SelectEventsHGT(c):
    s = 'SON' #for s in ['JJA','SON']:
    try:
        aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
            .rename({'__xarray_dataarray_variable__': 'index'})
    except:
        aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
            .rename({'sst': 'index'})

    data_hgt_s = xr.open_dataset(data_dir + 'hgt_' + s.lower() + '.nc')

    case_events = SelectVariables(aux_cases, data_hgt_s)
    case_events.to_netcdf(out_dir + 'hgt_' + c + '_' + s + '.nc')

# Multiprocess #########################################################################################################
from multiprocessing import Process
processes = [Process(target=SelectEventsHGT, args=(c,)) for c in cases]
for process in processes:
    process.start()