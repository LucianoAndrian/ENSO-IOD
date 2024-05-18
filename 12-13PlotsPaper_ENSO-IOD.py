"""
Figuras 12-13 paper ENSO-IOD
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
dates_dir = '/pikachu/datos/luciano.andrian/DMI_N34_Leads_r/'

if save:
    dpi = 300
else:
    dpi = 100

################################################################################
cases_magnitude = [None, None, 's_en', 's_en-m_iodp', 's_en-s_iodp',
                   None, None, 'm_en', 'm_en-m_iodp', 'm_en-s_iodp',
                   's_iodn', 'm_iodn', 'clima', 'm_iodp', 's_iodp',
                   'm_ln-s_iodn', 'm_ln-m_iodn', 'm_ln', None, None,
                   's_ln-s_iodn', 's_ln-m_iodn', 's_ln', None, None]

cases = ['dmi_puros_pos', 'dmi_puros_neg',
        'n34_puros_pos', 'n34_puros_neg',
        'sim_pos', 'sim_neg',
        'dmi_neg_n34_pos', 'dmi_pos_n34_neg',
        'neutros']



bin_limits = [[-4.5,-1], [-1, -0.5], #1
              [-0.5, 0.5], #2
              [0.5, 1], [1, 4.5]] #4


bins_by_cases_dmi = [[3, 4], [0, 1],
                     [2], [2],
                     [3, 4], [0, 1],
                     [0, 1],[3, 4],
                     [2]]

bins_by_cases_n34 = [[2], [2],
                     [3, 4], [0,1],
                     [3, 4], [0, 1],
                     [3, 4], [0, 1],
                     [2]]
from ENSO_IOD_Funciones import BinsByCases


bin_names = ['s', 'm', '', 'm', 's']
cases_names = []
for c_count, c  in enumerate(cases):
    aux_h = '-'
    for d in bins_by_cases_dmi[c_count]:
        d_aux = sum(bin_limits[d])
        d_aux_mag_name = bin_names[d]
        d_aux_h = '_'
        if d_aux>0:
            d_aux_name = 'iodp'
        elif d_aux<0:
            d_aux_name = 'iodn'
        elif d_aux==0:
            d_aux_name = ''
            d_aux_mag_name = ''
            d_aux_h = ''
            aux_h = ''

        iod_name = f"{d_aux_mag_name}{d_aux_h}{d_aux_name}"
        #print(iod_name)

        for n in bins_by_cases_n34[c_count]:
            n_aux = sum(bin_limits[n])
            n_aux_mag_name = bin_names[n]
            n_aux_h = '_'

            if n_aux > 0:
                n_aux_name = 'en'
            elif n_aux < 0:
                n_aux_name = 'ln'
            elif n_aux == 0:
                n_aux_name = ''
                n_aux_mag_name=''
                n_aux_h = ''
                aux_h = ''

            enso_name = f"{n_aux_mag_name}{n_aux_h}{n_aux_name}"

            case_name = f"{enso_name}{aux_h}{iod_name}"

            cases_names.append(case_name)


aux_comps = {}
n_count = 0
for c_count, c  in enumerate(cases):
    cases_bin, num_bin, aux = BinsByCases(v='hgt', v_name='hgt', fix_factor=9.8,
                                          s='SON', mm=10, c=c,
                                          c_count=c_count,
                                          bin_limits=bin_limits,
                                          bins_by_cases_dmi=bins_by_cases_dmi,
                                          bins_by_cases_n34=bins_by_cases_n34,
                                          snr=False,
                                          cases_dir=cases_dir,
                                          dates_dir=dates_dir,
                                          obsdates=False)

    bins_aux_dmi = bins_by_cases_dmi[c_count]
    bins_aux_n34 = bins_by_cases_n34[c_count]

    for b_dmi in range(0, len(bins_aux_dmi)):
        for b_n34 in range(0, len(bins_aux_n34)):
            aux_comps[cases_names[n_count]] = cases_bin[b_dmi][b_n34]

            n_count += 1