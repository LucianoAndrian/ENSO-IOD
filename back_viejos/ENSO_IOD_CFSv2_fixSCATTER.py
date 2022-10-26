"""
Scatter por seasons para ENSO vs IOD
"""
########################################################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
#from ENSO_IOD_Funciones import SelectVariables
#######################################################################################################################
dates_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/'
cases_date_dir = '/pikachu/datos/luciano.andrian/cases_dates/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Scatter/'

save = True
dpi=400
########################################################################################################################
indices = ['DMI', 'N34']
seasons = ['JJA', 'JAS', 'ASO', 'SON']
cases = ['puros_pos', 'puros_neg', 'sim_pos', 'sim_neg', 'neutros', 'pos', 'neg']
colors = ['firebrick', 'lime', 'red', 'lightseagreen', 'gray']
#----------------------------------------------------------------------------------------------------------------------#
for s in seasons:
    data_dates_dmi_or = xr.open_dataset(dates_dir + 'DMI_' + s + '_Leads_r_CFSv2.nc')
    data_dates_dmi = data_dates_dmi_or.where(np.abs(data_dates_dmi_or) > 0.5 * data_dates_dmi_or.mean('r').std())
    data_dates_dmi /= data_dates_dmi_or.std(['time', 'r'])
    data_dates_dmi_or /= data_dates_dmi_or.std(['time', 'r'])

    data_dates_n34_or = xr.open_dataset(dates_dir + 'N34_' + s + '_Leads_r_CFSv2.nc')
    data_dates_n34 = data_dates_n34_or.where(np.abs(data_dates_n34_or) > 0.5 * data_dates_n34_or.mean('r').std())
    data_dates_n34 /= data_dates_n34_or.std(['time', 'r'])
    data_dates_n34_or /= data_dates_n34_or.std(['time', 'r'])

    # que elegancia la de francia... ###################################################################################

    # dmi_puros_pos ---------------------------------------------------------------------------------------------------#
    aux_cases = xr.open_dataset(cases_date_dir + 'dmi_puros_pos_f_' + s + '.nc') \
        .rename({'__xarray_dataarray_variable__': 'index'})

    dmi_puros_pos = SelectVariables(aux_cases, data_dates_dmi)
    dmi_puros_pos_n34 = SelectVariables(aux_cases, data_dates_n34_or)

    # dmi_puros_neg ---------------------------------------------------------------------------------------------------#
    aux_cases = xr.open_dataset(cases_date_dir + 'dmi_puros_neg_f_' + s + '.nc') \
        .rename({'__xarray_dataarray_variable__': 'index'})

    dmi_puros_neg = SelectVariables(aux_cases, data_dates_dmi)
    dmi_puros_neg_n34 = SelectVariables(aux_cases, data_dates_n34_or)

    # n34_puros_pos ---------------------------------------------------------------------------------------------------#
    aux_cases = xr.open_dataset(cases_date_dir + 'n34_puros_pos_f_' + s + '.nc') \
        .rename({'__xarray_dataarray_variable__': 'index'})

    n34_puros_pos = SelectVariables(aux_cases, data_dates_n34)
    n34_puros_pos_dmi = SelectVariables(aux_cases, data_dates_dmi_or)

    # n34_puros_neg ---------------------------------------------------------------------------------------------------#
    aux_cases = xr.open_dataset(cases_date_dir + 'n34_puros_neg_f_' + s + '.nc') \
        .rename({'__xarray_dataarray_variable__': 'index'})

    n34_puros_neg = SelectVariables(aux_cases, data_dates_n34)
    n34_puros_neg_dmi = SelectVariables(aux_cases, data_dates_dmi_or)

    # sim_pos ---------------------------------------------------------------------------------------------------------#
    aux_cases = xr.open_dataset(cases_date_dir + 'sim_pos_f_' + s + '.nc') \
        .rename({'__xarray_dataarray_variable__': 'index'})

    sim_pos_dmi = SelectVariables(aux_cases, data_dates_dmi)
    sim_pos_n34 = SelectVariables(aux_cases, data_dates_n34)

    # sim_neg ---------------------------------------------------------------------------------------------------------#
    aux_cases = xr.open_dataset(cases_date_dir + 'sim_neg_f_' + s + '.nc') \
        .rename({'__xarray_dataarray_variable__': 'index'})

    sim_neg_dmi = SelectVariables(aux_cases, data_dates_dmi)
    sim_neg_n34 = SelectVariables(aux_cases, data_dates_n34)

    # neutros ---------------------------------------------------------------------------------------------------------#
    aux_cases = xr.open_dataset(cases_date_dir + 'neutros_f_' + s + '.nc') \
        .rename({'__xarray_dataarray_variable__': 'index'})

    neutro_dmi = SelectVariables(aux_cases, data_dates_dmi_or)
    neutro_n34 = SelectVariables(aux_cases, data_dates_n34_or)

    # dmi_neg_n34_pos_ ------------------------------------------------------------------------------------------------#
    aux_cases = xr.open_dataset(cases_date_dir + 'dmi_neg_n34_pos_f_' + s + '.nc') \
        .rename({'time': 'index'})  # por algun motivo esto tiene los valores de algun indice que no deberia estar

    dmi_neg_n34_pos_dmi = SelectVariables(aux_cases, data_dates_dmi)
    dmi_neg_n34_pos_n34 = SelectVariables(aux_cases, data_dates_n34)

    # dmi_pos_n34_neg -------------------------------------------------------------------------------------------------#
    aux_cases = xr.open_dataset(cases_date_dir + 'dmi_pos_n34_neg_f_' + s + '.nc') \
        .rename({'time': 'index'})

    dmi_pos_n34_neg_dmi = SelectVariables(aux_cases, data_dates_dmi)
    dmi_pos_n34_neg_n34 = SelectVariables(aux_cases, data_dates_n34)

    # PLOT #############################################################################################################

    fig, ax = plt.subplots(dpi=dpi)
    plt.scatter(y=dmi_puros_pos_n34.sst.values, x=dmi_puros_pos.sst.values, marker='>',
                s=20, edgecolor='k', color='firebrick', alpha=.5)
    plt.scatter(y=dmi_puros_neg_n34.sst.values, x=dmi_puros_neg.sst.values, marker='<',
                s=20, edgecolor='k', color='lime', alpha=.5)

    plt.scatter(y=n34_puros_pos.sst.values, x=n34_puros_pos_dmi.sst.values, marker='^',
                s=20, edgecolor='k', color='darkorange', alpha=.5)
    plt.scatter(y=n34_puros_neg.sst.values, x=n34_puros_neg_dmi.sst.values, marker='v',
                s=20, edgecolor='k', color='blue', alpha=.5)

    plt.scatter(y=sim_pos_n34.sst.values, x=sim_pos_dmi.sst.values, marker='o',
                s=20, edgecolor='k', color='red', alpha=.5)
    plt.scatter(y=sim_neg_n34.sst.values, x=sim_neg_dmi.sst.values, marker='o',
                s=20, edgecolor='k', color='lightseagreen', alpha=.5)

    plt.scatter(y=neutro_n34.sst.values, x=neutro_dmi.sst.values, marker='o',
                s=20, edgecolor='k', color='gray', alpha=.5)

    plt.scatter(y=dmi_pos_n34_neg_n34.sst.values, x=dmi_pos_n34_neg_dmi.sst.values, marker='o',
                s=20, edgecolor='k', color='purple', alpha=.5)
    plt.scatter(y=dmi_neg_n34_pos_n34.sst.values, x=dmi_neg_n34_pos_dmi.sst.values, marker='o',
                s=20, edgecolor='k', color='orange', alpha=.5)

    plt.ylim((-4, 4));
    plt.xlim((-4, 4))
    plt.axhspan(-.5, .5, alpha=0.2, color='black', zorder=0)
    plt.axvspan(-.5, .5, alpha=0.2, color='black', zorder=0)
    # ax.grid(True)
    fig.set_size_inches(6, 6)
    plt.xlabel('IOD', size=15)
    plt.ylabel('Niño 3.4', size=15)

    plt.text(-3.8, 3.6, 'EN/IOD-', dict(size=15))
    plt.text(-.3, 3.6, 'EN', dict(size=15))
    plt.text(+2.3, 3.6, 'EN/IOD+', dict(size=15))
    plt.text(+3, -.1, 'IOD+', dict(size=15))
    plt.text(+2.6, -3.9, 'LN/IOD+', dict(size=15))
    plt.text(-.3, -3.9, 'LN', dict(size=15))
    plt.text(-3.8, -3.9, ' LN/IOD-', dict(size=15))
    plt.text(-3.8, -.1, 'IOD-', dict(size=15))
    plt.title(s)
    plt.tight_layout()
    if save:
        plt.savefig(out_dir + 'ENSO_IOD_CFSv2_Scatter_' + s + '.jpg')
    else:
        plt.show()
########################################################################################################################