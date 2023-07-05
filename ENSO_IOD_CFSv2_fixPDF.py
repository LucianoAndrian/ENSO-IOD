"""
PDF de las anomalías de PP (y T cuando ande...) en distans regiones de
SA para cada caso de IOD, ENSO y ENSO-IOD
"""
################################################################################
import xarray as xr
from matplotlib import colors
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
warnings.filterwarnings('ignore')
import warnings
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
from sklearn.cluster import KMeans
################################################################################
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/PDFs/'

save = False
dpi = 100
# Funciones ####################################################################
def best_fit_distribution(data, size, start, end):
    """Model data by finding best fit distribution to data"""
    # Get histogram of original data
    #y, x = np.histogram(data, bins=bins, density=True, weights=np.ones(len(data))/len(data)*100)
    y, x = np.histogram(data)
    x = (x + np.roll(x, -1))[:-1] / 2.0

    # Best holders
    best_distributions = []

    # Estimate distribution parameters from data
    # for ii, distribution_name in enumerate([d for d in _distn_names if not d in ['arcsine', 'beta', 'laplace', 'rdist', 'reciprocal' ,
    #                                                                         'halfcauchy' ,'bradford', 'genpareto' ,'gennorm',
    #                                                                         'genexpon', 'dgamma', 'foldcauchy', 'expon' ,'dweibull',
    #                                                                         'cauchy' ,'levy_stable', 'studentized_range', 'wrapcauchy',
    #                                                                         'wald', 'uniform', 'tukeylambda', 'truncnorm', 'truncexpon',
    #                                                                         'triang', 'semicircular', 'reciprocal', 'rdist', 'powerlaw',
    #                                                                         'pareto', 'lomax', 'lognorm', 'loglaplace', 'levy_l','levy',
    #                                                                              'gausshyper', 'foldnorm']]):

                                                                            # algunas distribuciones que no se tienen en cuenta
    for distribution_name in ['norm']:
        distribution = getattr(st, distribution_name)

        # Try to fit the distribution
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                # fit dist to data
                params = distribution.fit(data)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]

                # Calculate fitted PDF and error with fit in distribution
                pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                sse = np.sum(np.power(y - pdf, 2.0))

                # identify if this distribution is better
                best_distributions.append((distribution, params, sse, distribution_name))

        except Exception:
            pass

    best_dist = sorted(best_distributions, key=lambda x: x[2])[0]
    dist = best_dist[0]
    params =  best_dist[1]

    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)

    return pdf, best_dist[-1]

def MakeMask(DataArray, dataname='mask'):
    import regionmask
    mask = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.mask(DataArray)
    mask = xr.where(np.isnan(mask), mask, 1)
    mask = mask.to_dataset(name=dataname)
    return mask

def CFSv2Data(v):
    out_data_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
    dir_hc = '/pikachu/datos/osman/nmme/monthly/hindcast/'
    try:
        data_cfsv2 = xr.open_dataset(out_data_dir + v + '_data_cfsv2_noRoll.nc')
    except:
        def fix_calendar(ds, timevar='time'):
            """
            agrega los dias a los archivos nc de NMME
            """
            if ds[timevar].attrs['calendar'] == '360':
                ds[timevar].attrs['calendar'] = '360_day'
            return ds

        from ENSO_IOD_Funciones import SelectNMMEFiles
        files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                                dir=dir_hc, All=True)
        files = sorted(files, key=lambda x: x.split()[0])

        # abriendo todos los archivos
        data_cfsv2 = xr.open_mfdataset(files, decode_times=False).sel(
            L=[0.5, 1.5, 2.5, 3.5])  # xr no entiende la codificacion de Leads, r y las fechas
        data_cfsv2 = data_cfsv2.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
        data_cfsv2['L'] = [0, 1, 2, 3]
        data_cfsv2 = xr.decode_cf(fix_calendar(data_cfsv2))  # corrigiendo fechas
        data_cfsv2 = data_cfsv2.sel(lon=slice(275, 330), lat=slice(-60, 15))
        data_cfsv2 = data_cfsv2.sel(r=slice(1, 24))
        data_cfsv2.to_netcdf(out_data_dir + v + '_data_cfsv2_noRoll.nc')
        data_cfsv2 = data_cfsv2.compute()
    return data_cfsv2
################################################################################
seasons = ['JJA', 'SON']

cases = ['dmi_puros_pos', 'dmi_puros_neg',
        'n34_puros_pos', 'n34_puros_neg',
        'sim_pos', 'sim_neg',
        'dmi_neg_n34_pos', 'dmi_pos_n34_neg']

positive_cases = ['dmi_puros_pos', 'n34_puros_pos', 'sim_pos']
negative_cases = ['dmi_puros_neg', 'n34_puros_neg', 'sim_neg']

title_case = ['DMI pure - positive',
              'DMI pure - negative',
              'El Niño pure', 'La Niña pure',
              'DMI positive - El Niño',
              'DMI negative - La Niña']
# cajas ------------------------------------------------------------------------
box_name = ['S_SESA', 'N-SESA', 'SESA']
box_lats = [[-39,-29], [-29,-17], [-39,-17]]
box_lons = [[296, 315], [296, 315], [296,315]]
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
v = 'prec' # <<<<<<solo precpitacion>>>>>>
data_cfsv2 = CFSv2Data(v)

fix_factor = 30
units = '[mm]'
startend = -60
ylim = 0.020

#ylim = [[0.020,0.06], [0.7,0.7]]

for s in seasons:
    for bl, bt, name in zip(box_lons, box_lats, box_name):
        # Climatologia
        neutro = xr.open_dataset(
            cases_dir + v + '_neutros_' + s + '.nc').rename(
            {v: 'var'})
        neutro *= fix_factor
        mask = MakeMask(neutro, 'var')
        neutro *= mask
        neutro = neutro.sel(lon=slice(bl[0], bl[1]),
                            lat=slice(bt[0], bt[1]))

        for c, c_count in zip(cases, range(0, len(cases))):
            case = xr.open_dataset(
                cases_dir + v + '_' + c + '_' + s + '.nc').rename({v: 'var'})
            case *= mask
            case *= fix_factor
            case = case.sel(lat=slice(bl[0], bl[1]),
                            lon=slice(bt[0], bt[1]))

            # creando la climatologia con TODOS los eventos
            # (esto es mas rapido que hacer otro SELECT_EVENTS
            # if c_count == 0:
            #     clim_full = xr.concat([neutro, case], dim='time')
            #     c_count = 1
            # else:
            #     clim_full = xr.concat([clim_full, case], dim='time')

        #clim_full = clim_full.mean(['lon', 'lat'])
        clim_full = neutro.mean(['lon', 'lat'])

        # plot
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.5), dpi=dpi)
        # ax[0].set_xlim(-startend, startend)
        # ax[1].set_xlim(-startend, startend)
        # ax[0].set_ylim(0, ylim)
        # ax[1].set_ylim(0, ylim)
        ax[0].set_xlabel(v + '- anomaly ' + units)
        ax[1].set_xlabel(v + '- anomaly ' + units)
        ax[0].set_title('Clim - Positive Phases', color='firebrick')
        ax[1].set_title('Clim. - Negative Phases', color='royalblue')

        cases_color = ['darkorange', 'blue', 'firebrick', 'cyan', 'forestgreen',
                       'purple']

        positive_cases_colors = ['red', '#F47B00', 'forestgreen']
        negative_cases_colors = ['blue', '#55F400', 'magenta']

        aux_clim_full = clim_full - clim_full.mean('time')
        aux_clim_full = aux_clim_full['var'].values

        pdf_clim_full, name_clim_full_dist = best_fit_distribution(
            np.nan_to_num(aux_clim_full), len(aux_clim_full),
            -1 * startend, startend)

        ax[0].plot(pdf_clim_full, lw=1.5, color='k', label='climatology')
        ax[1].plot(pdf_clim_full, lw=1.5, color='k', label='climatology')

        for c, c_count in zip(positive_cases, range(0, len(positive_cases))):
            case = xr.open_dataset(
                cases_dir + v + '_' + c + '_' + s + '.nc').rename({v: 'var'})
            case *= fix_factor
            case = case.sel(lon=slice(bl[0], bl[1]),
                            lat=slice(bt[0], bt[1]))
            case = case.mean(['lon', 'lat'])

            # signal
            case_anom = case - clim_full.mean('time')
            case_anom = case_anom['var'].values

            pdf_case, name_case_dist = best_fit_distribution(
                np.nan_to_num(case_anom), len(case_anom),
                -1 * startend, startend)
            ax[0].plot(pdf_case, lw=2.5, color=positive_cases_colors[c_count],
                       label=c)

        for c, c_count in zip(negative_cases, range(0, len(negative_cases))):
            case = xr.open_dataset(
                cases_dir + v + '_' + c + '_' + s + '.nc').rename({v: 'var'})
            case *= fix_factor
            case = case.sel(lon=slice(bl[0], bl[1]),
                            lat=slice(bt[0], bt[1]))
            case = case.mean(['lon', 'lat'])

            # signal
            case_anom = case - clim_full.mean('time')
            case_anom = case_anom['var'].values

            pdf_case, name_case_dist = best_fit_distribution(
                np.nan_to_num(case_anom), len(case_anom),
                -1 * startend, startend)
            ax[1].plot(pdf_case, lw=2.5, color=negative_cases_colors[c_count],
                       label=c)

        ax[0].grid(alpha=0.5)
        ax[1].grid(alpha=0.5)
        ax[0].legend(loc='best')
        ax[1].legend(loc='best')

        fig.suptitle(v + ' - ' + name + ' - ' + s)
        # fig.legend(bbox_to_anchor=(1, 0.5), loc="center right")  # , bbox_transform=fig.transFigure, ncol=len(cases2))
        # plt.subplots_adjust(right=0.85)
        plt.yticks(size=10)
        plt.xticks(size=10)
        # plt.tight_layout(rect=[0, 0, 0.7, 0])
        if save:
            plt.savefig(
                out_dir + v + '_PDF_' + 'cluster_name[cs]' + '_' + s + '.jpg',
                box_inches="tight")
        else:
            plt.show()